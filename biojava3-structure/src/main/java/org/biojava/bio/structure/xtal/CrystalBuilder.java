package org.biojava.bio.structure.xtal;


import java.util.ArrayList;
import java.util.Iterator;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3i;
import javax.vecmath.Vector3d;

import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.PDBCrystallographicInfo;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.contact.AtomContactSet;
import org.biojava.bio.structure.contact.StructureInterface;
import org.biojava.bio.structure.contact.StructureInterfaceList;



/**
 * A class containing methods to find interfaces in a given crystallographic Structure by
 * reconstructing the crystal lattice through application of symmetry operators
 *  
 * @author duarte_j
 *
 */
public class CrystalBuilder {
	
	// Default number of cell neighbors to try in interface search (in 3 directions of space). 
	// In the search, only bounding box overlaps are tried, thus there's not so much overhead in adding 
	// more cells. We actually tested it and using numCells from 1 to 10 didn't change runtimes at all.
	// Examples with interfaces in distant neighbor cells:
	//   2nd neighbors: 3hz3, 1wqj, 2de3, 1jcd
	//   3rd neighbors: 3bd3, 1men, 2gkp, 1wui
	//   5th neighbors: 2ahf, 2h2z
	//   6th neighbors: 1was (in fact interfaces appear only at 5th neighbors for it) 
	// Maybe this could be avoided by previously translating the given molecule to the first cell,
	// BUT! some bona fide cases exist, e.g. 2d3e: it is properly placed at the origin but the molecule 
	// is enormously long in comparison with the dimensions of the unit cell, some interfaces come at the 7th neighbor.
	// After a scan of the whole PDB (Oct 2013) using numCells=50, the highest one was 4jgc with 
	// interfaces up to the 11th neighbor. Other high ones (9th neighbors) are 4jbm and 4k3t.
	// We set the default value to 12 based on that (having not seen any difference in runtime)
	private static final int DEF_NUM_CELLS = 12;
	
	/**
	 * Whether to consider HETATOMs in contact calculations
	 */
	private static final boolean INCLUDE_HETATOMS = true;
	
	public static final Matrix4d IDENTITY = new Matrix4d(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1);

	private Structure structure;
	private PDBCrystallographicInfo crystallographicInfo;
	private int numChainsAu;
	private int numOperatorsSg;
	private int numOperatorsNcs;
	
	private Matrix4d[] operators;
	
	private boolean verbose;
	
	private int numCells;
	
	private ArrayList<Matrix4d> visited;
	

	
	public CrystalBuilder(Structure structure) {
		this.structure = structure;
		this.crystallographicInfo = structure.getCrystallographicInfo();
		
		this.numChainsAu = structure.getChains().size();		
		this.numOperatorsSg = 1;
		if (structure.isCrystallographic()) {
			this.numOperatorsSg = this.crystallographicInfo.getSpaceGroup().getMultiplicity();
		}
		this.numOperatorsNcs = 1;
		this.verbose = false;
		this.numCells = DEF_NUM_CELLS;				
		
	}
	
	/**
	 * Set the verbose flag for verbose output of search algorithm to stdout
	 * @param verbose
	 */
	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}
	
	/**
	 * Set the number of neighboring crystal cells that will be used in the search for contacts 
	 * @param numCells
	 */
	public void setNumCells(int numCells) {
		this.numCells = numCells;
	}
	
	private void initialiseVisited() {
		visited = new ArrayList<Matrix4d>();
	}
	
	/**
	 * Returns the list of unique interfaces that the given Structure has upon 
	 * generation of all crystal symmetry mates. An interface is defined as any pair of chains 
	 * that contact, i.e. for which there is at least a pair of atoms (one from each chain) within 
	 * the given cutoff distance.
	 * @param cutoff the distance cutoff for 2 chains to be considered in contact
	 * @return
	 */
	public StructureInterfaceList getUniqueInterfaces(double cutoff) {	

		
		StructureInterfaceList set = new StructureInterfaceList();
		
		// initialising the visited ArrayList for keeping track of symmetry redundancy
		initialiseVisited();
		
		
		
		// the isCrystallographic() condition covers 3 cases:
		// a) entries with expMethod X-RAY/other diffraction and defined crystalCell (most usual case)
		// b) entries with expMethod null but defined crystalCell (e.g. PDB file with CRYST1 record but no expMethod annotation) 
		// c) entries with expMethod not X-RAY (e.g. NMR) and defined crystalCell (NMR entries do have a dummy CRYST1 record "1 1 1 90 90 90 P1")



		calcInterfacesCrystal(set, cutoff, structure.isCrystallographic());

		
		return set;
	}
	
	/**
	 * Calculate interfaces between original asymmetric unit and neighboring 
	 * whole unit cells, including the original full unit cell i.e. i=0,j=0,k=0 
	 * @param set
	 * @param cutoff
	 */
	private void calcInterfacesCrystal(StructureInterfaceList set, double cutoff, boolean isCrystallographic) {

		
		// initialising debugging vars
		long start = -1; 
		long end = -1;
		int trialCount = 0;
		int skippedRedundant = 0;
		int skippedAUsNoOverlap = 0;
		int skippedChainsNoOverlap = 0;
		int skippedSelfEquivalent = 0;

		

		// generate complete unit cell, by applying all SG operators (array is of size numOperatorsSg)
		Structure[] cell = getUnitCell();

		// The bounding boxes of all AUs of the unit cell		
		UnitCellBoundingBox bbGrid = new UnitCellBoundingBox(numOperatorsSg, numOperatorsNcs, numChainsAu);;
		// we calculate all the bounds of each of the asym units, those will then be reused and translated
		bbGrid.setAllBbs(cell, INCLUDE_HETATOMS);
		
		// if not crystallographic there's no search to do in other cells, only chains within "AU" will be checked
		if (!isCrystallographic) numCells = 0;
		
		if (verbose) {
			trialCount = 0;
			start= System.currentTimeMillis();
			int neighbors = (2*numCells+1)*(2*numCells+1)*(2*numCells+1)-1;
			int auTrials = (numChainsAu*(numChainsAu-1))/2;
			int trials = numChainsAu*numOperatorsSg*numChainsAu*neighbors;
			System.out.println("Chain clash trials within original AU: "+auTrials);
			System.out.println(
					"Chain clash trials between the original AU and the neighbouring "+neighbors+
					" whole unit cells ("+numCells+" neighbours)" +
					"(2x"+numChainsAu+"chains x "+numOperatorsSg+"AUs x "+neighbors+"cells) : "+trials);
			System.out.println("Total trials: "+(auTrials+trials));
		}


		for (int a=-numCells;a<=numCells;a++) {
			for (int b=-numCells;b<=numCells;b++) {
				for (int c=-numCells;c<=numCells;c++) {
					
					Point3i trans = new Point3i(a,b,c);
					Vector3d transOrth = new Vector3d(a,b,c);
					if (a!=0 || b!=0 || c!=0)
						// we avoid doing the transformation for 0,0,0 (in case it's not crystallographic) 
						this.crystallographicInfo.getCrystalCell().transfToOrthonormal(transOrth);
					UnitCellBoundingBox bbGridTrans = bbGrid.getTranslatedBbs(transOrth);

					for (int csIdx = 0; csIdx<numOperatorsSg; csIdx++) {
					for (int ncsIdx = 0; ncsIdx<numOperatorsNcs; ncsIdx++) {
					
					//for (int n=0;n<numOperatorsSg*numOperatorsNcs;n++) { 

						int n = csIdx + ncsIdx;
						// short-cut strategies
						// 1) we skip first of all if the bounding boxes of the AUs don't overlap
						if (!bbGrid.getAuBoundingBox(0).overlaps(bbGridTrans.getAuBoundingBox(n), cutoff)) {
							if (verbose) skippedAUsNoOverlap++;
							continue;
						}

						// 2) we check if we didn't already see its equivalent symmetry operator partner 													
						//CrystalTransform tt = new CrystalTransform(this.crystallographicInfo.getSpaceGroup(), n);
						//tt.translate(trans);
						Matrix4d tt = new Matrix4d(operators[n]);
						crystalTranslate(tt, transOrth);
						if (isRedundant(tt)) { 								
							if (verbose) {
								System.out.println("Skipping redundant transformation: "+transfString(csIdx, ncsIdx, a, b, ncsIdx));
								skippedRedundant++;								
							}
							continue;
						}
						addVisited(tt);
						
						
						boolean selfEquivalent = false;
						
						// now we copy and actually translate the AU if we saw it does overlap and the sym op was not redundant
						Structure jAsym = null;
						if (n==0 && a==0 && b==0 && c==0) {
							// special case: for original AU there's no need to clone and no need to translate 
							jAsym = cell[0];
						} else { 
							jAsym = cell[n].clone();
							Calc.translate(jAsym, transOrth);
						}
						

						// 3) an operator can be "self redundant" if it is the inverse of itself (involutory, e.g. all pure 2-folds with no translation)						
						if (isEquivalent(tt,tt)) { 
							if (verbose) 
								System.out.println("Transform "+transfString(csIdx, ncsIdx, a, b, ncsIdx)+
										" is equivalent to itself, will skip half of i-chains to j-chains comparisons");
							// in this case we can't skip the operator, but we can skip half of the matrix comparisons e.g. j>i
							// we set a flag and do that within the loop below
							selfEquivalent = true;
						}
						
						if (verbose) System.out.print(transfString(csIdx, ncsIdx, a, b, ncsIdx)+" ");
						
						// Now that we know that boxes overlap and operator is not redundant, we have to go to the details 
						int contactsFound = 0;
												
						int j = -1;
						for (Chain chainj:jAsym.getChains()) {
							j++;
							int i = -1;
							for (Chain chaini:structure.getChains()) { // we only have to compare the original asymmetric unit to every full cell around
								i++;
								if(selfEquivalent && (j>i)) {
									// in case of self equivalency of the operator we can safely skip half of the matrix
									skippedSelfEquivalent++;
									continue;
								}
								// special case of original AU, we don't compare a chain to itself
								if (n==0 && a==0 && b==0 && c==0 && i==j) continue;
								
								// before calculating the AtomContactSet we check for overlap, then we save putting atoms into the grid
								if (!bbGrid.getChainBoundingBox(0,i).overlaps(bbGridTrans.getChainBoundingBox(n,j),cutoff)) {
									if (verbose) {
										skippedChainsNoOverlap++;
										System.out.print(".");
									}
									continue;
								}
								if (verbose) trialCount++;

								// finally we've gone through all short-cuts and the 2 chains seem to be close enough:
								// we do the calculation of contacts
								StructureInterface interf = calcContacts(chaini, chainj, cutoff, csIdx, ncsIdx, trans);
								
								if (interf!=null) {
									contactsFound++;
									set.add(interf);
								}
							}
						}
						if (verbose) {
							if (a==0 && b==0 && c==0 && n==0) 
								System.out.println(" "+contactsFound+"("+(numChainsAu*(numChainsAu-1))/2+")");
							else if (selfEquivalent) 								
								System.out.println(" "+contactsFound+"("+(numChainsAu*(numChainsAu+1))/2+")");								
							else
								System.out.println(" "+contactsFound+"("+numChainsAu*numChainsAu+")");
						}
					}
				}
				}
			}
		}
		
		if (verbose) {
			end = System.currentTimeMillis();
			System.out.println("\n"+trialCount+" chain-chain clash trials done. Time "+(end-start)/1000+"s");
			System.out.println("  skipped (not overlapping AUs)       : "+skippedAUsNoOverlap);
			System.out.println("  skipped (not overlapping chains)    : "+skippedChainsNoOverlap);
			System.out.println("  skipped (sym redundant op pairs)    : "+skippedRedundant);
			System.out.println("  skipped (sym redundant self op)     : "+skippedSelfEquivalent);

			System.out.println("Found "+set.size()+" interfaces.");
		}
	}

	private StructureInterface calcContacts(Chain chaini, Chain chainj, double cutoff, int csIdx, int ncsIdx, Point3i trans) {
		
		// note that we don't consider hydrogens when calculating contacts
		AtomContactSet graph = StructureTools.getAtomsInContact(chaini, chainj, cutoff, INCLUDE_HETATOMS);
		
		if (graph.size()>0) {
			if (verbose) System.out.print("x");
			
			
			String chainIdi = chaini.getChainID();
			String chainIdj = chainj.getChainID();
			
			if (numOperatorsNcs>1)
				chainIdj = chainIdj+":"+ncsIdx;
			
			CrystalTransform transfi = new CrystalTransform(this.crystallographicInfo.getSpaceGroup());
			
			CrystalTransform transfj = new CrystalTransform(this.crystallographicInfo.getSpaceGroup(), csIdx);
			transfj.translate(trans);
			
			StructureInterface interf = new StructureInterface(
					StructureTools.getAllAtomArray(chaini), StructureTools.getAllAtomArray(chainj),
					chainIdi, chainIdj,
					graph,
					transfi, transfj);

			return interf;
			
		} else {
			if (verbose) System.out.print("o");
			return null;
		}		
	}
	
	private void addVisited(Matrix4d tt) {
		visited.add(tt);
	}
	
	/**
	 * Checks whether given transformId/translation is symmetry redundant 
	 * Two transformations are symmetry redundant if their matrices (4d) multiplication gives the identity, i.e.
	 * if one is the inverse of the other.
	 * @param tt
	 * @return
	 */
	private boolean isRedundant(Matrix4d tt) {
		
		Iterator<Matrix4d> it = visited.iterator();
		while (it.hasNext()) {
			Matrix4d v = it.next();
			
			if (isEquivalent(v, tt)) {
				//if (verbose) System.out.println("Skipping redundant transformation: "+tt+", equivalent to "+v);
				
				// there's only 1 possible equivalent partner for each visited matrix 
				// (since the equivalent is its inverse matrix and the inverse matrix is unique)
				// thus once the partner has been seen, we don't need to check it ever again
				it.remove();
				return true;
			}
			
		}
		
		return false;
	}
	
	/**
	 * Generates all symmetry-related objects from this asym unit and returns the whole
	 * unit cell (this asymmetric unit plus the symmetry-related objects). 
	 * @return
	 */
	private Structure[] getUnitCell() {
		
		Matrix4d[] ncsOps = null;
		
		if (this.crystallographicInfo.getNcsOperators()!=null &&
			this.crystallographicInfo.getNcsOperators().length>0) {
					
			ncsOps = this.crystallographicInfo.getNcsOperators();

			if (verbose) 
				System.out.println(ncsOps.length+" NCS operators found, adding them...");
			
			this.numOperatorsNcs = ncsOps.length;
		}

		Structure[] aus = new Structure[numOperatorsSg*numOperatorsNcs];
		aus[0] = structure;

		operators = new Matrix4d[numOperatorsSg*numOperatorsNcs];
		
		operators[0] = IDENTITY;
		
		// for non-crystallographic cases we can return here already
		if (numOperatorsSg*numOperatorsNcs==1) return aus;
		
		Matrix4d[] csOps = this.crystallographicInfo.getTransformationsOrthonormal();				
		
		if (ncsOps == null) {			
			// 1) normal case: there's no NCS operators
			// for 0 we already set it to original AU above
			for (int csIdx=1;csIdx<numOperatorsSg;csIdx++) {
				operators[csIdx] = csOps[csIdx];				
				Structure sym = structure.clone();
				Calc.transform(sym, csOps[csIdx]);
				aus[csIdx] = sym;
			}
		} else {
			// 2) we have NCS operators, normally viral capsids
			for (int csIdx=0; csIdx<numOperatorsSg; csIdx++) {
				for (int ncsIdx=0; ncsIdx<numOperatorsNcs; ncsIdx++) {					
					Matrix4d composedOp = new Matrix4d();
					// note they are not conmutative: order matters
					// we apply first NCS operator and then the CS operator
					composedOp.mul(csOps[csIdx], ncsOps[ncsIdx]);					
					if ((csIdx+ncsIdx)!=0) {
						// for 0 (corresponding to composition of identity with identity) 
						// we already set it to original AU above						
						operators[csIdx+ncsIdx] = composedOp;
						Structure sym = structure.clone();
						Calc.transform(sym, composedOp);
						aus[csIdx+ncsIdx] = sym;
					}
				}

			}
		}

		return aus;
	}
	
	private static void crystalTranslate (Matrix4d m, Vector3d translation) {

		m.m03 = m.m03+(double)translation.x;
		m.m13 = m.m13+(double)translation.y;
		m.m23 = m.m23+(double)translation.z;

	}

	/**
	 * Returns true if the given transformation matrices are equivalent.
	 * Two transformations are equivalent if one is the inverse of the other, i.e.
	 * their matrices multiplication is equal to the identity.
	 * @param m1
	 * @param m2
	 * @return
	 */
	public static boolean isEquivalent(Matrix4d m1, Matrix4d m2) {
		Matrix4d mul = new Matrix4d();
		mul.mul(m1, m2);

		if (mul.epsilonEquals(IDENTITY, 0.0001)) {
			return true;
		}
		return false;
	}
	
	private String transfString(int csIdx, int ncsIdx, int a, int b, int c) {
		if (numOperatorsNcs==1) {
			return String.format("[%2d-(%2d,%2d,%2d)]", csIdx,a,b,c);
		} else {
			return String.format("[%2d:%2d-(%2d,%2d,%2d)]", csIdx,ncsIdx,a,b,c);
		}
	}
}
