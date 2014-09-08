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
	
	public static final Matrix4d IDENTITY = new Matrix4d(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1);

	
	/**
	 * Whether to consider HETATOMs in contact calculations
	 */
	private static final boolean INCLUDE_HETATOMS = true;
	
	private Structure structure;
	private PDBCrystallographicInfo crystallographicInfo;
	private int numChainsAu;
	private int numOperatorsSg;
	private int numOperatorsNcs;
	//private int numChainsSuperAu;
	
	private Matrix4d[] ops;
	private Matrix4d[] ncsOps;
	
	private boolean verbose;
	
	private int numCells;
	
	private ArrayList<CrystalTransform> visited;
	
	private UnitCellBoundingBox bbGrid;
	
	private long start; 
	private long end;
	private int trialCount;
	private int skippedRedundant;
	private int skippedAUsNoOverlap;
	private int skippedChainsNoOverlap;
	private int skippedSelfEquivalent;

	
	public CrystalBuilder(Structure structure) {
		this.structure = structure;
		this.crystallographicInfo = structure.getCrystallographicInfo();
		
		this.numChainsAu = structure.getChains().size();		
		
		if (structure.isCrystallographic()) {
			this.ops = structure.getCrystallographicInfo().getTransformationsOrthonormal();
			this.numOperatorsSg = ops.length;
		} else {
			this.ops = new Matrix4d[1];
			this.ops[0] = IDENTITY;
			this.numOperatorsSg = 1;
		}
		
		this.ncsOps = this.crystallographicInfo.getNcsOperators();
		if (this.ncsOps==null) {
			ncsOps = new Matrix4d[1];
			ncsOps[0] = IDENTITY;
			this.numOperatorsNcs = 1;
		} else {
			this.numOperatorsNcs = this.ncsOps.length;
		}
		
		//this.numChainsSuperAu = numChainsAu * numOperatorsNcs;
		
		this.verbose = false;
		this.numCells = DEF_NUM_CELLS;
		
		// if not crystallographic there's no search to do in other cells, only chains within "AU" will be checked
		if (!structure.isCrystallographic()) numCells = 0;
				
		
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
		visited = new ArrayList<CrystalTransform>();
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

		if (verbose) {
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
		start = -1; 
		end = -1;
		trialCount = 0;
		skippedRedundant = 0;
		skippedAUsNoOverlap = 0;
		skippedChainsNoOverlap = 0;
		skippedSelfEquivalent = 0;

			
		// The bounding boxes of all AUs of the unit cell
		// we calculate all the bounds of each of the asym units, those will then be reused and translated
		bbGrid = new UnitCellBoundingBox(structure, ops, ncsOps, INCLUDE_HETATOMS);
		
		
		if (verbose) {
			trialCount = 0;
			start= System.currentTimeMillis();
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

					for (int n=0;n<numOperatorsSg;n++) { 

						// short-cut strategies
						// 1) we skip first of all if the bounding boxes of the AUs don't overlap
						if (!bbGrid.getAuBoundingBox(0).overlaps(bbGridTrans.getAuBoundingBox(n), cutoff)) {
							if (verbose) skippedAUsNoOverlap++;
							continue;
						}

						// 2) we check if we didn't already see its equivalent symmetry operator partner 													
						CrystalTransform tt = new CrystalTransform(this.crystallographicInfo.getSpaceGroup(), n);
						tt.translate(trans);
						if (isRedundant(tt)) { 								
							if (verbose) skippedRedundant++;								
							continue;
						}
						addVisited(tt);
						
						
						checkChains(set, cutoff, bbGridTrans, transOrth, tt, n, a, b, c); 
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
	
	private void checkChains(StructureInterfaceList set, double cutoff,
			UnitCellBoundingBox bbGridTrans, Vector3d transOrth, CrystalTransform tt,
			int n, int a, int b, int c) {
		
		boolean selfEquivalent = false;

		// 3) an operator can be "self redundant" if it is the inverse of itself (involutory, e.g. all pure 2-folds with no translation)						
		if (tt.isEquivalent(tt)) { 
			if (verbose) 
				System.out.println("Transform "+tt+" is equivalent to itself, will skip half of i-chains to j-chains comparisons");
			// in this case we can't skip the operator, but we can skip half of the matrix comparisons e.g. j>i
			// we set a flag and do that within the loop below
			selfEquivalent = true;
		}
		
		if (verbose) System.out.print(tt+" ");
		
		
		// Now that we know that boxes overlap and operator is not redundant, we have to go to the details 
		int contactsFound = 0;

		for (int iAuIdx=0;iAuIdx<numChainsAu;iAuIdx++) {
			for (int iNcsIdx=0;iNcsIdx<numOperatorsNcs;iNcsIdx++) {
				int i = iNcsIdx*numChainsAu + iAuIdx;
				for (int jAuIdx=0;jAuIdx<numChainsAu;jAuIdx++) {
					for (int jNcsIdx=0;jNcsIdx<numOperatorsNcs;jNcsIdx++) {
						int j = jNcsIdx*numChainsAu + jAuIdx ;


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
						//Chain chainj = null;

						// we only have to compare the original asymmetric unit to every full cell around
						Chain chaini = (Chain)structure.getChain(iAuIdx).clone();
						Matrix4d iNcsOp = ncsOps[iNcsIdx];
						Calc.transform(chaini,iNcsOp);

						//if (n==0 && a==0 && b==0 && c==0) {
						//	// special case of original AU
						//	chainj = structure.getChain(j);
						//} else {

						Chain chainj = (Chain)structure.getChain(jAuIdx).clone();
						Matrix4d jNcsOp = ncsOps[jNcsIdx];
						Matrix4d csOp = new Matrix4d(ops[n]);
						translate(csOp, transOrth);
						Matrix4d composed = new Matrix4d();
						composed.mul(csOp, jNcsOp);
						Calc.transform(chainj,composed);

						//}

						StructureInterface interf = calcContacts(chaini, iNcsIdx, chainj, jNcsIdx, cutoff, tt);

						if (interf!=null) {
							contactsFound++;
							set.add(interf);
						}
					}
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

	private StructureInterface calcContacts(Chain chaini, int iNcsIdx, Chain chainj, int jNcsIdx, double cutoff, CrystalTransform tt) {
		
		// note that we don't consider hydrogens when calculating contacts
		AtomContactSet graph = StructureTools.getAtomsInContact(chaini, chainj, cutoff, INCLUDE_HETATOMS);
		
		if (graph.size()>0) {
			if (verbose) System.out.print("x");
			
			String iChainId = chaini.getChainID();
			String jChainId = chainj.getChainID();
			if (iNcsIdx>0) iChainId += ":" + iNcsIdx;
			if (jNcsIdx>0) jChainId += ":" + jNcsIdx;
			
			CrystalTransform transf = new CrystalTransform(this.crystallographicInfo.getSpaceGroup());
			StructureInterface interf = new StructureInterface(
					StructureTools.getAllAtomArray(chaini), StructureTools.getAllAtomArray(chainj),
					iChainId, jChainId,
					graph,
					transf, tt);

			return interf;
			
		} else {
			if (verbose) System.out.print("o");
			return null;
		}		
	}
	
	private void addVisited(CrystalTransform tt) {
		visited.add(tt);
	}
	
	/**
	 * Checks whether given transformId/translation is symmetry redundant 
	 * Two transformations are symmetry redundant if their matrices (4d) multiplication gives the identity, i.e.
	 * if one is the inverse of the other.
	 * @param tt
	 * @return
	 */
	private boolean isRedundant(CrystalTransform tt) {
		
		Iterator<CrystalTransform> it = visited.iterator();
		while (it.hasNext()) {
			CrystalTransform v = it.next();
			
			if (tt.isEquivalent(v)) {

				if (verbose) System.out.println("Skipping redundant transformation: "+tt+", equivalent to "+v);
				
				// there's only 1 possible equivalent partner for each visited matrix 
				// (since the equivalent is its inverse matrix and the inverse matrix is unique)
				// thus once the partner has been seen, we don't need to check it ever again
				it.remove();
				
				return true;
			}
		}
		
		return false;
	}
	
	public void translate(Matrix4d m, Vector3d translation) {
		m.m03 = m.m03+(double)translation.x;
		m.m13 = m.m13+(double)translation.y;
		m.m23 = m.m23+(double)translation.z;
		
	}
	
//	/**
//	 * If NCS operators are given in MTRIX records, a bigger AU has to be constructed based on those.
//	 * Later they have to be removed with {@link #removeExtraChains()}
//	 */
//	private void constructFullStructure() {
//		
//		if (this.crystallographicInfo.getNcsOperators()==null ||
//			this.crystallographicInfo.getNcsOperators().length==0) {
//			// normal case: nothing to do			
//			return;
//		}
//				
//		// first we store the original chains in a new list to be able to restore the structure to its original state afterwards
//		origChains = new ArrayList<Chain>();
//		for (Chain chain:structure.getChains()) {
//			origChains.add(chain);
//		}
//		
//		// if we are here, it means that the NCS operators exist and we have to complete the given AU by applying them
//		Matrix4d[] ncsOps = this.crystallographicInfo.getNcsOperators();
//
//		if (verbose) 
//			System.out.println(ncsOps.length+" NCS operators found, generating new AU...");
//
//		
//		for (int i=0;i<ncsOps.length;i++) {
//			Structure transformedStruct = (Structure)structure.clone();			   
//			Calc.transform(transformedStruct, ncsOps[i]);
//			
//			for (Chain chain: transformedStruct.getChains()) {
//				// we assign a new AU id (chain ID) consisting in original chain ID + an operator index from 1 to n
//				chain.setChainID(chain.getChainID()+(i+1));
//				structure.addChain(chain);
//			}
//		}
//		
//		// now we have more chains in AU, we have to update the value 
//		this.numChainsAu = structure.getChains().size();
//	}
//	
//	/**
//	 * Removes the extra chains that were added to the original structure in {@link #constructFullStructure()}
//	 */
//	private void removeExtraChains() {
//		structure.setChains(origChains);
//	}
}
