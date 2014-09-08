package org.biojava.bio.structure.xtal;

import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;

import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.contact.BoundingBox;

/**
 * A class to contain the BoundingBoxes of all molecules in a full unit cell
 * 
 * @author duarte_j
 *
 */
public class UnitCellBoundingBox {

	/**
	 * An array with dimensions numOperatorsSg x numChainsSuperAu to contain all 
	 * bounding boxes of all chains of all AUs in unit cell
	 * e.g. chainBbs[0] would be the bounding boxes for all chains in the original AU
	 */
	private BoundingBox[][] chainBbs;
	
	/**
	 * An array with dimensions numOperatorsSg to contain all bounding boxes of
	 * all AUs in unit cell
	 */
	private BoundingBox[] auBbs;
	
	private int numOperatorsSg; // i.e. multiplicity of space group
	private int numOperatorsNcs;
	private int numChainsAu;
	private int numChainsSuperAu;
	
	public UnitCellBoundingBox(Structure structure, Matrix4d[] ops, Matrix4d[] ncsOps, boolean includeHetAtoms) {
		this.numOperatorsSg = ops.length;
		this.numChainsAu = structure.getChains().size();
		this.numOperatorsNcs = ncsOps.length;
		this.numChainsSuperAu = numOperatorsNcs * numChainsAu;
		
		initArrays();
		
		setBbs(structure, ops, ncsOps, includeHetAtoms);
	}
	
	private UnitCellBoundingBox(int numOperatorsSg, int numOperatorsNcs, int numChainsAu) {
		this.numChainsAu = numChainsAu;
		this.numOperatorsNcs = numOperatorsNcs;
		this.numOperatorsSg = numOperatorsSg;
		this.numChainsSuperAu = numOperatorsNcs * numChainsAu;
		initArrays();
	}
	
	private void initArrays() {
		this.chainBbs = new BoundingBox[numOperatorsSg][numChainsSuperAu];
		for (int cellIdx=0;cellIdx<numOperatorsSg;cellIdx++) {
			this.chainBbs[cellIdx] = new BoundingBox[numChainsSuperAu];
		}
		this.auBbs = new BoundingBox[numOperatorsSg];
	}
	
	private void setBbs(Structure structure, Matrix4d[] ops, Matrix4d[] ncsOps, boolean includeHetAtoms) {

		//setBb(structure, includeHetAtoms, 0);
		
		for (int csOpIdx=0;csOpIdx<numOperatorsSg;csOpIdx++) {
			for (int ncsOpIdx=0;ncsOpIdx<numOperatorsNcs;ncsOpIdx++) {
				Structure mate = structure.clone();	
				Matrix4d composed = new Matrix4d();
				composed.mul(ops[csOpIdx], ncsOps[ncsOpIdx]);
				Calc.transform(mate, ops[csOpIdx]); 
				
				for (int j = 0;j<numChainsAu; j++) {
					chainBbs[csOpIdx][ncsOpIdx*numChainsAu+j] = 
							new BoundingBox(StructureTools.getAllNonHAtomArray(mate.getChain(j), includeHetAtoms));
				}
			}
		}
		
		for (int csOpIdx=0;csOpIdx<numOperatorsSg;csOpIdx++) {	
			auBbs[csOpIdx] = new BoundingBox(chainBbs[csOpIdx]);
		}

	}
		
	/**
	 * Get the chain BoundingBox for the given cell index (cellIdx=0 would be original AU)
	 * and chain index
	 * @param cellIdx
	 * @param chainIdx
	 * @return
	 */
	public BoundingBox getChainBoundingBox(int cellIdx, int chainIdx) {
		return chainBbs[cellIdx][chainIdx];
	}
	
	/**
	 * Get the AU BoundingBox for the given cell index (cellIdx=0 would be original AU)
	 * The AU BoundingBox is the BoundingBox that bounds all chains belonging to the AU
	 * @param cellIdx
	 * @return
	 */
	public BoundingBox getAuBoundingBox(int cellIdx) {		
		return auBbs[cellIdx];		
	}
	
	/**
	 * Returns a new BoundingBoxes object containing the same bounds as this 
	 * BoundingBoxes object translated by the given translation
	 * @param translation
	 * @return
	 */
	public UnitCellBoundingBox getTranslatedBbs(Vector3d translation) {
		UnitCellBoundingBox translatedBbs = new UnitCellBoundingBox(numOperatorsSg, numOperatorsNcs, numChainsAu);
		
		for (int i=0; i<numOperatorsSg; i++) {
			for (int j = 0;j<numChainsSuperAu; j++) {
				translatedBbs.chainBbs[i][j] = new BoundingBox(this.chainBbs[i][j]);
				translatedBbs.chainBbs[i][j].translate(translation);
			}
			translatedBbs.auBbs[i] = new BoundingBox(translatedBbs.chainBbs[i]);
		}
		
		return translatedBbs;
	}
	
}
