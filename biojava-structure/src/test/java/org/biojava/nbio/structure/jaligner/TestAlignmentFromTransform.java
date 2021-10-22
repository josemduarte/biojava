package org.biojava.nbio.structure.jaligner;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.StructureAlignmentFactory;
import org.biojava.nbio.structure.align.fatcat.FatCatRigid;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.jama.Matrix;
import org.junit.Test;

import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;
import java.util.Arrays;
import java.util.List;

public class TestAlignmentFromTransform {

    private static final Matrix4d ID_MATRIX = new Matrix4d(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);

    @Test
    public void testAlignerAgainstFatcatAlignment() throws Exception {

        String pdbId1 = "4HGX";
        String asymId1 = "A";
        String pdbId2 = "3NGF";
        String asymId2 = "A";

        Structure s1 = StructureIO.getStructure(pdbId1);
        Structure s2 = StructureIO.getStructure(pdbId2);

        Chain c1 = s1.getPolyChain(asymId1);
        Chain c2 = s2.getPolyChain(asymId2);

        StrucAlignment strucAli = computeFatcatAlignment(c1, c2);

        System.out.println("Fatcat alignment:");
        System.out.println(strucAli.seqs.get(0));
        System.out.println(strucAli.seqs.get(1));

        Calc.transform(c1, strucAli.matrices.get(0));
        Calc.transform(c2, strucAli.matrices.get(1));

        GroupSequence atoms1 = new GroupSequence(c1.getAtomGroups());
        GroupSequence atoms2 = new GroupSequence(c2.getAtomGroups());
        GroupSequencePairScorer scorer = new GroupSequencePairScorer(atoms1, atoms2, 8.0);

        long start = System.currentTimeMillis();
        Alignment ali = NeedlemanWunschGotoh.align(atoms1, atoms2, scorer,10f, 1f);
        long end = System.currentTimeMillis();

        // the aligner will swap to have longer sequence first, here we unswap to have the original order
        char[] seq1 = ali.getSequence1();
        char[] seq2 = ali.getSequence2();
        if (atoms1.length() < atoms2.length()) {
            seq1 = ali.getSequence2();
            seq2 = ali.getSequence1();
        }
        System.out.println("Alignment from transform (calculated in " + (end-start) + "ms)");
        System.out.println(seq1);
        System.out.println(ali.getMarkupLine());
        System.out.println(seq2);
    }

    private static class StrucAlignment {
        List<Matrix4d> matrices;
        List<String> seqs;
        private StrucAlignment(List<Matrix4d> matrices, List<String> seqs) {
            this.matrices = matrices;
            this.seqs = seqs;
        }
    }

    private StrucAlignment computeFatcatAlignment(Chain c1, Chain c2) throws StructureException {
        StructureAlignment algorithm  = StructureAlignmentFactory.getAlgorithm(FatCatRigid.algorithmName);
        Atom[] ca1 = StructureTools.getAtomCAArray(c1);
        Atom[] ca2 = StructureTools.getAtomCAArray(c2);
        // Perform the alignment
        AFPChain afpChain = algorithm.align(ca1,ca2);

        // Print text output
        //System.out.println(afpChain.toFatcat(ca1,ca2));

        // assuming 1 block
        Matrix m = afpChain.getBlockRotationMatrix()[0];
        Matrix4d transform = new Matrix4d();
        for (int i=0; i<3; i++) {
            for (int j=0; j<3; j++) {
                transform.setElement(i, j, m.get(j, i));
            }
        }
        Atom[] translation = afpChain.getBlockShiftVector();
        // assuming 1 block
        transform.setTranslation(new Vector3d(translation[0].getX(), translation[0].getY(), translation[0].getZ()));
        transform.setElement(3,3,1);

        List<Matrix4d> transforms = Arrays.asList(ID_MATRIX, transform);
        // for some reason getAlnseq1/2 have a lot of empty 0 chars at end, we trim with the function
        List<String> seqs = Arrays.asList(charArrayToStringTrimmed(afpChain.getAlnseq1()), charArrayToStringTrimmed(afpChain.getAlnseq2()));

        return new StrucAlignment(transforms, seqs);
    }

    private String charArrayToStringTrimmed(char[] array) {
        StringBuilder sb = new StringBuilder();
        for (char c : array) {
            if (c == 0) break;
            sb.append(c);
        }
        return sb.toString();
    }
}
