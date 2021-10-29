package org.biojava.nbio.structure.jaligner;

import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;

import java.io.Serializable;

/**
 * Note this is similar to SubstitutionMatrixScorer in the core package.
 * TODO review and unify if possible
 */
public class AminoacidSequencePairScorer implements Serializable, SequencePairScorer {

    private static final AminoAcidCompoundSet AMINO_ACID_COMPOUND_SET = new AminoAcidCompoundSet();

    private final AminoacidSequence s1;
    private final AminoacidSequence s2;
    private SubstitutionMatrix<AminoAcidCompound> matrix;

    public AminoacidSequencePairScorer(AminoacidSequence s1, AminoacidSequence s2, SubstitutionMatrix<AminoAcidCompound> matrix) {
        this.s1 = s1;
        this.s2 = s2;
        this.matrix = matrix;
    }

    public void setSimilarityMatrix(SubstitutionMatrix<AminoAcidCompound> matrix) {
        this.matrix = matrix;
    }

    @Override
    public double score(int i, int j) {
        char iChar = s1.getCharAt(i);
        char jChar = s2.getCharAt(j);
        AminoAcidCompound amino1 = AMINO_ACID_COMPOUND_SET.getCompoundForString(String.valueOf(iChar));
        AminoAcidCompound amino2 = AMINO_ACID_COMPOUND_SET.getCompoundForString(String.valueOf(jChar));
        return matrix.getValue(amino1, amino2);
    }
}
