package org.biojava.nbio.structure.jaligner;

/**
 * A generalised sequence of objects for which there can be similarities defined. To be used in
 * dynamic programming alignment (Smith-Waterman and Needleman-Wunsch).
 */
public interface Sequence<T> {

    T getElement(int i);

    char getCharAt(int i);

    /**
     * Returns the length of the sequence
     * @return sequence length
     */
    int length();

    /**
     * Retrieve a string identifier for this {@link Sequence}
     * @return the identifier
     */
    String getId();

}
