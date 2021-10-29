package org.biojava.nbio.structure.jaligner;

import java.io.Serializable;

public class AminoacidSequence implements Serializable, Sequence<Character> {

    char[] charArray;

    public AminoacidSequence(String sequence) {
        charArray = sequence.toCharArray();
    }


    @Override
    public Character getElement(int i) {
        return charArray[i];
    }

    @Override
    public char getCharAt(int i) {
        return getElement(i);
    }

    @Override
    public int length() {
        return charArray.length;
    }

    @Override
    public String getId() {
        return null;
    }
}
