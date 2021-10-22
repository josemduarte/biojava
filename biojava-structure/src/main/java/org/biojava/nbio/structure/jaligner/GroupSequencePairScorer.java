package org.biojava.nbio.structure.jaligner;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.contact.AtomContactSet;
import org.biojava.nbio.structure.contact.Grid;

public class GroupSequencePairScorer implements SequencePairScorer{


    private final GroupSequence s1;
    private final GroupSequence s2;
    private final double cutoff;
    private Atom[] atoms1;
    private Atom[] atoms2;

    private AtomContactSet atomContacts;

    public GroupSequencePairScorer(GroupSequence s1, GroupSequence s2, double cutoff) {
        this.s1 = s1;
        this.s2 = s2;
        this.cutoff = cutoff;
        calcDistanceMatrix(cutoff);
    }

    private void calcDistanceMatrix(double cutoff) {
        Grid grid = new Grid(cutoff);
        atoms1 = getRepresentativeAtomArray(s1);
        atoms2 = getRepresentativeAtomArray(s2);

        grid.addAtoms(atoms1, atoms2);

        atomContacts = grid.getAtomContacts();
    }

    private static Atom[] getRepresentativeAtomArray(GroupSequence s) {
        // TODO make it more general (P or CA)
        return s.getAtomGroups().stream().map(g->g.getAtom("CA")).toArray(Atom[]::new);
    }

    @Override
    public double score(int i, int j) {

        if (atomContacts.hasContact(atoms1[i], atoms2[j])) {
            // TODO find a more solid inversion procedure
            // TODO is it ok that the max is the cutoff (after that score is 0)
            return cutoff - atomContacts.getContact(atoms1[i], atoms2[j]).getDistance();
        }

        // TODO should we penalise this a bit more and make it negative?
        return 0;
    }
}
