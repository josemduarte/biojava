/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava.nbio.structure.contact;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.EntityInfo;
import org.biojava.nbio.structure.Group;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.function.Function;

/**
 * A set of residue-residue contacts.
 * Relies on residue indices to store the pairs and to match contacts.
 *
 * Defining an integer index from a Group can be customised by passing a lambda that computes index from Group.
 *
 * @author Jose Duarte
 */
public class GroupContactSet implements Serializable, Iterable<GroupContact>{

	private static final long serialVersionUID = 1L;

	private static final Logger logger = LoggerFactory.getLogger(GroupContactSet.class);

	private Map<Pair<Integer>, GroupContact> contacts;

	private Function<Group, Integer> indexFunction;

	public GroupContactSet() {
		contacts = new HashMap<>();
		indexFunction = getDefaultIndexFunction();
	}

	/**
	 * Constructs a <code>GroupContactSet</code> by collapsing the given <code>AtomContactSet</code> into
	 * residue-residue (group-group) contacts.
	 * Uses {@link #getDefaultIndexFunction()} to compute indices from Groups
 	 * @param atomContacts a set of atom-atom contacts
	 */
	public GroupContactSet(AtomContactSet atomContacts) {
		contacts = new HashMap<>();
		atoms2groups(atomContacts);
		indexFunction = getDefaultIndexFunction();
	}

	/**
	 * Constructs a <code>GroupContactSet</code> by collapsing the given <code>AtomContactSet</code> into
	 * residue-residue (group-group) contacts.
	 * @param atomContacts a set of atom-atom contacts
	 * @param indexFunction the function to compute indices from Groups
	 */
	public GroupContactSet(AtomContactSet atomContacts, Function<Group, Integer> indexFunction) {
		contacts = new HashMap<>();
		atoms2groups(atomContacts);
		this.indexFunction = indexFunction;
	}

	public void setGroupToIndexFunction(Function<Group, Integer> indexFunction) {
		this.indexFunction = indexFunction;
	}

	private void atoms2groups(AtomContactSet atomContacts) {

		for (AtomContact atomContact:atomContacts) {

			Pair<Atom> atomPair = atomContact.getPair();

			Group iGroup = atomPair.getFirst().getGroup();
			Group jGroup = atomPair.getSecond().getGroup();

			// we skip the self-residue contacts
			if (iGroup.equals(jGroup)) continue;

			Pair<Group> residuePair = new Pair<> (iGroup, jGroup);
			Pair<Integer> pair = new Pair<>(indexFunction.apply(iGroup), indexFunction.apply(jGroup));

			GroupContact groupContact = contacts.computeIfAbsent(pair, k -> new GroupContact(residuePair));

			groupContact.addAtomContact(atomContact);
		}
	}

	public void add(GroupContact groupContact) {
		contacts.put(getResIdPairFromContact(groupContact), groupContact);
	}

	/**
	 * Tell whether the given pair is a contact in this GroupContactSet,
	 * in a chain-identifier independent way: contacts happening between different copies of
	 * the same Compound(Entity) will be considered equal as long as they have the same
	 * residue numbers.
	 * @param group1 first group of contact
	 * @param group2 second group of contact
	 * @return the contact
	 */
	public boolean hasContact(Group group1, Group group2) {
		return contacts.containsKey(new Pair<>(indexFunction.apply(group1), indexFunction.apply(group2)));
	}

	public int size() {
		return contacts.size();
	}

	@Override
	public Iterator<GroupContact> iterator() {
		return contacts.values().iterator();
	}

	private Pair<Integer> getResIdPairFromContact(GroupContact groupContact) {
		return new Pair<>(
				indexFunction.apply(groupContact.getPair().getFirst()),
				indexFunction.apply(groupContact.getPair().getSecond()) );
	}

	/**
	 * Function to get indices based on SEQRES all-chains-within-entity alignment (1-based indices).
	 * @return the function to compute an index from a group
	 */
	public static Function<Group, Integer> getDefaultIndexFunction() {

		return group -> {
			Chain c = group.getChain();
			if (c == null) {
				logger.warn("Chain is not available for group {}. Contact comparison will not work for this residue", group.toString());
				return -1;
			} else {
				EntityInfo comp = c.getEntityInfo();
				if (comp == null) {
					logger.warn("Entity is not available for group {}. Contact comparison will not work for this residue", group.toString());
					return -1;
				} else {
					return comp.getAlignedResIndex(group, c);
				}

			}
		};
	}

}
