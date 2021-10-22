/*
 * $Id: Sequence.java,v 1.6 2006/07/27 16:28:24 ahmed Exp $
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

package org.biojava.nbio.structure.jaligner;


import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.StructureTools;

import java.io.Serializable;
import java.util.List;

/**
 * A sequence of atom groups to be used in dynamic programming alignment.
 * 
 * @author Jose Duarte
 */

public class GroupSequence implements Serializable, Sequence<Group> {

	private static final long serialVersionUID = 3256721801357898297L;

	private List<Group> atomGroups;
	private String id;

	/**
	 * Constructor
	 * 
	 * @param atomGroups
	 */
	public GroupSequence(List<Group> atomGroups) {
		this.atomGroups = atomGroups;
	}

	/**
	 * Returns the sequence string
	 * 
	 * @return Returns the sequence
	 */
	public List<Group> getAtomGroups() {
		return atomGroups;
	}

	/**
	 * Sets the sequence string
	 * 
	 * @param atomGroups
	 *            The sequence to set
	 */
	public void setAtomGroups(List<Group> atomGroups) {
		this.atomGroups = atomGroups;
	}

	@Override
	public String getId() {
		return id;
	}

	/**
	 * Sets the sequence id
	 * 
	 * @param id
	 *            The id to set
	 */
	public void setId(String id) {
		this.id = id;
	}

	@Override
	public int length() {
		return this.atomGroups.size();
	}

	/**
	 * Returns the atom groups as an array of CA Atoms
	 * 
	 * @return array of chars.
	 */
	public Atom[] toAtomArray() {
		return atomGroups.stream().map(g->g.getAtom("CA")).toArray(Atom[]::new);
	}

    public String toString() {
        return id + Commons.TAB + atomGroups;
    }

	@Override
	public Group getElement(int i) {
		return atomGroups.get(i);
	}

	@Override
	public char getCharAt(int i) {
		return StructureTools.get1LetterCodeAmino(getElement(i).getPDBName());
	}
}