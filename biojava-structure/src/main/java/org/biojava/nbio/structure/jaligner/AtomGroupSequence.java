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

import java.io.Serializable;
import java.util.List;

/**
 * A sequence of atom groups to be used in dynamic programming alignment.
 * 
 * @author Ahmed Moustafa
 */

public class AtomGroupSequence implements Serializable {

	private static final long serialVersionUID = 3256721801357898297L;

	/**
	 * Sequence
	 */
	private List<Group> atomGroups;

	/**
	 * Sequence id.
	 */
	private String id = null;

	/**
	 * Constructor
	 */
	public AtomGroupSequence() {
		super();
	}

	/**
	 * Constructor
	 * 
	 * @param atomGroups
	 */
	public AtomGroupSequence(List<Group> atomGroups) {
		super();
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

	/**
	 * Returns the sequence id
	 * 
	 * @return Returns the id
	 */
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

	/**
	 * Returns the length of the sequence
	 * 
	 * @return sequence length
	 */
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
    
    /**
     * Returns the sequence id and the sequence string
     * 
     * @return Returns the sequence id and the sequence string
     */
    public String toString() {
        return id + Commons.TAB + atomGroups;
    }
}