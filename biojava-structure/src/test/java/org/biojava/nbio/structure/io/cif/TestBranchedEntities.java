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
package org.biojava.nbio.structure.io.cif;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.io.InputStream;
import java.util.zip.GZIPInputStream;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.EntityType;
import org.biojava.nbio.structure.Structure;
import org.junit.jupiter.api.Test;

/**
 * Tests that chains of entity type 'branched' (carbohydrates) are parsed from mmCIF
 * and kept separate from non-polymer chains.
 */
public class TestBranchedEntities {

	/**
	 * 1jpc: one protein chain (A), three branched mannose chains (B, C, D), one water chain (E)
	 */
	private Structure parse1jpc() throws IOException {
		InputStream inStream = new GZIPInputStream(getClass().getResourceAsStream("/org/biojava/nbio/structure/io/1jpc.cif.gz"));
		return CifStructureConverter.fromInputStream(inStream);
	}

	@Test
	public void testBranchedChainsSeparated() throws IOException {
		Structure s = parse1jpc();

		assertEquals(1, s.getPolyChains().size());
		assertEquals(0, s.getNonPolyChains().size());
		assertEquals(3, s.getBranchedChains().size());
		assertEquals(1, s.getWaterChains().size());
		assertEquals(5, s.getChains().size());

		for (Chain c : s.getBranchedChains()) {
			assertEquals(EntityType.BRANCHED, c.getEntityInfo().getType());
			assertTrue(c.getSeqResGroups().isEmpty());
		}

		Chain b = s.getBranchedChain("B");
		assertNotNull(b);
		assertEquals(3, b.getAtomGroups().size());
		assertEquals("MAN", b.getAtomGroup(0).getPDBName());
		assertEquals(2, s.getBranchedChain("C").getAtomGroups().size());
		assertNull(s.getNonPolyChain("B"));

		assertEquals(1, s.getBranchedChainsByPDB("D").size());
		assertTrue(s.getBranchedChainsByPDB("A").isEmpty());
	}

	@Test
	public void testBranchedChainsCloned() throws IOException {
		Structure s = parse1jpc().clone();
		assertEquals(3, s.getBranchedChains().size());
		assertEquals(0, s.getNonPolyChains().size());
	}

	@Test
	public void testBranchedChainsWrittenToPdb() throws IOException {
		Structure s = parse1jpc();
		long manAtoms = s.getBranchedChains().stream().mapToLong(c -> c.getAtomGroups().stream().mapToLong(g -> g.getAtoms().size()).sum()).sum();
		long manLines = s.toPDB().lines().filter(l -> l.startsWith("HETATM") && l.substring(17, 20).equals("MAN")).count();
		assertTrue(manAtoms > 0);
		assertEquals(manAtoms, manLines);
	}
}
