package org.biojava.nbio.structure;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.junit.Test;

public class TestStructureIdentifierIssue {

	@Test
	public void testFileNameAsPdbId() throws IOException, StructureException {
		
		// let's create a file name in working dir with pdb id
		File file = new File("1smt");
		file.deleteOnExit();
		PrintWriter pw = new PrintWriter(file);
		pw.println("Nonsense data");
		pw.close();
		
		// now let's request the same PDB id
		Structure s = StructureIO.getStructure("1smt");
		
		// these shouldn't fail
		assertEquals("1SMT", s.getPDBCode());
		assertEquals(2, s.getChains().size());
		
		
		
	}

}
