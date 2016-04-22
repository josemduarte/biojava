package org.biojava.nbio.structure.test.xtal;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.contact.StructureInterface;
import org.biojava.nbio.structure.contact.StructureInterfaceCluster;
import org.biojava.nbio.structure.contact.StructureInterfaceList;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.xtal.CrystalBuilder;
import org.junit.Test;

public class TestNcsOpsXtalBuilder {

	@Test
	public void test1AUY() throws IOException, StructureException {


		AtomCache cache = new AtomCache();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		cache.setFileParsingParams(params);
		cache.setUseMmCif(true);

		StructureIO.setAtomCache(cache);

		Structure s = StructureIO.getStructure("1AUY");
		
		StructureTools.expandNcsOps(s);

		CrystalBuilder cb = new CrystalBuilder(s);
		StructureInterfaceList interfaces = cb.getUniqueInterfaces();
		interfaces.calcAsas(1000, 4, 0);
		interfaces.removeInterfacesBelowArea();

		List<StructureInterfaceCluster> clusters = interfaces.getClusters();

		System.out.println("number of interfaces: "+interfaces.size()); 
		System.out.println("number of interface clusters "+clusters.size()); 
		
		for (StructureInterface interf:interfaces) {
			System.out.printf("%3d %3d %15s %25s %7.2f\n", interf.getId(), interf.getCluster().getId(), interf.getParentChains().getFirst().getChainID()+"+"+interf.getParentChains().getSecond().getChainID(), interf.getTransforms().getSecond(), interf.getTotalArea());
		}
		
		assertTrue(interfaces.size()>3);


	}

}
