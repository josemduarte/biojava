package org.biojava.nbio.structure.contact;


import java.util.List;

import javax.vecmath.Point3d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.StructureTools;
import org.junit.Test;

public class TestSpatialHashingPerformance {

	public static void main(String[] args) throws Exception {



		double cutoff = 8;

		Structure structure = StructureIO.getStructure("3j3q");		

		System.out.println("Structure loaded");
		
		Atom[] atoms = StructureTools.getAllNonHAtomArray(structure,false);
		//Atom[] atoms = StructureTools.getAtomCAArray(structure);
		
		//Point3d[] points = Calc.atomsToPoints(StructureTools.getAllNonHAtomArray(structure,false));
		Point3d[] points = Calc.atomsToPoints(StructureTools.getAtomCAArray(structure));

		System.out.println("All atoms array created, size: "+atoms.length);
		
		long start = System.currentTimeMillis();

		Grid grid = new Grid(cutoff);
		//grid.addAtoms(atoms);
		//AtomContactSet atomContacts = grid.getAtomContacts();		

		grid.addCoords(points);
		List<Contact> atomContacts = grid.getIndicesContacts();
		
		long end = System.currentTimeMillis();
		
		System.out.println("Found atom contacts: " + atomContacts.size());

		System.out.printf("Calculated contacts in %.3f s\n",((end-start)/1000.0));




	}

}
