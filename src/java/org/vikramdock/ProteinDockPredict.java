package org.vikramdock;

import java.io.*;
import java.util.*;
import java.util.zip.*;
import java.lang.*;
import java.lang.reflect.*;

public class ProteinDockPredict{
	ProteinStruct ps1;
	ProteinStruct ps2;
	ArrayList<TestCase> cases;
	public ProteinDockPredict(String id1, String id2) {
		ps1 = new ProteinStruct(id1);
		ps2 = new ProteinStruct(id2);
		cases = new ArrayList<TestCase>();
		ps1 = ps1.transrot(-ps1.getXCoordCent(), -ps1.getYCoordCent(), -ps1.getZCoordCent(), 0, 0);
		ps2 = ps2.transrot(-ps2.getXCoordCent(), -ps2.getYCoordCent(), -ps2.getZCoordCent(), 0, 0);
	}
	public void genTestCases() {
		for (double i = -100; i < 100; i += 5) {
			for (double j = -100; j < 100; j += 5) {
				for (double k = -100; k < 100; k += 5) {
					for (double theta = 0; theta < 2*Math.PI; theta += Math.PI/5) {
						for (double phi = 0; phi < Math.PI; phi += Math.PI/5) {
							TestCase next = new TestCase(ps1, ps2, i, j, k, theta, phi);
							System.out.println(i + " " + j + " " + k + " " + theta + " " + phi);
							if (next.score() <= Constants.SCORETHRES) {
								cases.add(next);
							}
						}
					}
				}
			}
		}
	}
	public static void main(String[] args) {
		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
			String id1 = br.readLine();
			String id2 = br.readLine();
			long start = System.currentTimeMillis();
			ProteinDockPredict pdp = new ProteinDockPredict(id1, id2);
			//ProteinStruct newps = pdp.ps1;
			//newps.printSurface();
			pdp.genTestCases();
			long end = System.currentTimeMillis();
			System.out.println(end - start);
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}
}
