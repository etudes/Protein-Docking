/*
 * Copyright (c) 2011-2012 Vikram Sundar.
 * All Rights Reserved.
 */
package org.vikramdock;

import java.io.*;
import java.util.*;
import java.util.zip.*;
import java.lang.*;
import java.lang.reflect.*;

public class ProteinDockPredict{
	ProteinStruct ps1;
	ProteinStruct ps2;
	ProteinStruct origps1;
	ProteinStruct origps2;
	TestCaseStore cases;
	int numthread;
	Thread[] ths;
	public ProteinDockPredict(String file1, String file2) {
		ps1 = new ProteinStruct(file1);
		ps2 = new ProteinStruct(file2);
		TestCaseComparator tcc = new TestCaseComparator();
		cases = new TestCaseStore(Constants.NUMCASES, Constants.NUMCASES, tcc);
		ps1 = ps1.trans(-ps1.getXCoordCent(), -ps1.getYCoordCent(), -ps1.getZCoordCent());
		ps2 = ps2.trans(-ps2.getXCoordCent(), -ps2.getYCoordCent(), -ps2.getZCoordCent());
	}
	public ProteinDockPredict(ProteinStruct ps1, ProteinStruct ps2) {
		TestCaseComparator tcc = new TestCaseComparator();
		cases = new TestCaseStore(Constants.NUMCASES, Constants.NUMCASES, tcc);
		this.origps1 = ps1;
		this.origps2 = ps2;
		this.ps1 = ps1.trans(-ps1.getXCoordCent(), -ps1.getYCoordCent(), -ps1.getZCoordCent());
		System.out.println("TRANSLATED AND ROTATED FIRST PROTEIN");
		this.ps2 = ps2.trans(-ps2.getXCoordCent(), -ps2.getYCoordCent(), -ps2.getZCoordCent());
		System.out.println("TRANSLATED AND ROTATED SECOND PROTEIN");
	}
	public void genTestCases() {
		try {
			ths = new Thread[numthread];
			for (int i = 0; i < numthread; i++) {
				TestCaseGenerator gen = new TestCaseGenerator(this, 0 + 2*Math.PI*i/numthread, 0 + 2*Math.PI*(i+1)/numthread, i);
				Thread th = new Thread(gen);
				ths[i] = th;
				th.start();
			}
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}
	public static void main(String[] args) {
		benchParse(args[0], args[1], args[2], args[3], args[4], args[5], Integer.parseInt(args[6]));
	}
	public static void PDBParse() {
		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
			String id1 = br.readLine();
			String id2 = br.readLine();
			String file1 = "E:\\Research\\pdb\\wwpdb\\pdb\\".concat(id1.substring(1,3)).concat("\\pdb").concat(id1).concat(".ent.gz");
			String file2 = "E:\\Research\\pdb\\wwpdb\\pdb\\".concat(id2.substring(1,3)).concat("\\pdb").concat(id2).concat(".ent.gz");
			long start = System.currentTimeMillis();
			ProteinDockPredict pdp = new ProteinDockPredict(file1, file2);
			pdp.genTestCases();
			for (int i = 0; i < pdp.numthread; i++) {
				pdp.ths[i].join();
			}
			pdp.printCases();
			long end = System.currentTimeMillis();
			System.out.println(end - start);
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}
	public static void benchParse(String sourcepath, String pdbpath, String id1, String chain1, String id2, String chain2, int numthread) {
		try {
			String filesep = System.getProperty("file.separator");
			boolean chainreq1 = true;
			if (chain1.equals("ALL")) {
				chainreq1 = false;
			}
			boolean chainreq2 = true;
			if (chain2.equals("ALL")) {
				chainreq2 = false;
			}
			long start = System.currentTimeMillis();
			PrintWriter prot1 = new PrintWriter(new BufferedWriter(new FileWriter("firstprot.txt")));
			PrintWriter prot2 = new PrintWriter(new BufferedWriter(new FileWriter("secondprot.txt")));
			String filepath1 = pdbpath.concat(id1.substring(1,3).concat(filesep).concat("pdb").concat(id1).concat(".ent.gz"));
			BufferedReader br2 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(filepath1))));
			String filepath2 = pdbpath.concat(id2.substring(1,3).concat(filesep).concat("pdb").concat(id2).concat(".ent.gz"));
			BufferedReader br3 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(filepath2))));
			ArrayList chain1a = new ArrayList();
			ArrayList chain2a = new ArrayList();
			if (chainreq1) {
				for (int i = 0; i < chain1.length(); i++) {
					chain1a.add(chain1.charAt(i));
				}
			}
			if (chainreq2) {
				for (int i = 0; i < chain2.length(); i++) {
					chain2a.add(chain2.charAt(i));
				}
			}
			while(true) {
				String s = br2.readLine();
				if(s != null) {
					String[] ssplit = new String[20];
					ssplit = s.split(" ");
					if (ssplit[0].equals("ATOM")) {
						if (!chainreq1 || chain1a.contains(s.charAt(21))) {
							prot1.println(s);
						}
					}
					if (ssplit[0].equals("SEQRES")) {
						if (!chainreq1 || chain1a.contains(s.charAt(11))) {
							prot1.println(s);
						}
					}
				} else {
					break;
				}
			}
			while(true) {
				String s = br3.readLine();
				if(s != null) {
					String[] ssplit = new String[20];
					ssplit = s.split(" ");
					if (ssplit[0].equals("ATOM")) {
						if (!chainreq2 || chain2a.contains(s.charAt(21))) {
							prot2.println(s);
						}
					}
					if (ssplit[0].equals("SEQRES")) {
						if (!chainreq2 || chain2a.contains(s.charAt(11))) {
							prot2.println(s);
						}
					}
				} else {
					break;
				}
			}
			prot1.close();
			prot2.close();
			ProteinStruct ps1 = new ProteinStruct(sourcepath.concat("firstprot.txt"));
			System.out.println("DONE WITH FIRST PROTEIN");
			ProteinStruct ps2 = new ProteinStruct(sourcepath.concat("secondprot.txt"));
			System.out.println("DONE WITH SECOND PROTEIN");
			ProteinDockPredict pdp = new ProteinDockPredict(ps1, ps2);
			pdp.numthread = numthread;
			System.out.println(pdp.numthread + " NUMTHREAD");
			pdp.genTestCases();
			for (int i = 0; i < pdp.numthread; i++) {
				pdp.ths[i].join();
			}
			pdp.printCases();
			long end = System.currentTimeMillis();
			System.out.println(end - start);
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}
	public static void capriParse(String sourcepath, String capripath) {
		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
			String filename = br.readLine();
			String filepath = capripath.concat("capri_").concat(filename).concat(".brk");
			long start = System.currentTimeMillis();
			PrintWriter prot1 = new PrintWriter(new BufferedWriter(new FileWriter("firstprot.txt")));
			PrintWriter prot2 = new PrintWriter(new BufferedWriter(new FileWriter("secondprot.txt")));
			FileInputStream fis = new FileInputStream(filepath);
			InputStreamReader isr = new InputStreamReader(fis);
			BufferedReader br2 = new BufferedReader(isr);
			BufferedReader br3 = new BufferedReader(new InputStreamReader(new FileInputStream(sourcepath.concat("chaindata").concat(filename).concat(".txt"))));
			String chain1 = br3.readLine();
			String chain2 = br3.readLine();
			if (br3.readLine() != null) {
				System.err.println("TOO MANY PROTEINS");
				System.exit(-1);
			}
			ArrayList chain1a = new ArrayList();
			ArrayList chain2a = new ArrayList();
			for (int i = 0; i < chain1.length(); i++) {
				chain1a.add(chain1.charAt(i));
			}
			for (int i = 0; i < chain2.length(); i++) {
				chain2a.add(chain2.charAt(i));
			}
			while(true) {
				String s = br2.readLine();
				if(s != null) {
					String[] ssplit = new String[20];
					ssplit = s.split(" ");
					if (ssplit[0].equals("ATOM")) {
						if (chain1a.contains(s.charAt(21))) {
							prot1.println(s);
						} else if (chain2a.contains(s.charAt(21))) {
							prot2.println(s);
						} else {
							System.err.println("UNKNOWN CHAIN FOUND " + s.charAt(21));
							System.exit(-1);
						}
					}
				} else {
					break;
				}
			}
			FileInputStream fis2 = new FileInputStream(filepath);
			InputStreamReader isr2 = new InputStreamReader(fis2);
			BufferedReader br4 = new BufferedReader(isr2);
			while(true) {
				String s = br4.readLine();
				if(s != null) {
					String[] ssplit = new String[20];
					ssplit = s.split(" ");
					if (ssplit[0].equals("SEQRES")) {
						if (chain1a.contains(s.charAt(11))) {
							prot1.println(s);
						} else if (chain2a.contains(s.charAt(11))) {
							prot2.println(s);
						} else {
							System.err.println("UNKNOWN CHAIN FOUND " + s.charAt(11));
							System.exit(-1);
						}
					}
				} else {
					break;
				}
			}
			prot1.close();
			prot2.close();
			ProteinStruct ps1 = new ProteinStruct(sourcepath.concat("firstprot.txt"));
			System.out.println("DONE WITH FIRST PROTEIN");
			ProteinStruct ps2 = new ProteinStruct(sourcepath.concat("secondprot.txt"));
			System.out.println("DONE WITH SECOND PROTEIN");
			ProteinDockPredict pdp = new ProteinDockPredict(ps1, ps2);
			pdp.numthread = Integer.parseInt(br.readLine());
			System.out.println(pdp.numthread + " NUMTHREAD");
			pdp.genTestCases();
			for (int i = 0; i < pdp.numthread; i++) {
				pdp.ths[i].join();
			}
			pdp.printCases();
			long end = System.currentTimeMillis();
			System.out.println(end - start);
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}
	public void printCases() {
		for (int i = 0; i < Math.min(cases.size(),10); i++) {
			TestCase current = (TestCase)cases.poll();
			double score = current.getScore();
			if (i < 9) {
				System.out.println("MODEL        " + (i+1) + " " + score);
			} else {
				System.out.println("MODEL       10 " + score);
			}
			current.printInfo();
			ProteinStruct newps1 = origps1.transall(-origps1.getXCoordCent(), -origps1.getYCoordCent(), -origps1.getZCoordCent()).transrotnewall(0, 0, 0, current.getAlphamov(), current.getBetamov(), current.getGammamov());
			newps1.printStructurePDB();
			ProteinStruct newps2 = origps2.transall(-origps2.getXCoordCent(), -origps2.getYCoordCent(), -origps2.getZCoordCent()).transrotnewall(current.getRmov(), 0, 0, current.getAlpha(), Math.PI - current.getBeta(), current.getGamma());
			newps2.printStructurePDB();
			System.out.println("ENDMDL");
		}
		System.out.println("MODELS PASSED " + cases.size());
	}
	public synchronized void add(TestCase worked) {
		cases.addCap(worked);
	}
	public static String removeSpace(String test) {
		String answer = "";
		for (int i = 0; i < test.length(); i++) {
			if (test.charAt(i) != ' ') {
				answer = answer + test.charAt(i);
			}
		}
		return answer;
	}
}
