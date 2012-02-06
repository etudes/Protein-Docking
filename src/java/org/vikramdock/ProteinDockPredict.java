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
	TestCaseStore cases;
	int numthread;
	Thread[] ths;
	public ProteinDockPredict(String file1, String file2) {
		ps1 = new ProteinStruct(file1);
		ps2 = new ProteinStruct(file2);
		TestCaseComparator tcc = new TestCaseComparator();
		cases = new TestCaseStore(Constants.NUMCASES, Constants.NUMCASES, tcc);
		ps1 = ps1.transrot(-ps1.getXCoordCent(), -ps1.getYCoordCent(), -ps1.getZCoordCent(), 0, 0);
		ps2 = ps2.transrot(-ps2.getXCoordCent(), -ps2.getYCoordCent(), -ps2.getZCoordCent(), 0, 0);
	}
	public ProteinDockPredict(ProteinStruct ps1, ProteinStruct ps2) {
		TestCaseComparator tcc = new TestCaseComparator();
		cases = new TestCaseStore(Constants.NUMCASES, Constants.NUMCASES, tcc);
		this.ps1 = ps1.transrot(-ps1.getXCoordCent(), -ps1.getYCoordCent(), -ps1.getZCoordCent(), 0, 0);
		this.ps2 = ps2.transrot(-ps2.getXCoordCent(), -ps2.getYCoordCent(), -ps2.getZCoordCent(), 0, 0);
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
		capriParse(args[0], args[1]);
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
			int curatom = Integer.MAX_VALUE;
			int counter = 0;
			char chain = ' ';
			while(true) {
				String s = br2.readLine();
				if(s != null) {
					String[] ssplit = new String[20];
					ssplit = s.split(" ");
					if (ssplit[0].equals("ATOM")) {
						if (Integer.parseInt(removeSpace(s.substring(6,11))) < curatom) {
							counter++;
							if (counter == 2) {
								chain = s.charAt(21);
							}
						}
						curatom = Integer.parseInt(removeSpace(s.substring(6,11)));
						if (counter == 1) {
							prot1.println(s);
						} else if (counter == 2) {
							prot2.println(s);
						}
					}
				} else {
					break;
				}
			}
			boolean change = false;
			FileInputStream fis2 = new FileInputStream(filepath);
			InputStreamReader isr2 = new InputStreamReader(fis2);
			BufferedReader br3 = new BufferedReader(isr2);
			while(true) {
				String s = br3.readLine();
				if(s != null) {
					String[] ssplit = new String[20];
					ssplit = s.split(" ");
					if (ssplit[0].equals("SEQRES")) {
						if (s.charAt(11) == chain) {
							change = true;
						}
						if (!change) {
							prot1.println(s);
						} else if (change) {
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
			ProteinStruct ps2 = new ProteinStruct(sourcepath.concat("secondprot.txt"));
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
			ps1.printStructurePDB();
			ProteinStruct newps = ps2.transrotpolarall(current.getRmov(), current.getThetamov(), current.getPhimov(), current.getTheta(), current.getPhi());
			newps.printStructurePDB();
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
