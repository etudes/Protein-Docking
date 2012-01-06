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
	int numthread;
	Thread[] ths;
	public ProteinDockPredict(String file1, String file2) {
		ps1 = new ProteinStruct(file1);
		ps2 = new ProteinStruct(file2);
		cases = new ArrayList<TestCase>();
		ps1 = ps1.transrot(-ps1.getXCoordCent(), -ps1.getYCoordCent(), -ps1.getZCoordCent(), 0, 0);
		ps2 = ps2.transrot(-ps2.getXCoordCent(), -ps2.getYCoordCent(), -ps2.getZCoordCent(), 0, 0);
	}
	public ProteinDockPredict(ProteinStruct ps1, ProteinStruct ps2) {
		this.ps1 = ps1;
		this.ps2 = ps2;
		cases = new ArrayList<TestCase>();
		ps1 = ps1.transrot(-ps1.getXCoordCent(), -ps1.getYCoordCent(), -ps1.getZCoordCent(), 0, 0);
		ps2 = ps2.transrot(-ps2.getXCoordCent(), -ps2.getYCoordCent(), -ps2.getZCoordCent(), 0, 0);
	}
	public void genTestCases() {
		try {
			ths = new Thread[numthread];
			for (int i = 0; i < numthread; i++) {
				TestCaseGenerator gen = new TestCaseGenerator(this, -100 + 200*i/numthread, -100 + 200*(i+1)/numthread, 5, i);
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
			//ProteinStruct newps = pdp.ps1;
			//newps.printSurface();
			pdp.genTestCases();
			for (int i = 0; i < pdp.cases.size(); i++) {
				TestCase current = (TestCase)pdp.cases.get(i);
				current.printSurfacebya();
			}
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
			int counter = 0;
			char chain = ' ';
			while(true) {
				String s = br2.readLine();
				if(s != null) {
					String[] ssplit = new String[20];
					ssplit = s.split(" ");
					if (ssplit[0].equals("ATOM")) {
						if (Integer.parseInt(removeSpace(s.substring(6,11))) == 1) {
							counter++;
							if (counter == 2) {
								chain = s.charAt(21);
							}
						}
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
			for (int i = 0; i < pdp.cases.size(); i++) {
				TestCase current = (TestCase)pdp.cases.get(i);
				current.printSurfacebya();
			}
			long end = System.currentTimeMillis();
			System.out.println(end - start);
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}
	public synchronized void add(TestCase worked) {
		cases.add(worked);
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
