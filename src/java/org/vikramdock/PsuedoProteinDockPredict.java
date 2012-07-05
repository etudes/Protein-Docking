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

public class PsuedoProteinDockPredict{
	PsuedoProteinStruct ps1;
	PsuedoProteinStruct ps2;
	PrintWriter out;
	int counter = 0;
	String id;
	public PsuedoProteinDockPredict(PsuedoProteinStruct ps1, PsuedoProteinStruct ps2) {
		this.ps1 = ps1;
		this.ps1 = ps2;
		this.out = out;
		PsuedoTestCase tc = new PsuedoTestCase(ps1, ps2);
		System.out.println(tc.getScore());
	}
	public static void main(String[] args) throws IOException {
		nativeParse(args[0], args[1], args[2], args[3], args[4]);
	}
	public static void nativeParse(String sourcepath, String pdbpath, String complex, String chain1, String chain2) throws IOException {
		try {
			String filesep = System.getProperty("file.separator");
			long start = System.currentTimeMillis();
			PrintWriter prot1 = new PrintWriter(new BufferedWriter(new FileWriter("firstprot.txt")));
			PrintWriter prot2 = new PrintWriter(new BufferedWriter(new FileWriter("secondprot.txt")));
			String filepath = pdbpath.concat(complex.substring(1,3).concat(filesep).concat("pdb").concat(complex).concat(".ent.gz"));
			BufferedReader br2 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(filepath))));
			BufferedReader br3 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(filepath))));
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
						}
					}
					if (ssplit[0].equals("SEQRES")) {
						if (chain1a.contains(s.charAt(11))) {
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
						if (chain2a.contains(s.charAt(21))) {
							prot2.println(s);
						}
					}
					if (ssplit[0].equals("SEQRES")) {
						if (chain2a.contains(s.charAt(11))) {
							prot2.println(s);
						}
					}
				} else {
					break;
				}
			}
			prot1.flush();
			prot2.flush();
			prot1.close();
			prot2.close();
			PsuedoProteinStruct ps1 = new PsuedoProteinStruct(sourcepath.concat("firstprot.txt"));
			System.out.println("DONE WITH FIRST PROTEIN");
			PsuedoProteinStruct ps2 = new PsuedoProteinStruct(sourcepath.concat("secondprot.txt"));
			System.out.println("DONE WITH SECOND PROTEIN");
			PsuedoProteinDockPredict pdp = new PsuedoProteinDockPredict(ps1, ps2);
		} catch (Exception ex) {
			ex.printStackTrace();
		}
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