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
	public PsuedoProteinDockPredict(PsuedoProteinStruct ps1, PsuedoProteinStruct ps2, String complex) throws IOException {
		this.ps1 = ps1;
		this.ps1 = ps2;
		this.out = out;
		PsuedoTestCase tc = new PsuedoTestCase(ps1, ps2);
		out = new PrintWriter(new BufferedWriter(new FileWriter(complex + "decoys.txt", true)));
		out.println(tc.getScore());
		out.flush();
		out.close();
	}
	public static void main(String[] args) throws IOException {
		for (int i = 1; i <= 1000; i++) {
			nativeParse(args[0], args[1], args[2], args[3], args[4], i);
			System.out.println(i);
		}
	}
	public static void nativeParse(String sourcepath, String pdbpath, String complex, String chain1, String chain2, int num) throws IOException {
		try {
			String filesep = System.getProperty("file.separator");
			long start = System.currentTimeMillis();
			PrintWriter prot1 = new PrintWriter(new BufferedWriter(new FileWriter("firstprot.txt")));
			PrintWriter prot2 = new PrintWriter(new BufferedWriter(new FileWriter("secondprot.txt")));
			String iString = Integer.toString(num);
			for (int j = iString.length(); j < 4; j++) {
				iString = "0" + iString;
			}
			String filepath = pdbpath.concat(complex.concat(filesep).concat("aa").concat(complex).concat(".ppk_").concat(iString).concat(".pdb"));
			BufferedReader br2 = new BufferedReader(new InputStreamReader(new FileInputStream(filepath)));
			BufferedReader br3 = new BufferedReader(new InputStreamReader(new FileInputStream(filepath)));
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
			PsuedoProteinStruct ps2 = new PsuedoProteinStruct(sourcepath.concat("secondprot.txt"));
			PsuedoProteinDockPredict pdp = new PsuedoProteinDockPredict(ps1, ps2, complex);
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
