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

public class Runner {
	public static void main(String[] args) throws Exception {
		String filesep = System.getProperty("file.separator");
		ProteinDockPredict pdp = new ProteinDockPredict(args[0].concat(filesep + "base" + filesep), args[1], "base", "ABCDEFGHIJKLMNOPQRSTUVWXYZ", args[2], args[3], args[4], args[5], 4);
		TestCaseStore<TestCase> cases = pdp.getCases();
		TestCase[] casesA = (TestCase[])cases.toArray();
		Arrays.sort(casesA);
		ArrayList<ArrayList<Atom>> interfacesBase = new ArrayList<ArrayList<Atom>>();
		double sum = 0;
		for (int i = 0; i < Math.min(casesA.length, Constants.NUMCASES); i++) {
			if (casesA[i].getProbability() > 1) {
				interfacesBase.add(casesA[i].getInterfacea());
				sum += casesA[i].getProbability();
			}
		}
		BufferedReader br = new BufferedReader(new FileReader("drugtest.txt"));
		PriorityQueue<Drug> drugqueue = new PriorityQueue<Drug>();
		while (true) {
			String code = br.readLine();
			if (code == null) break;
			BufferedReader br2 = new BufferedReader(new FileReader(args[0].concat(filesep + "drugbank" + filesep + code + ".txt")));
			String id2 = null;
			while (id2 != null) {
				String s = br2.readLine();
				if (s == null) {
					System.err.println("NO PDB ID FOR " + code);
					System.exit(-1);
				}
				StringTokenizer st = new StringTokenizer(s);
				if (st.hasMoreTokens()) st.nextToken();
				if (st.hasMoreTokens() && st.nextToken().equals("PDB_Experimental_ID:")) {
					s = br2.readLine();
					if (s.equals("Not Available")) {
						System.err.println("NO PDB ID FOR " + code);
						System.exit(-1);
					}
					id2 = s;
				}
			}
			Drug drugtest = new Drug(code, id2, args[2], args[3], args[0], args[1], interfacesBase, casesA, sum);
			if (drugtest.getSimilarity() > 0) drugqueue.add(drugtest);
		}
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("drugresults.txt")));
		out.println(drugqueue.size());
		while (!drugqueue.isEmpty()) {
  			out.println(drugqueue.poll());
		}
		out.flush();
		out.close();
		System.exit(0);
	}
}