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
		ProteinDockPredict pdp = new ProteinDockPredict(args[0], args[1], "base", "ABCDEFGHIJKLMNOPQRSTUVWXYZ", args[2], args[3], args[4], args[5], 4, Constants.DRUGNUMCASES, new PrintWriter(new BufferedWriter(new FileWriter("basesummary.txt"))));
		pdp = null;
		ArrayList<ArrayList<Atom>> interfacesBase = new ArrayList<ArrayList<Atom>>();
		ArrayList<Double> probability = new ArrayList<Double>();
		ArrayList<Atom> curInt = new ArrayList<Atom>();
		interfacesBase.add(curInt);
		BufferedReader br = new BufferedReader(new FileReader("basesummary.txt"));
		double sum = 0;
		while (true) {
			String s = br.readLine();
			if (s == null) break;
			StringTokenizer st = new StringTokenizer(s);
			String first = st.nextToken();
			if (first.equals("SUM")) sum = Double.parseDouble(st.nextToken());
			else if (first.equals("PROBABILITY")) probability.add(Double.parseDouble(st.nextToken()));
			else if (first.equals("MODEL")) {
				curInt = new ArrayList<Atom>();
				interfacesBase.add(curInt);
			}
			else if (first.equals("INTERFACE")) continue;
			else {
				double xcoord = Double.parseDouble(removeSpace(s.substring(30,38)));
				double ycoord = Double.parseDouble(removeSpace(s.substring(38,46)));
				double zcoord = Double.parseDouble(removeSpace(s.substring(46,54)));
				char element = ' ';
				if (s.length() > 77) element = s.charAt(77);
				int resnum = Integer.parseInt(removeSpace(s.substring(22,26)));
				int atomnum = Integer.parseInt(removeSpace(s.substring(6,11)));
				char chainnum = s.charAt(21);
				String eType = removeSpace(s.substring(12,16));
				String AA = removeSpace(s.substring(17,20));
				if (element == ' ') element = eType.charAt(0);
				Atom next = new Atom(xcoord, ycoord, zcoord, element, resnum, atomnum, chainnum, eType, AA);
				curInt.add(next);
			}
		}
		br = new BufferedReader(new FileReader("SAMPdrugtest.txt"));
		PriorityQueue<Drug> drugqueue = new PriorityQueue<Drug>();
		while (true) {
			String code = br.readLine();
			if (code == null) break;
			BufferedReader br2 = new BufferedReader(new FileReader(args[0].concat(filesep + "drugbank" + filesep + code + ".txt")));
			String id2 = null;
			while (id2 == null) {
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
					id2 = s.toLowerCase();
					break;
				}
			}
			Drug drugtest = new Drug(code, id2, args[2], args[3], args[0], args[1], interfacesBase, probability, sum);
			/*if (drugtest.getSimilarity() > 0)*/ drugqueue.add(drugtest);
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
