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

public class Drug implements Comparable<Drug> {
	private String code;
	private String id;
	private String target;
	private String targetChain;
	private ArrayList<ArrayList<Atom>> interfacesBase;
	private ArrayList<Double> probabilityBase;
	private double baseSum;
	private ArrayList<ArrayList<Atom>> interfacesTest;
	private ArrayList<Double> probability;
	private double testSum;
	private double similarity;
	public Drug(String code, String id, String target, String targetChain, String sourcepath, String pdbpath, ArrayList<ArrayList<Atom>> interfacesBase, ArrayList<Double> probabilityBase, double baseSum) throws Exception {
		this.code = code;
		this.id = id;
		this.target = target;
		this.targetChain = targetChain;
		this.interfacesBase = interfacesBase;
		this.probabilityBase = probabilityBase;
		this.baseSum = baseSum;
		ProteinDockPredict pdp = new ProteinDockPredict(sourcepath, pdbpath, "prot" + target + "to" + code, "ABCDEFGHIJKLMNOPQRSTUVWXYZ", target, targetChain, id, "ALL", 4, Constants.DRUGNUMCASES, new PrintWriter(new BufferedWriter(new FileWriter(code + "summary.txt"))));
		interfacesTest = new ArrayList<ArrayList<Atom>>();
		probability = new ArrayList<Double>();
		ArrayList<Atom> curInt = new ArrayList<Atom>();
		interfacesTest.add(curInt);
		BufferedReader br = new BufferedReader(new FileReader(code + "summary.txt"));
		while (true) {
			String s = br.readLine();
			if (s == null) break;
			StringTokenizer st = new StringTokenizer(s);
			String first = st.nextToken();
			if (first.equals("SUM")) testSum = Double.parseDouble(st.nextToken());
			else if (first.equals("PROBABILITY")) probability.add(Double.parseDouble(st.nextToken()));
			else if (first.equals("MODEL")) {
				curInt = new ArrayList<Atom>();
				interfacesTest.add(curInt);
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
		similarity = similarity();
	}
	public double similarity() {
		double similarity = 0;
		for (int i = 0; i < interfacesBase.size(); i++) {
			for (int j = 0; j < interfacesTest.size(); j++) {
				if (probability.get(j) < 1) continue;
				if (probabilityBase.get(i) < 1) continue;
				similarity += probabilityBase.get(i)/baseSum * probability.get(j)/testSum * similarityOne(interfacesTest.get(j), interfacesBase.get(i));
			}
		}
		return similarity;
	}
	public double similarityOne(ArrayList<Atom> interfaceOne, ArrayList<Atom> interfaceTwo) {
		if (interfaceOne.size() == 0) return 0;
		int count = 0;
		for (int i = 0; i < interfaceOne.size(); i++) {
			Atom test = interfaceOne.get(i);
			for (int j = 0; j < interfaceTwo.size(); j++) {
				if (test.same(interfaceTwo.get(j))) {
					count++;
					break;
				}
			}
		}
		return ((double)count)/interfaceOne.size();
	}
	public double getSimilarity() {
		return similarity;
	}
	public int compareTo(Drug other) {
		if (similarity > other.getSimilarity()) return 1;
		else if (similarity == other.getSimilarity()) return 0;
		else return -1;
	}
	public String toString() {
		return code + " " + id + " " + similarity;
	}
	public String removeSpace(String test) {
		String answer = "";
		for (int i = 0; i < test.length(); i++) {
			if (test.charAt(i) != ' ') {
				answer = answer + test.charAt(i);
			}
		}
		return answer;
	}
}
