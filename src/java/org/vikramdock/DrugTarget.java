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

public class DrugTarget implements Comparable<DrugTarget> {
	private String code;
	private String id;
	private String target;
	private String targetChain;
	private String target2;
	private String target2Chain;
	private double testSum;
	private double test2Sum;
	private double similarity;
	public DrugTarget(String code, String id, String target, String targetChain, String target2, String target2Chain, String sourcepath, String pdbpath, int numthread) throws Exception {
		this.code = code;
		this.id = id;
		this.target = target;
		this.targetChain = targetChain;
		this.target2 = target2;
		this.target2Chain = target2Chain;
		ProteinDockPredict pdp = new ProteinDockPredict(sourcepath, pdbpath, "prot" + target + "to" + code, "ABCDEFGHIJKLMNOPQRSTUVWXYZ", target, targetChain, id, "ALL", numthread, Constants.DRUGNUMCASES, new PrintWriter(new BufferedWriter(new FileWriter(code + "summary.txt"))));
		BufferedReader br = new BufferedReader(new FileReader(code + "summary.txt"));
		while (true) {
			String s = br.readLine();
			if (s == null) break;
			StringTokenizer st = new StringTokenizer(s);
			String first = st.nextToken();
			if (first.equals("SUM")) testSum = Double.parseDouble(st.nextToken());
		}
		if (testSum < 1000) {
			similarity = 0;
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(code + "summaryall.txt", true)));
			out.println("SIMILARITY " + similarity);
			out.flush();
			out.close();
		} else {
			pdp = new ProteinDockPredict(sourcepath, pdbpath, "prot" + target2 + "to" + code, "ABCDEFGHIJKLMNOPQRSTUVWXYZ", target2, target2Chain, id, "ALL", numthread, Constants.DRUGNUMCASES, new PrintWriter(new BufferedWriter(new FileWriter(code + "summary2.txt"))));
			br = new BufferedReader(new FileReader(code + "summary2.txt"));
			while (true) {
				String s = br.readLine();
				if (s == null) break;
				StringTokenizer st = new StringTokenizer(s);
				String first = st.nextToken();
				if (first.equals("SUM")) testSum = Double.parseDouble(st.nextToken());
			}
			similarity = testSum - test2Sum;
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(code + "summaryall.txt", true)));
			out.println("SIMILARITY " + similarity);
			out.flush();
			out.close();
		}
	}
	public double getSimilarity() {
		return similarity;
	}
	public int compareTo(DrugTarget other) {
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
