/*
 * Copyright (c) 2011-2012 Vikram Sundar.
 * All Rights Reserved.
 */
//package org.vikramdock;

import java.io.*;
import java.util.*;
import java.lang.*;

public class DecoyKBParser {
	public static void main(String[] args) throws IOException {
		PrintWriter out = new PrintWriter(new FileWriter("decoyKBSummary.txt"));
		String[] decoys = {"1ATN", "1BVK", "1DFJ", "1DQJ", "1FBI", "1FIN", "1IAI", "1IGC", "1JHL", "1MEL", "1WEJ", "2PCC", "2SIC", "2SNI", "2TEC", "2VIR"};
		BufferedReader sum = new BufferedReader(new FileReader("decoySummary.txt"));
		for (int i = 0; i < decoys.length; i++) {
			BufferedReader br = new BufferedReader(new FileReader("results\\natives\\" + decoys[i].toLowerCase() + "native.txt"));
			double nativeScore = Double.parseDouble(br.readLine());
			ArrayList<Integer> badDecoys = new ArrayList<Integer>();
			br = new BufferedReader(new FileReader("results\\decoys\\" + decoys[i].toUpperCase() + "decoyskb.txt"));
			int counter = 0;
			String s = br.readLine();
			while (s != null) {				
				if (Double.parseDouble(s) < nativeScore) badDecoys.add(counter);
				counter++;
				s = br.readLine();
			}
			out.print("DECOY " + decoys[i] + ": ");
			if (badDecoys.size() > 0) {
				s = sum.readLine();
				while (s != null) {
					if (s.startsWith("DECOY " + decoys[i])) break;
					s = sum.readLine();
				}
				ArrayList<Integer> map = new ArrayList<Integer>();
				StringTokenizer st = new StringTokenizer(s);
				st.nextToken(); st.nextToken();
				while (st.hasMoreTokens()) {
					String test = st.nextToken();
					if (test.equals("FAILED")) break;
					map.add(Integer.parseInt(test));
				}
				for (int j = 0; j < badDecoys.size(); j++) {
					out.print(map.get(badDecoys.get(j)) + " ");
				}
				out.println("FAILED " + badDecoys.size());
			}
			else {
				out.println("SUCCESS");
			}
		}
		out.flush();
		out.close();
		System.exit(0);
	}
}
