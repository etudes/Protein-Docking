/*
 * Copyright (c) 2011-2012 Vikram Sundar.
 * All Rights Reserved.
 */
//package org.vikramdock;

import java.io.*;
import java.util.*;
import java.lang.*;

public class DecoyParser {
	public static void main(String[] args) throws IOException {
		PrintWriter out = new PrintWriter(new FileWriter("decoySummary.txt"));
		String[] decoys = {"1ACB", "1ATN", "1AVW", "1AVZ", "1BQL", "1BRC", "1BRS", "1BTH", "1BVK", "1CHO", "1CSE", "1DFJ", "1DQJ", "1EFU", "1EO8", "1FBI", "1FIN", "1FQ1", "1FSS", "1GLA", "1GOT", "1IAI", "1IGC", "1JHL", "1MAH", "1MDA", "1MEL", "1MLC", "1NCA", "1NMB", "1PPE", "1QFU", "1SPB", "1STF", "1TAB", "1TGS", "1UDI", "1UGH", "1WEJ", "1WQ1", "2BTF", "2JEL", "2KAI", "2PCC", "2PTC", "2SIC", "2SNI", "2TEC", "2VIR", "3HHR", "4HTC"};
		for (int i = 0; i < decoys.length; i++) {
			BufferedReader br = new BufferedReader(new FileReader("results\\natives\\" + decoys[i].toLowerCase() + "native.txt"));
			double nativeScore = Double.parseDouble(br.readLine());
			ArrayList<Integer> badDecoys = new ArrayList<Integer>();
			br = new BufferedReader(new FileReader("results\\decoys\\" + decoys[i].toUpperCase() + "decoys.txt"));
			int counter = 0;
			String s = br.readLine();
			while (s != null) {				
				if (Double.parseDouble(s) < nativeScore) badDecoys.add(counter);
				counter++;
				s = br.readLine();
			}
			out.print("DECOY " + decoys[i] + ": ");
			if (badDecoys.size() > 0) {
				for (int j = 0; j < badDecoys.size(); j++) {
					out.print(badDecoys.get(j) + " ");
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
