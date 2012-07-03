/*
 * Copyright (c) 2011-2012 Vikram Sundar.
 * All Rights Reserved.
 */
//package org.vikramdock;

import java.io.*;
import java.util.*;
import java.lang.*;

public class RMSEScoreGen {
	public static void main(String[] args) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(args[0]));
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(args[2])));
		for (int i = 0; i < 1000; i++) {
			double curscore = -1;
			double currms = -1;
			BufferedReader br2 = new BufferedReader(new FileReader(args[1]+"model"+i+".txt"));
			while (true) {
				String s = br2.readLine();
				if (s == null) break;
				if (s.startsWith("SCORE ")) {
					StringTokenizer st = new StringTokenizer(s);
					curscore = Double.parseDouble(st.nextToken());
					break;
				}
			}
			while (true) {
				String s = br.readLine();
				if (s == null) break;
				StringTokenizer st = new StringTokenizer(s);
				if (!st.nextToken().equals("The")) continue;
				st.nextToken();
				if (!st.nextToken().equals("atoms")) continue;
				if (!st.nextToken().equals("have")) continue;
				if (!st.nextToken().equals("an")) continue;
				if (!st.nextToken().equals("RMS")) continue; 
				if (!st.nextToken().equals("distance")) continue;
				if (!st.nextToken().equals("of")) continue;
				currms = Double.parseDouble(st.nextToken());
				
			}
		}
		out.println(curscore + " " + currms);
		out.flush();
	}
}
