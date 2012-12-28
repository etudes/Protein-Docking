/*
 * Copyright (c) 2011-2012 Vikram Sundar.
 * All Rights Reserved.
 */	
//package org.vikramdock;

import java.io.*;
import java.util.*;
import java.lang.*;

public class MinRMS {
	public static void main(String[] args) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(args[0]));
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(args[1], true)));
		double minrms = Double.POSITIVE_INFINITY;
		for (int i = 1; i <= 1000; i++) {
			double currms = Double.POSITIVE_INFINITY;
			boolean first = false;
			while (true) {
				String s = br.readLine();
				if (s == null) break;
				StringTokenizer st = new StringTokenizer(s);
				if (!st.hasMoreTokens()) continue;
				if (!st.nextToken().equals("The")) continue;
				st.nextToken();
				if (!st.hasMoreTokens()) continue;
				if (!st.nextToken().equals("atoms")) continue;
				if (!st.hasMoreTokens()) continue;
				if (!st.nextToken().equals("have")) continue;
				if (!st.hasMoreTokens()) continue;
				if (!st.nextToken().equals("an")) continue;
				if (!st.hasMoreTokens()) continue;
				if (!st.nextToken().equals("RMS")) continue;
				if (!st.hasMoreTokens()) continue; 
				if (!st.nextToken().equals("distance")) continue;
				if (!st.hasMoreTokens()) continue;
				if (!st.nextToken().equals("of")) continue;
				if (!st.hasMoreTokens()) continue;
				currms = Double.parseDouble(st.nextToken());
				if (!first) first = true;
				else break;
			}
			if (i == 1) {
				out.print(currms);
				minrms = currms;
			} else {
				minrms = Math.min(minrms, currms);
			}
		}
		out.print(" " + minrms);
		out.flush();
		out.close();
		System.exit(0);
	}
}
