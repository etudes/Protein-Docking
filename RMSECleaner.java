/*
 * Copyright (c) 2011-2012 Vikram Sundar.
 * All Rights Reserved.
 */
//package org.vikramdock;

import java.io.*;
import java.util.*;
import java.lang.*;

public class RMSECleaner {
	public static void main(String[] args) throws IOException {
		int num = Integer.parseInt(args[0]);
		BufferedReader br = new BufferedReader(new FileReader("p"+num+"resultstmp.txt"));
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("p"+num+"results.txt", true)));
		boolean read = false;
		while (true) {
			String s = br.readLine();
			if (s == null) break;
			String begin = " LSQMAN >  Spawn system command : ";
			String end = " *** LSQMAN *** LSQMAN *** LSQMAN *** LSQMAN *** ";
			if (s.startsWith(begin)) {
				read = true;
			}
			if (s.startsWith(end)) {
				read = false;
			}
			if (s.startsWith("Result for")) {
				out.println(s);
			}
			if (read) {
				out.println(s);
			}
		}
		out.flush();
		out.close();
	}
}
