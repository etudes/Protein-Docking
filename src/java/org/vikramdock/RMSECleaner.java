/*
 * Copyright (c) 2011-2012 Vikram Sundar.
 * All Rights Reserved.
 */
//package org.vikramdock;

import java.io.*;
import java.util.*;
import java.lang.*;

public class RMSECleaner {
	public static void main(String[] args) {
		int num = Integer.parseInt(args[0]);
		BufferedReader br = new BufferedReader(new FileReader("p"+num+"resultstmp.txt"));
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("p"+num+"results.txt")));
		boolean read = false;
		while (true) {
			String s = br.readLine();
			if (s == null) break;
			if (s.startsWith(" LSQMAN >  Spawn system command : ")) {
				read = true;
			}
			if (s.startsWith(" *** LSQMAN *** LSQMAN *** LSQMAN *** LSQMAN *** LSQMAN *** LSQMAN *** ")) {
				read = false;
			}
			if (read) {
				out.println(s);
			}
		}
	}
}
