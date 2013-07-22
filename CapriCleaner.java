/*
 * Copyright (c) 2011-2012 Vikram Sundar.
 * All Rights Reserved.
 */
//package org.vikramdock;

import java.io.*;
import java.util.*;
import java.lang.*;

public class CapriCleaner {
	public static void main(String[] args) throws Exception {
		BufferedReader br = new BufferedReader(new FileReader("capri_" + args[0] + ".brk"));
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("capri_" + args[0] + ".brkFIXED")));
		while (true) {
			String s = br.readLine();
			if (s == null) break;
			if (s.length() > 77 && s.substring(0,4).equals("ATOM") && s.charAt(77) == 'H') continue;
			else out.println(s);
		}
		out.flush();
		System.exit(0);
	}
}

