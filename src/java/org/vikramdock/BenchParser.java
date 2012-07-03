/*
 * Copyright (c) 2011-2012 Vikram Sundar.
 * All Rights Reserved.
 */
//package org.vikramdock;

import java.io.*;
import java.util.*;
import java.lang.*;

public class BenchParser {
	public static void main(String[] args) throws Exception {
		BufferedReader br = new BufferedReader(new FileReader("result"+args[0]+".txt"));
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("result"+args[0]+"model1.txt")));
		int counter = 1;
		while (true) {
			String s = br.readLine();
			if (s == null) {
				break;
			}
			if (s.length() > 5) {
				if (s.substring(0,4).equals("ATOM") || s.substring(0,3).equals("TER")) {
					out.println("");
					out.print(s);
				}
				if (s.substring(0,5).equals("MODEL") && !s.substring(0,6).equals("MODELS")) {
					counter++;
					out = new PrintWriter(new BufferedWriter(new FileWriter("result"+args[0]+"model"+counter+".txt")));
					out.print(s.substring(0,14));
				}
			}
			if (s.length() == 5) {
				if (s.substring(0,4).equals("ATOM") || s.substring(0,3).equals("TER")) {
					out.println("");
					out.print(s);
				}
			}
			if (s.length() == 4) {
				if (s.substring(0,3).equals("TER")) {
					out.println("");
					out.print(s);
				}
			}
		}
		out.flush();
		out.close();
	}
}
