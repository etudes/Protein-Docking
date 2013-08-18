/*
 * Copyright (c) 2011-2012 Vikram Sundar.
 * All Rights Reserved.
 */
//package org.vikramdock;

import java.io.*;
import java.util.*;
import java.lang.*;
import java.nio.file.*;
import java.nio.charset.*;

public class DrugResultParser implements Comparator<String[]>{
	public static void main(String[] args) throws IOException {
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("results\\drugtnf\\summary.txt")));
		Path drugbank = FileSystems.getDefault().getPath("results\\drugtnf");
		DirectoryStream<Path> files = Files.newDirectoryStream(drugbank);
		ArrayList<String[]> drugs = new ArrayList<String[]>();
		for (Path entry: files) {
			BufferedReader br = Files.newBufferedReader(entry, StandardCharsets.ISO_8859_1);
			System.out.println(entry.toString());
			if (entry.toString().substring(16,18).equals("DB")) {
				while (true) {
					String s = br.readLine();
					if (s == null) break;
					StringTokenizer st = new StringTokenizer(s);
					if (st.nextToken().equals("SIMILARITY")) {
						String[] next = {entry.toString().substring(18,23), st.nextToken()};
						drugs.add(next);
					}
				}
			}
		}
		DrugResultParser drp = new DrugResultParser();
		Collections.sort(drugs, drp);
		for (int i = 0; i < drugs.size(); i++) {
			out.println(drugs.get(i)[0] + " " + drugs.get(i)[1]);
		}
		out.flush();
		out.close();
	}
	public DrugResultParser() {}
	public int compare(String[] s1, String[] s2) {
		double a1 = Double.parseDouble(s1[1]);
		double a2 = Double.parseDouble(s2[1]);
		if (a1 > a2) return 1;
		if (a1 == a2) return 0;
		else return -1;
	}
}
