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

public class DrugBankParser {
	public static void main(String[] args) throws IOException {
		/**BufferedReader br = new BufferedReader(new FileReader("C:\\Users\\Vikram\\Documents\\GitHub\\Protein-Docking\\drugbank.txt"));
		String current = "";
		PrintWriter out = null;
		while (true) {
			String s = br.readLine();
			if (s == null) break;
			if (current.equals("")) {
				if (out != null) {
					out.flush();
					out.close();
				}
				StringTokenizer st = new StringTokenizer(s);
				st.nextToken();
				current = "C:\\Users\\Vikram\\Documents\\GitHub\\Protein-Docking\\drugbank\\" + st.nextToken() + ".txt";
				out = new PrintWriter(new BufferedWriter(new FileWriter(current)));
			}
			out.println(s);
			out.flush();
			StringTokenizer st = new StringTokenizer(s);
			if (st.hasMoreTokens() && st.nextToken().equals("#END_DRUGCARD")) current = "";
		}**/
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("C:\\Users\\Vikram\\Documents\\GitHub\\Protein-Docking\\drugtest.txt")));
		Path drugbank = FileSystems.getDefault().getPath("drugbank");
		DirectoryStream<Path> files = Files.newDirectoryStream(drugbank);
		for (Path entry: files) {
			BufferedReader br = Files.newBufferedReader(entry, StandardCharsets.ISO_8859_1);
			boolean approved = false;
			boolean pdb = false;
			while (true) {
				String s = br.readLine();
				if (s == null) break;
				StringTokenizer st = new StringTokenizer(s);
				if (st.hasMoreTokens() && st.nextToken().equals("Approved")) approved = true;
				if (st.hasMoreTokens() && st.nextToken().equals("PDB_Experimental_ID:")) {
					s = br.readLine();
					if (!s.equals("Not Available")) pdb = true;
				}
			}
			if (approved && pdb) out.println(entry.toString().substring(9,16));
			out.flush();
		}
		out.close();
	}
}
