/*
 * Copyright (c) 2011-2012 Vikram Sundar.
 * All Rights Reserved.
 */
package org.vikramdock;

import java.io.*;
import java.util.*;
import java.util.zip.*;
import java.lang.*;
import java.lang.reflect.*;	

public class RunnerTarget {
	public static void main(String[] args) throws Exception {
		String filesep = System.getProperty("file.separator");
		BufferedReader br = new BufferedReader(new FileReader("drugtest.txt"));
		PriorityQueue<DrugTarget> drugqueue = new PriorityQueue<DrugTarget>();
		while (true) {
			String code = br.readLine();
			if (code == null) break;
			BufferedReader br2 = new BufferedReader(new FileReader(args[0].concat(filesep + "drugbank" + filesep + code + ".txt")));
			String id2 = null;
			boolean error = false;
			while (id2 == null) {
				String s = br2.readLine();
				if (s == null) {
					System.err.println("NO PDB ID FOR " + code);
					error = true;
				}
				StringTokenizer st = new StringTokenizer(s);
				if (st.hasMoreTokens()) st.nextToken();
				if (st.hasMoreTokens() && st.nextToken().equals("PDB_Experimental_ID:")) {
					s = br2.readLine();
					if (s.equals("Not Available")) {
						System.err.println("NO PDB ID FOR " + code);
						error = true;
					}
					id2 = s.toLowerCase();
					break;
				}
			}
			String filepath2 = args[1].concat(id2.substring(1,3).concat(filesep).concat("pdb").concat(id2).concat(".ent.gz"));
			File test = new File(filepath2);
			if (!test.exists()) error = true;
			if (error) continue;
			DrugTarget drugtest = new DrugTarget(code, id2, args[2], args[3], args[4], args[5], args[0], args[1], Integer.parseInt(args[6]));
			if (drugtest.getSimilarity() > 0) drugqueue.add(drugtest);
		}
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("drugresults.txt")));
		out.println(drugqueue.size());
		while (!drugqueue.isEmpty()) {
  			out.println(drugqueue.poll());
		}
		out.flush();
		out.close();
		System.exit(0);
	}
	public static String removeSpace(String test) {
		String answer = "";
		for (int i = 0; i < test.length(); i++) {
			if (test.charAt(i) != ' ') {
				answer = answer + test.charAt(i);
			}
		}
		return answer;
	}
}
