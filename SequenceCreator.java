/*
 * Copyright (c) 2011-2012 Vikram Sundar.
 * All Rights Reserved.
 */
//package org.vikramdock;

import java.io.*;
import java.util.*;
import java.lang.*;

public class SequenceCreator {
	public static void main(String[] args) throws IOException {
		for (int i = 1; i <= 1000; i++) {
			String iString = Integer.toString(i);
			for (int j = iString.length(); j < 4; j++) {
				iString = "0" + iString;
			}
			String filename = "../../../disk2/bound_perturb/".concat(args[0]).concat("/aa").concat(args[0]).concat(".ppk_").concat(iString).concat(".pdb");
			BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(filename)));
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(filename.concat("2"), true)));
			String print = "";
			int curnum = 0;
			char curChainnum = ' ';
			String catom;
			int curresnum = -1;
			while ((catom = br.readLine()) != null) {
				String[] catoms = catom.split(" ");
				if (catoms[0].equals("ATOM")) {
					out.println(catom);
					char chainnum = catom.charAt(21);
					int resnum = Integer.parseInt(removeSpace(catom.substring(22,26)));
					if (chainnum != curChainnum && curChainnum != ' ') {
						print += "\nSEQRES     " + chainnum + "       ";
						curChainnum = chainnum;
						curnum = 0;
					} else if (curChainnum == ' ') {
						print += "SEQRES     " + chainnum + "       ";
						curChainnum = chainnum;
						curnum = 0;
					}
					if (curnum == 13) {
						print += "\nSEQRES     " + chainnum + "       ";
						curnum = 0;
					}
					if (resnum != curresnum) {
						curnum++;
						String AA = removeSpace(catom.substring(17,20));
						print += AA + " ";
						curresnum = resnum;
					}
				}
			}
			out.print(print);
			out.flush();
			out.close();
			File origFile = new File(filename);
			File newFile = new File(filename.concat("2"));
			//origFile.delete();
			newFile.renameTo(origFile);
		}
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
