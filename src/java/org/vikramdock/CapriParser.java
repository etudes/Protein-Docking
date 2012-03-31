/*
 * Copyright (c) 2011-2012 Vikram Sundar.
 * All Rights Reserved.
 */
package org.vikramdock;

import java.io.*;
import java.util.*;
import java.lang.*;

public class CapriParser {
	public static void main(String[] args) throws Exception {
		String[] targets = {"02","03","04","05","06","07","08","09","10","13","14","15","16","17","18","19","20","21","22","23","25","26","27","29","30","32","36","39","41","46","50"};
		for (int i = 0; i < 29; i++) {
			String filereader = "C:\\Old Computer\\E\\Docking\\Protein-Docking\\results\\" + targets[i] + "\\results" + targets[i] + ".txt";/*"*/
			String filewriter = "C:\\Old Computer\\E\\Docking\\Protein-Docking\\results\\" + targets[i] + "\\results" + targets[i] + "capri.txt";/*"*/
			BufferedReader br = new BufferedReader(new FileReader(filereader));
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(filewriter)));
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
						out.println("");
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
}
