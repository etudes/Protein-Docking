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

public class StatReader {
	public static void main(String[] args) throws IOException {
		/**PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("results\\drugtnf\\stats.txt")));
		Path drugbank = FileSystems.getDefault().getPath("results\\drugtnf");
		DirectoryStream<Path> files = Files.newDirectoryStream(drugbank);
		ArrayList<Double> probs = new ArrayList<Double>();
		for (Path entry: files) {
			BufferedReader br = Files.newBufferedReader(entry, StandardCharsets.ISO_8859_1);
			System.out.println(entry.toString());
			if (entry.toString().substring(16,18).equals("DB")) {
				while (true) {
					String s = br.readLine();
					if (s == null) break;
					StringTokenizer st = new StringTokenizer(s);
					if (st.nextToken().equals("PROBABILITY")) probs.add(Double.parseDouble(st.nextToken()));
				}
			}
		}
		BufferedReader br = new BufferedReader(new FileReader("results\\drugtnf\\basesummary.txt"));
		ArrayList<Double> baseprobs = new ArrayList<Double>();
		ArrayList<Integer> basesizes = new ArrayList<Integer>();
		double basesum = 0;
		int cursize = 0;
		while (true) {
			String s = br.readLine();
			if (s == null) break;
			StringTokenizer st = new StringTokenizer(s);
			String first = st.nextToken();
			if (first.equals("SUM")) basesum = Double.parseDouble(st.nextToken());
			else if (first.equals("PROBABILITY")) {
				if (cursize != 0) basesizes.add(cursize);
				baseprobs.add(Double.parseDouble(st.nextToken()));
				cursize = 0;
			}
			else if (first.equals("ATOM")) cursize++;
		}
		basesizes.add(cursize);
		for (int i = 0; i < 100000; i++) {
			double[] testprobs = new double[700];
			double[] testsizes = new double[700];
			double testsum = 0;
			for (int j = 0; j < 700; j++) {
				testprobs[j] = probs.get(random(probs.size()));
				testsizes[j] = random(1223/3);
				testsum += testprobs[j];
			}
			double similarity = 0;
			for (int j = 0; j < 700; j++) {
				for (int k = 0; k < baseprobs.size(); k++) {
					similarity += baseprobs.get(k)/basesum * testprobs[j]/testsum * testsizes[j]/((1223./3.)*(1223./3.));
				}
			}
			out.println(similarity);
			out.flush();
		}
		out.close();
		System.exit(0);**/
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("results\\drugtnf\\statsordered.txt")));
		BufferedReader br = new BufferedReader(new FileReader("results\\drugtnf\\stats.txt"));
		ArrayList<Double> simils = new ArrayList<Double>();
		while (true) {
			String s = br.readLine();
			if (s == null) break;
			simils.add(Double.parseDouble(s));
		}
		Collections.sort(simils);
		for (int i = 0; i < simils.size(); i++) {
			out.println(simils.get(i));
		}
		out.flush();
		out.close();
		System.exit(0);
	}
	public static int random(int size) {
		return (int)(Math.random()*size);
	}
}
