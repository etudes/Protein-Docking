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

public class TestCaseGenerator implements Runnable {
	ProteinDockPredict pdp;
	ProteinStruct ps1;
	ProteinStruct ps2;
	double alphamin;
	double alphamax;
	int num;
	public TestCaseGenerator(ProteinDockPredict pdp, double alphamin, double alphamax, int num) {
		this.pdp = pdp;
		this.ps1 = new ProteinStruct(pdp.ps1);
		this.ps2 = new ProteinStruct(pdp.ps2);
		this.alphamin = alphamin;
		this.alphamax = alphamax;
		this.num = num;
	}
	public void run() {
		try {
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("thread" + num + ".txt")));
			for (double i = alphamin; i < alphamax; i += Constants.ALPHAINC) {
				for (double j = 0; j < 2*Math.PI; j += Constants.BETAINC) {
					for (double k = 0; k < 2*Math.PI; k += Constants.GAMMAINC) {
						out.println(i + " " + j + " " + k);
						out.flush();
						for (double alpha = 0; Math.abs(2*Math.PI - alpha) > Constants.FPPRECISION; alpha += Constants.SALPHAINC) {
							for (double beta = 0; Math.abs(2*Math.PI - beta) > Constants.FPPRECISION; beta += Constants.SBETAINC) {
								for (double gamma = 0; Math.abs(2*Math.PI - gamma) > Constants.FPPRECISION; gamma += Constants.SGAMMAINC) {
									out.println(i + " " + j + " " + k + " " + alpha + " " + beta + " " + gamma);
									out.flush();
									TestCase next = new TestCase(ps1, ps2, i, j, k, alpha, beta, gamma);
									if (next.score() != Double.POSITIVE_INFINITY) {
										pdp.add(next);
									} else {
										next = null;
									}
									//if (next.score() <= 0) {
									//	pdp.printCase(next);
									//} else {
									//	next = null;
									//}
								}
							}
						}
					}
				}
			}
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}
}
