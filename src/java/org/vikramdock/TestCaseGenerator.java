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
	double thetamin;
	double thetamax;
	int num;
	public TestCaseGenerator(ProteinDockPredict pdp, double alphamin, double alphamax, int num) {
		this.pdp = pdp;
		this.ps1 = new ProteinStruct(pdp.ps1);
		this.ps2 = new ProteinStruct(pdp.ps2);
		this.alphamin = alphamin;
		this.alphamax = alphamax;
		this.num = num;
	}
	public void runOld() {
		try {
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("thread" + num + ".txt")));
			for (double i = thetamin; i < thetamax; i += Constants.THETAINC) {
				for (double j = 0; j < Math.PI; j += Constants.PHIINC) {
					for (double k = 0; k < Constants.MAXCLASH; k += Constants.CLASHINC) {
						out.println(i + " " + j + " " + k);
						out.flush();
						for (double theta = 0; theta < 2*Math.PI; theta += Constants.STHETAINC) {
							for (double phi = 0; phi < Math.PI; phi += Constants.SPHIINC) {
								TestCase next = new TestCase(ps1, ps2, i, j, k, theta, phi);
								if (next.score() != Double.POSITIVE_INFINITY) {
									pdp.add(next);
								} else {
									next = null;
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
	public void run() {
		try {
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("thread" + num + ".txt")));
			for (double i = alphamin; i < alphamax; i += Constants.ALPHAINC) {
				for (double j = 0; j < 2*Math.PI; j += Constants.BETAINC) {
					for (double k = 0; k < 2*Math.PI; k += Constants.GAMMAINC) {
						out.println(i + " " + j + " " + k);
						out.flush();
						for (double alpha = 0; alpha < 2*Math.PI; alpha += Constants.SALPHAINC) {
							for (double beta = 0; beta < 2*Math.PI; beta += Constants.SBETAINC) {
								for (double gamma = 0; gamma < 2*Math.PI; gamma += Constants.SGAMMAINC) {
									TestCase next = new TestCase(ps1, ps2, i, j, k, alpha, beta, gamma);
									if (next.score() != Double.POSITIVE_INFINITY) {
										pdp.add(next);
									} else {
										next = null;
									}
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
