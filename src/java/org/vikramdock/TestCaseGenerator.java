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
	int xmin;
	int xmas;
	int xinc;
	int num;
	public TestCaseGenerator(ProteinDockPredict pdp, int xmin, int xmas, int xinc, int num) {
		this.pdp = pdp;
		this.ps1 = new ProteinStruct(pdp.ps1);
		this.ps2 = new ProteinStruct(pdp.ps2);
		this.xmin = xmin;
		this.xmas = xmas;
		this.xinc = xinc;
		this.num = num;
	}
	public void run() {
		try {
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("thread" + num + ".txt")));
			for (double i = xmin; i < xmas; i += xinc) {
				for (double j = -100; j < 100; j += 5) {
					for (double k = -100; k < 100; k += 5) {
						out.println(i + " " + j + " " + k);
						out.flush();
						for (double theta = 0; theta < 2*Math.PI; theta += Math.PI/5) {
							for (double phi = 0; phi < Math.PI; phi += Math.PI/5) {
								TestCase next = new TestCase(ps1, ps2, i, j, k, theta, phi);
								if (next.score() <= Constants.SCORETHRES) {
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
}
