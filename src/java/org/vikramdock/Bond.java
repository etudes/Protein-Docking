package org.vikramdock;

import java.io.*;
import java.util.*;
import java.util.zip.*;
import java.lang.*;
import java.lang.reflect.*;

public class Bond {
	Atom first;
	Atom second;
	double equildist;
	double k;
	public Bond(Atom first, Atom second, double k) {
		if (first != null && second != null) {
			this.first = first;
			this.second = second;
			equildist = first.distance(second);
			this.k = k;
			first.addBond(this);
			second.addBond(this);
			if (equildist > Constants.BONDDISTHRES) {
				System.err.println("BOND TOO LARGE");
				first.printAtom();
				second.printAtom();
			}
		} else if (first == null && second != null) {
			System.err.println("NULL ATOM IN BOND");
			second.printAtom();
		} else if (second == null && first != null) {
			System.err.println("NULL ATOM IN BOND");
			first.printAtom();
		} else {
			System.err.println("DOUBLE NULL ATOM IN BOND");
		}
	}
	public Atom getFirst() {
		return first;
	}
	public Atom getSecond() {
		return second;
	}
	public double getEquilDist() {
		return equildist;
	}
	public double getK() {
		return k;
	}
	public double bondE(double realdist) {
		return 1/2 * k * Math.pow(realdist - equildist, 2);
	}
	public double torsAng(Bond firstb, Bond last) {
		double answer = 0;
		Atom prev = null;
		if (firstb.getFirst() != first && firstb.getFirst() != second) {
			prev = firstb.getFirst();
		} else {
			prev = firstb.getSecond();
		}
		Atom next = null;
		if (last.getFirst() != first && last.getFirst() != second) {
			next = last.getFirst();
		} else {
			next = last.getSecond();
		}
		prev.setCartesian();
		first.setCartesian();
		second.setCartesian();
		next.setCartesian();
		double[] prevc = prev.getCoords();
		double[] firstc = first.getCoords();
		double[] secondc = second.getCoords();
		double[] nextc = next.getCoords();
		double[] norm1 = planeNorm(prevc, firstc, secondc);
		double[] norm2 = planeNorm(nextc, firstc, secondc);
		answer = Math.acos((norm1[0] * norm2[0] + norm1[1] * norm2[1] + norm1[2] * norm2[2])/(Math.sqrt(Math.pow(norm1[0],2) + Math.pow(norm1[1],2) + Math.pow(norm1[2],2)) * Math.sqrt(Math.pow(norm2[0],2) + Math.pow(norm2[1],2) + Math.pow(norm2[2],2))));
		if (answer > Math.PI/2) {
			answer = Math.PI - answer;
		}
		return answer;
	}
	public double[] planeNorm(double[] prevc, double[] firstc, double[] secondc) {
		double[] answer = new double[3];
		double x1 = prevc[0];
		double y1 = prevc[1];
		double z1 = prevc[2];
		double x2 = firstc[0];
		double y2 = firstc[1];
		double z2 = firstc[2];
		double x3 = secondc[0];
		double y3 = secondc[1];
		double z3 = secondc[2];
		answer[0] = y3*(z1 - z2) + y1*(z2 - z3) + y2*(-z1 + z3); 
		answer[1] = x3*(-z1 + z2) + x2*(z1 - z3) + x1*(-z2 + z3); 
		answer[2] = x3*(y1 - y2) + x1*(y2 - y3) + x2*(-y1 + y3);
		return answer;
	}
}
