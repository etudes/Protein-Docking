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
}
