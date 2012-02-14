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

public class TestCase {
	ProteinStruct ps1;
	ProteinStruct ps2;
	ProteinStruct newps;
	ProteinStruct newps1;
	ProteinStruct newps2;
	double rmov;
	double thetamov;
	double phimov;
	double alphamov;
	double betamov;
	double gammamov;
	double theta;
	double clash;
	double phi;
	double alpha;
	double beta;
	double gamma;
	ArrayList<Atom> ps1struct;
	ArrayList<Atom> ps2struct;
	ArrayList<Atom> surfacebya;
	ArrayList<ProteinStruct> surfacebyps;
	double score;
	double[][] ps1sizes;
	double[][] ps2sizes;
	double[][][] ps1newsizes;
	double[][][] ps2newsizes;
	public TestCase(ProteinStruct ps1, ProteinStruct ps2, double thetamov, double phimov, double clash, double theta, double phi) {
		surfacebya = new ArrayList<Atom>();
		surfacebyps = new ArrayList<ProteinStruct>();
		this.ps1 = ps1;
		this.ps2 = ps2;
		this.thetamov = thetamov;
		this.phimov = phimov;
		this.clash = clash;
		this.theta = theta;
		this.phi = phi;
		ps1struct = ps1.getSurface();
		surfacebyps.add(ps1);
		ps1sizes = ps1.getSizes();
		ps2sizes = ps2.getSizes();
		rmov = Math.abs(ps1sizes[(int)(thetamov/Constants.THETAINC)][(int)(phimov/Constants.PHIINC)] + ps2sizes[(int)(theta/Constants.THETAINC)][(int)(phi/Constants.PHIINC)] - clash);
		newps = ps2.transrotpolar(rmov, thetamov + Constants.THETAINC/2, phimov + Constants.PHIINC/2, theta, phi);
		ps2struct = newps.getSurface();
		surfacebyps.add(newps);
		surfacebya.addAll(ps1.getSurface());
		surfacebya.addAll(newps.getSurface());
		score = -1;
	}
	public TestCase(ProteinStruct ps1, ProteinStruct ps2, double alphamov, double betamov, double gammamov, double alpha, double beta, double gamma) {
		surfacebya = new ArrayList<Atom>();
		surfacebyps = new ArrayList<ProteinStruct>();
		this.ps1 = ps1;
		this.ps2 = ps2;
		this.alphamov = alphamov;
		this.betamov = betamov;
		this.gammamov = gammamov;
		this.alpha = alpha;
		this.beta = beta;
		this.gamma = gamma;
		ps1newsizes = ps1.getNewSizes();
		ps2newsizes = ps2.getNewSizes();
		rmov = Math.abs(ps1newsizes[(int)(alphamov/Constants.ALPHAINC)][(int)(betamov/Constants.BETAINC)][(int)(gammamov/Constants.GAMMAINC)] + ps2newsizes[(int)(alpha/Constants.ALPHAINC)][(int)(beta/Constants.BETAINC)][(int)(gamma/Constants.GAMMAINC)]);
		this.newps1 = ps1.transrotnew(0, 0, 0, alphamov, betamov, gammamov);
		this.newps2 = ps2.transrotnew(rmov, 0, 0, alpha, Math.PI - beta, gamma);
		ps1struct = newps1.getSurface();
		ps2struct = newps2.getSurface();
		surfacebyps.add(newps1);
		surfacebyps.add(newps2);
		surfacebya.addAll(newps1.getSurface());
		surfacebya.addAll(newps2.getSurface());
		score = -1;
	}
	public double score() {
		double score;
		//if (!surfaceScore()) {
		//	score = Double.POSITIVE_INFINITY;
		//} else {
			score = energyScore();
		//}
		return score;
	}
	public boolean surfaceScore() {
		boolean answer = false;
		double surfacesize = Constants.SURFACESIZE;
		ArrayList structurenew2 = new ArrayList<Atom>();
		for (int i = 0; i < ps2struct.size(); i++) {
			Atom current = (Atom)ps2struct.get(i);
			current.setSpherical();
			structurenew2.add(current);
		}
		Atom second = null;
		for (double i = Constants.THETAINC/2; i < 2*Math.PI; i += Constants.THETAINC) {
			for (double j = Constants.PHIINC/2; j < Math.PI; j += Constants.PHIINC) {
				for (int k = 0; k < structurenew2.size(); k++) {
					Atom current = (Atom)structurenew2.get(k);
					if (Math.abs(current.getYcoord() - i) <= Constants.THETAINC/2 && Math.abs(current.getZcoord() - j) <= Constants.PHIINC/2) {
						second = current;
						break;
					}
				}
				double firstrad = ps1sizes[(int)((i-Constants.THETAINC/2)/Constants.THETAINC)][(int)((phimov-Constants.PHIINC/2)/Constants.PHIINC)];
				if (firstrad != Double.POSITIVE_INFINITY && second != null) {
					double secondrad = second.getXcoord();
					if (Math.abs(firstrad - secondrad) <= surfacesize) {
						answer = true;
					} else if (secondrad < firstrad - surfacesize) {
						return false;
					}
				}
			}
		}
		return answer;
	}
	public double energyScore() {
		double Etot = 0;
		double vanDerWaalsEtot = 0;
		double bStretchEtot = 0;
		double aBendEtot = 0;
		double torsEtot = 0;
		double[][][][] potentials = ps1.getPotentials();
		for (int i = 0; i < ps2struct.size(); i++) {
			Atom current = (Atom)ps2struct.get(i);
			char cel = current.getElement();
			double cx = current.getXcoord();
			double cy = current.getYcoord();
			double cz = current.getZcoord();
			double vanDerWaalsE = 0;
			double ps1size = ps1.getSize();
			if (Math.abs(cx) < ps1size + Constants.VDWDISTTHRESHOLD && Math.abs(cy) < ps1size + Constants.VDWDISTTHRESHOLD && Math.abs(cz) < ps1size + Constants.VDWDISTTHRESHOLD) {
				double[] rotatedcoords = deRotate(cx, cy, cz, alphamov, betamov, gammamov); 
				int tx = (int)((round(rotatedcoords[0], Constants.GRIDGRAINSIZE) + round(ps1size, Constants.GRIDGRAINSIZE) + Constants.VDWDISTTHRESHOLD)/Constants.GRIDGRAINSIZE);
				int ty = (int)((round(rotatedcoords[1], Constants.GRIDGRAINSIZE) + round(ps1size, Constants.GRIDGRAINSIZE) + Constants.VDWDISTTHRESHOLD)/Constants.GRIDGRAINSIZE);
				int tz = (int)((round(rotatedcoords[2], Constants.GRIDGRAINSIZE) + round(ps1size, Constants.GRIDGRAINSIZE) + Constants.VDWDISTTHRESHOLD)/Constants.GRIDGRAINSIZE);
				if (cel == 'C') {
					vanDerWaalsE = potentials[tx][ty][tz][0];
				} else if (cel == 'N') {
					vanDerWaalsE = potentials[tx][ty][tz][1];
				} else if (cel == 'O') {
					vanDerWaalsE = potentials[tx][ty][tz][2];
				} else if (cel == 'S') {
					vanDerWaalsE = potentials[tx][ty][tz][3];
				} else if (cel == 'H') {
					vanDerWaalsE = potentials[tx][ty][tz][4];
				}
			}
			if (vanDerWaalsE > Constants.ETHRESHOLD) {
				return Double.POSITIVE_INFINITY;
			}
			vanDerWaalsEtot += vanDerWaalsE;
		}
		//bond stretching
		ArrayList ps1bonds = ps1.getSurfaceBonds();
		for (int i = 0; i < ps1bonds.size(); i++) {
			double bStretchE = 0;
			Bond current = (Bond)ps1bonds.get(i);
			Atom first = ps1.getAtomByNum(current.getFirst().getAtomnum());
			Atom second = ps1.getAtomByNum(current.getSecond().getAtomnum());
			double distance = first.distance(second);
			bStretchE = current.bondE(distance);
			if (bStretchE > Constants.ETHRESHOLD) {
				return Double.POSITIVE_INFINITY;
			}
			bStretchEtot += bStretchE;
		}
		ArrayList ps2bonds = newps.getSurfaceBonds();
		for (int i = 0; i < ps2bonds.size(); i++) {
			double bStretchE = 0;
			Bond current = (Bond)ps2bonds.get(i);
			Atom first = newps.getAtomByNum(current.getFirst().getAtomnum());
			Atom second = newps.getAtomByNum(current.getSecond().getAtomnum());
			double distance = first.distance(second);
			bStretchE = current.bondE(distance);
			if (bStretchE > Constants.ETHRESHOLD) {
				return Double.POSITIVE_INFINITY;
			}
			bStretchEtot += bStretchE;
		}
		//angle bending
		ArrayList ps1backbone = ps1.getSurfaceBackbone();
		for (int i = 1; i < ps1backbone.size() - 1; i++) {
			double aBendE = 0;
			Atom current = (Atom)ps1backbone.get(i);
			Atom previous = (Atom)ps1backbone.get(i-1);
			Atom next = (Atom)ps1backbone.get(i+1);
			if (!current.getBonded().contains(previous) || !current.getBonded().contains(next)) {
				break;
			} else {
				double origangle = ps1.getAtomByNum(current.getAtomnum()).angle(ps1.getAtomByNum(previous.getAtomnum()), ps1.getAtomByNum(next.getAtomnum()));
				double newangle = current.angle(previous, next);
				aBendE = 1/2 * Constants.ABENDCONST * Math.pow((newangle - origangle),2);
				if (aBendE > Constants.ETHRESHOLD) {
					return Double.POSITIVE_INFINITY;
				}
				aBendEtot += aBendE;
			}
		}
		ArrayList ps2backbone = newps.getSurfaceBackbone();
		for (int i = 1; i < ps2backbone.size() - 1; i++) {
			double aBendE = 0;
			Atom current = (Atom)ps2backbone.get(i);
			Atom previous = (Atom)ps2backbone.get(i-1);
			Atom next = (Atom)ps2backbone.get(i+1);
			if (!current.getBonded().contains(previous) || !current.getBonded().contains(next)) {
				break;
			} else {
				double origangle = ps2.getAtomByNum(current.getAtomnum()).angle(ps2.getAtomByNum(previous.getAtomnum()), ps2.getAtomByNum(next.getAtomnum()));
				double newangle = current.angle(previous, next);
				aBendE = 1/2 * Constants.ABENDCONST * Math.pow((newangle - origangle),2);
				if (aBendE > Constants.ETHRESHOLD) {
					return Double.POSITIVE_INFINITY;
				}
				aBendEtot += aBendE;
			}
		}
		//torsion
		ArrayList ps1backbonebonds = ps1.getSurfaceBackboneBonds();
		for (int i = 1; i < ps1backbonebonds.size() - 1; i++) {
			double torsE = 0;
			Bond current = (Bond)ps1backbonebonds.get(i);
			Bond previous = (Bond)ps1backbonebonds.get(i-1);
			Bond next = (Bond)ps1backbonebonds.get(i+1);
			if ((current.getFirst() != previous.getFirst() && current.getFirst() != previous.getSecond() && current.getSecond() != previous.getFirst() && current.getSecond() != previous.getSecond()) || (current.getFirst() != next.getFirst() && current.getFirst() != next.getSecond() && current.getSecond() != next.getFirst() && current.getSecond() != next.getSecond())) {
				break;
			} else {
				double torsAngle = current.torsAng(previous, next);
				torsE = 1/2*Constants.TORSCONST*(1 + Math.cos(3*torsAngle));
				if (torsE > Constants.ETHRESHOLD) {
					return Double.POSITIVE_INFINITY;
				}
				torsEtot += torsE;
			}
		}
		ArrayList ps2backbonebonds = newps.getSurfaceBackboneBonds();
		for (int i = 1; i < ps2backbonebonds.size() - 1; i++) {
			double torsE = 0;
			Bond current = (Bond)ps2backbonebonds.get(i);
			Bond previous = (Bond)ps2backbonebonds.get(i-1);
			Bond next = (Bond)ps2backbonebonds.get(i+1);
			if ((current.getFirst() != previous.getFirst() && current.getFirst() != previous.getSecond() && current.getSecond() != previous.getFirst() && current.getSecond() != previous.getSecond()) || (current.getFirst() != next.getFirst() && current.getFirst() != next.getSecond() && current.getSecond() != next.getFirst() && current.getSecond() != next.getSecond())) {
				break;
			} else {
				double torsAngle = current.torsAng(previous, next);
				torsE = 1/2*Constants.TORSCONST*(1 + Math.cos(3*torsAngle));
				if (torsE > Constants.ETHRESHOLD) {
					return Double.POSITIVE_INFINITY;
				}
				torsEtot += torsE;
			}
		}
		Etot = Constants.VDWSCALE * vanDerWaalsEtot + Constants.BSTRETCHSCALE * bStretchEtot + Constants.ABENDSCALE * aBendEtot + Constants.TORSSCALE * torsEtot;
		return Etot;
	}
	public double getRmov() {
		return rmov;
	}
	public double getThetamov() {
		return thetamov;
	}
	public double getPhimov() {
		return phimov;
	}
	public double getClash() {
		return clash;
	}
	public double getTheta() {
		return theta;
	}
	public double getPhi() {
		return phi;
	}
	public ProteinStruct getPS1() {
		return ps1;
	}
	public ProteinStruct getPS2() {
		return ps2;
	}
	public ProteinStruct getNewPS() {
		return newps;
	}
	public ArrayList getSurfacebya() {
		return surfacebya;
	}
	public ArrayList getSurfacebyps() {
		return surfacebyps;
	}
	public double getScore() {
		if (score == -1) {
			score();
		}
		return score;
	}
	public double distance(double[] first, double[] second) {
		return Math.sqrt(Math.pow(first[0] - second[0], 2) + Math.pow(first[1] - second[1], 2) + Math.pow(first[2] - second[2], 2));
	}
	public double round(double number, double roundTo) {
		return (double)(Math.round(number/roundTo)*roundTo);
	}
	public void printSurfacebya() {
		for (int i = 0; i < surfacebya.size(); i++) {
			Atom current = (Atom)surfacebya.get(i);
			current.printAtomPDB();
		}
	}
	public void printInfo() {
		System.out.println(rmov + " " + thetamov + " " + phimov + " " + clash + " " + theta + " " + phi);
	}
	public double[] deRotate(double cx, double cy, double cz, double a, double b, double c) {
		double[] answer = new double[3];
		answer[0] = cx*Math.cos(b)*Math.cos(c) + cy*(Math.cos(c)*Math.sin(a)*Math.sin(b) - Math.cos(a)*Math.sin(c)) + cz*(Math.cos(a)*Math.cos(c)*Math.sin(b) + Math.sin(a)*Math.sin(c));
		answer[1] = cx*Math.cos(b)*Math.sin(c) + cz*(Math.cos(a)*Math.sin(b)*Math.sin(c) - Math.cos(c)*Math.sin(a)) + cy*(Math.cos(a)*Math.cos(c) + Math.sin(a)*Math.sin(b)*Math.sin(c));
		answer[2] = cz*Math.cos(a)*Math.cos(b) + cy*Math.sin(a)*Math.cos(b) - cx*Math.sin(b);
		return answer;
	}
}
