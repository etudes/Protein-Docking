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
	double xmov;
	double ymov;
	double zmov;
	double theta;
	double phi;
	ArrayList<Atom> ps1struct;
	ArrayList<Atom> ps2struct;
	ArrayList<Atom> surfacebya;
	ArrayList<ProteinStruct> surfacebyps;
	public TestCase(ProteinStruct ps1, ProteinStruct ps2, double xmov, double ymov, double zmov, double theta, double phi) {
		surfacebya = new ArrayList<Atom>();
		surfacebyps = new ArrayList<ProteinStruct>();
		this.ps1 = ps1;
		this.ps2 = ps2;
		this.xmov = xmov;
		this.ymov = ymov;
		this.zmov = zmov;
		this.theta = theta;
		this.phi = phi;
		ps1struct = ps1.getSurface();
		surfacebyps.add(ps1);
		newps = ps2.transrot(xmov, ymov, zmov, theta, phi);
		ps2struct = newps.getSurface();
		surfacebyps.add(newps);
		surfacebya.addAll(ps1.getSurface());
		surfacebya.addAll(newps.getSurface());
	}
	public double score() {
		if (!surfaceScore()) {
			return Double.POSITIVE_INFINITY;
		} else {
			return energyScore();
		}
	}
	public boolean surfaceScore() {
		boolean answer = false;
		int space = Constants.SPACE;
		double surfacesize = Constants.SURFACESIZE;
		if (distance(ps1.getCent(), ps2.getCent()) > ps1.getSize() + ps2.getSize() + surfacesize) {
			return false;
		}
		ArrayList structurenew1 = new ArrayList<Atom>();
		for (int i = 0; i < ps1struct.size(); i++) {
			Atom current = (Atom)ps1struct.get(i);
			current.setSpherical();
			structurenew1.add(current);
		}
		ArrayList structurenew2 = new ArrayList<Atom>();
		for (int i = 0; i < ps2struct.size(); i++) {
			Atom current = (Atom)ps2struct.get(i);
			current.setSpherical();
			structurenew2.add(current);
		}
		Atom first = null;
		Atom second = null;
		for (double i = 0; i < 2*Math.PI; i += space*Math.PI/180) {
			for (double j = space*Math.PI/360; j < Math.PI; j += space*Math.PI/180) {
				for (int k = 0; k < structurenew1.size(); k++) {
					Atom current = (Atom)structurenew1.get(k);
					if (Math.abs(current.getYcoord() - i) <= space*Math.PI/360 && Math.abs(current.getZcoord() - j) <= space*Math.PI/360 ) {
						first = current;
						break;
					}
				}
				for (int k = 0; k < structurenew2.size(); k++) {
					Atom current = (Atom)structurenew2.get(k);
					if (Math.abs(current.getYcoord() - i) <= space*Math.PI/360 && Math.abs(current.getZcoord() - j) <= space*Math.PI/360 ) {
						second = current;
						break;
					}
				}
				if (first != null && second != null) {
					double firstrad = first.getXcoord();
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
		double hBondEtot = 0;
		double bStretchEtot = 0;
		double aBendEtot = 0;
		double torsEtot = 0;
		for (int i = 0; i < ps1struct.size(); i++) {
			for (int j = 0; j < ps2struct.size(); j++) {
				Atom first = ps1struct.get(i);
				Atom second = ps2struct.get(j);
				double distance = first.distance(second);
				char firste = first.getElement();
				char seconde = second.getElement();
				if (distance <= Constants.VDWDISTTHRESHOLD) {
					double vanDerWaalsE = 0; 
					if (firste == 'C' && seconde == 'C') {
						vanDerWaalsE = Constants.C12_C_C/Math.pow(distance, 12) - Constants.C6_C_C/Math.pow(distance, 6);
					}
					if ((firste == 'C' && seconde == 'N') || (firste == 'N' && seconde == 'C')) {
						vanDerWaalsE = Constants.C12_C_N/Math.pow(distance, 12) - Constants.C6_C_N/Math.pow(distance, 6);
					}
					if ((firste == 'C' && seconde == 'O') || (firste == 'O' && seconde == 'C')) {
						vanDerWaalsE = Constants.C12_C_O/Math.pow(distance, 12) - Constants.C6_C_O/Math.pow(distance, 6);
					}
					if ((firste == 'C' && seconde == 'S') || (firste == 'S' && seconde == 'C')) {
						vanDerWaalsE = Constants.C12_C_S/Math.pow(distance, 12) - Constants.C6_C_S/Math.pow(distance, 6);
					}
					if ((firste == 'C' && seconde == 'H') || (firste == 'H' && seconde == 'C')) {
						vanDerWaalsE = Constants.C12_C_H/Math.pow(distance, 12) - Constants.C6_C_H/Math.pow(distance, 6);
					}
					if (firste == 'N' && seconde == 'N') {
						vanDerWaalsE = Constants.C12_N_N/Math.pow(distance, 12) - Constants.C6_N_N/Math.pow(distance, 6);
					}
					if ((firste == 'N' && seconde == 'O') || (firste == 'O' && seconde == 'N')) {
						vanDerWaalsE = Constants.C12_N_O/Math.pow(distance, 12) - Constants.C6_N_O/Math.pow(distance, 6);
					}
					if ((firste == 'N' && seconde == 'S') || (firste == 'S' && seconde == 'N')) {
						vanDerWaalsE = Constants.C12_N_S/Math.pow(distance, 12) - Constants.C6_N_S/Math.pow(distance, 6);
					}
					if ((firste == 'N' && seconde == 'H') || (firste == 'H' && seconde == 'N')) {
						vanDerWaalsE = Constants.C12_N_H/Math.pow(distance, 12) - Constants.C6_N_H/Math.pow(distance, 6);
					}
					if (firste == 'O' && seconde == 'O') {
						vanDerWaalsE = Constants.C12_O_O/Math.pow(distance, 12) - Constants.C6_O_O/Math.pow(distance, 6);
					}
					if ((firste == 'O' && seconde == 'S') || (firste == 'S' && seconde == 'O')) {
						vanDerWaalsE = Constants.C12_O_S/Math.pow(distance, 12) - Constants.C6_O_S/Math.pow(distance, 6);
					}
					if ((firste == 'O' && seconde == 'H') || (firste == 'H' && seconde == 'O')) {
						vanDerWaalsE = Constants.C12_O_H/Math.pow(distance, 12) - Constants.C6_O_H/Math.pow(distance, 6);
					}
					if (firste == 'S' && seconde == 'S') {
						vanDerWaalsE = Constants.C12_S_S/Math.pow(distance, 12) - Constants.C6_S_S/Math.pow(distance, 6);
					}
					if ((firste == 'S' && seconde == 'H') || (firste == 'H' && seconde == 'S')) {
						vanDerWaalsE = Constants.C12_S_H/Math.pow(distance, 12) - Constants.C6_S_H/Math.pow(distance, 6);
					}
					if (firste == 'H' && seconde == 'H') {
						vanDerWaalsE = Constants.C12_H_H/Math.pow(distance, 12) - Constants.C6_H_H/Math.pow(distance, 6);
					}
					if (vanDerWaalsE >= Constants.VANDERWAALSTHRESHOLD) {
						return Double.POSITIVE_INFINITY;
					}
					vanDerWaalsEtot += vanDerWaalsE; 
				}
				if (firste == 'H' && distance <= Constants.HBDISTTHRESHOLD) {
					double hBondE = 0;
					double mindist = Constants.BONDDISTHRES;
					int bondedto = -1;
					for (int k = 0; k < ps1struct.size(); k++) {
						Atom current = (Atom)ps1struct.get(k);
						if (current.getElement() != 'H') {
							double curdistance = first.distance(current);
							if (mindist >= curdistance) {
								bondedto = k;
								mindist = curdistance;
							}
						}
					}
					Atom third = (Atom)ps1struct.get(bondedto);
					char thirde = third.getElement();
					double angle = first.angle(second, third);
					if (angle >= Math.PI/2 && thirde != 'C') {
						if (seconde == 'N') {
							hBondE = Constants.H_C12_N_H/Math.pow(distance, 12) - Constants.H_C10_N_H/Math.pow(distance, 10);
						} 
						if (seconde == 'O') {
							hBondE = Constants.H_C12_O_H/Math.pow(distance, 12) - Constants.H_C10_O_H/Math.pow(distance, 10);
						}
						if (seconde == 'S') {
							hBondE = Constants.H_C12_S_H/Math.pow(distance, 12) - Constants.H_C10_S_H/Math.pow(distance, 10);
						}
						hBondE *= -Math.cos(angle);
						if (hBondE >= Constants.HBONDTHRESHOLD) {
							return Double.POSITIVE_INFINITY;
						}
						hBondEtot += hBondE;
					}
				} else if (seconde == 'H' && distance <= Constants.HBDISTTHRESHOLD) {
					double hBondE = 0;
					double mindist = Constants.BONDDISTHRES;
					int bondedto = -1;
					for (int k = 0; k < ps2struct.size(); k++) {
						Atom current = (Atom)ps2struct.get(k);
						if (current.getElement() != 'H') {
							double curdistance = first.distance(current);
							if (mindist >= curdistance) {
								bondedto = k;
								mindist = curdistance;
							}
						}
					}
					Atom third = (Atom)ps2struct.get(bondedto);
					char thirde = third.getElement();
					double angle = second.angle(first, third);
					if (angle >= Math.PI/2 && thirde != 'C') {
						if (firste == 'N') {
							hBondE = Constants.H_C12_N_H/Math.pow(distance, 12) - Constants.H_C10_N_H/Math.pow(distance, 10);
						} 
						if (firste == 'O') {
							hBondE = Constants.H_C12_O_H/Math.pow(distance, 12) - Constants.H_C10_O_H/Math.pow(distance, 10);
						}
						if (firste == 'S') {
							hBondE = Constants.H_C12_S_H/Math.pow(distance, 12) - Constants.H_C10_S_H/Math.pow(distance, 10);
						}
						//extremely dangerous assumption follows here
						hBondE *= -Math.cos(angle);
						if (hBondE >= Constants.HBONDTHRESHOLD) {
							return Double.POSITIVE_INFINITY;
						}
						hBondEtot += hBondE;
					}
				}
			}
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
			bStretchEtot += bStretchE;
		}
		ArrayList ps1backbone = ps1.getBackbone();
		for (int i = 0; i < ps1backbone.size(); i++) {
			double aBendE = 0;
		}
		Etot = Constants.VDWSCALE * vanDerWaalsEtot + Constants.HBONDSCALE * hBondEtot + Constants.BSTRETCHSCALE * bStretchEtot + Constants.ABENDSCALE * aBendEtot + Constants.TORSSCALE * torsEtot;
		return Etot;
	}
	public double getXmov() {
		return xmov;
	}
	public double getYmov() {
		return ymov;
	}
	public double getZmov() {
		return zmov;
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
	public double distance(double[] first, double[] second) {
		return Math.sqrt(Math.pow(first[0] - second[0], 2) + Math.pow(first[1] - second[1], 2) + Math.pow(first[2] - second[2], 2));
	}
	public void printSurfacebya() {
		for (int i = 0; i < surfacebya.size(); i++) {
			Atom current = (Atom)surfacebya.get(i);
			current.printAtomPDB();
		}
	}
}
