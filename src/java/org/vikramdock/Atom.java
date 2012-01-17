package org.vikramdock;

import java.io.*;
import java.util.*;
import java.util.zip.*;
import java.lang.*;
import java.lang.reflect.*;

public class Atom {
	double xcoord;
	double ycoord;
	double zcoord;
	char element;
	int resnum;
	int atomnum;
	int chainnum;
	boolean spherical;
	ArrayList<Bond> bonded;
	String eType;
	public Atom(double xcoord, double ycoord, double zcoord, char element, int resnum, int atomnum, int chainnum, String eType) {
		this.xcoord = xcoord;
		this.ycoord = ycoord;
		this.zcoord = zcoord;
		this.element = element;
		this.resnum = resnum;
		this.atomnum = atomnum;
		this.chainnum = chainnum;
		this.spherical = false;
		this.eType = eType;
	}
	public Atom(double xcoord, double ycoord, double zcoord, char element, int resnum, int atomnum, int chainnum, String eType, ArrayList bonded) {
		this.xcoord = xcoord;
		this.ycoord = ycoord;
		this.zcoord = zcoord;
		this.element = element;
		this.resnum = resnum;
		this.atomnum = atomnum;
		this.chainnum = chainnum;
		this.spherical = false;
		this.bonded = bonded;
		this.eType = eType;
	}
	public Atom(double xcoord, double ycoord, double zcoord, char element, int resnum, int atomnum, int chainnum, boolean spherical, String eType) {
		this.xcoord = xcoord;
		this.ycoord = ycoord;
		this.zcoord = zcoord;
		this.element = element;
		this.resnum = resnum;
		this.atomnum = atomnum;
		this.chainnum = chainnum;
		this.spherical = spherical;
		this.eType = eType;
	}
	public Atom(double xcoord, double ycoord, double zcoord, char element, int resnum, int atomnum, int chainnum, boolean spherical, String eType, ArrayList bonded) {
		this.xcoord = xcoord;
		this.ycoord = ycoord;
		this.zcoord = zcoord;
		this.element = element;
		this.resnum = resnum;
		this.atomnum = atomnum;
		this.chainnum = chainnum;
		this.spherical = spherical;
		this.bonded = bonded;
		this.eType = eType;
	}
	public Atom(Atom clone) {
		this.xcoord = clone.getXcoord();
		this.ycoord = clone.getYcoord();
		this.zcoord = clone.getZcoord();
		this.element = clone.getElement();
		this.resnum = clone.getResnum();
		this.atomnum = clone.getAtomnum();
		this.chainnum = clone.getChainnum();
		this.spherical = clone.getSpherical();
		this.bonded = clone.getBonded();
		this.eType = clone.getEtype();
	}
	public void addBond(Bond toBeBonded) {
		bonded.add(toBeBonded);
	}
	public double getXcoord() {
		return xcoord;
	}
	public double getYcoord() {
		return ycoord;
	}
	public double getZcoord() {
		return zcoord;
	}
	public double[] getCoords() {
		double[] coords = new double[3];
		coords[0] = xcoord;
		coords[1] = ycoord;
		coords[2] = zcoord;
		return coords;
	}
	public char getElement() {
		return element;
	}
	public int getResnum() {
		return resnum;
	}
	public int getAtomnum() {
		return atomnum;
	}
	public int getChainnum() {
		return chainnum;
	}
	public boolean getSpherical() {
		return spherical;
	}
	public ArrayList getBonded() {
		return bonded;
	}
	public String getEtype() {
		return eType;
	}
	public void printAtom() {
		setCartesian();
		System.out.println("ATOM NO " + atomnum + " FROM RES NO " + resnum + " FROM CHAIN NO " + chainnum + " ELEMENT " + element + " COORDS " + xcoord + " " + ycoord + " " + zcoord);
	}
	public void printAtomErr() {
		System.out.println("ATOM NO " + atomnum + " FROM RES NO " + resnum + " FROM CHAIN NO " + chainnum + " ELEMENT " + element + " COORDS " + xcoord + " " + ycoord + " " + zcoord + " SPHERICAL " + spherical);
	}
	public void printAtomPDB() {
		setCartesian();
		String toBePrinted = "ATOM  ";
		int length = Integer.toString(atomnum).length();
		if (length == 5) { 
			toBePrinted = toBePrinted.concat(Integer.toString(atomnum));
		} else if (length == 4) {
			toBePrinted = toBePrinted.concat(" ").concat(Integer.toString(atomnum));
		} else if (length == 3) {
			toBePrinted = toBePrinted.concat("  ").concat(Integer.toString(atomnum));
		} else if (length == 2) {
			toBePrinted = toBePrinted.concat("   ").concat(Integer.toString(atomnum));
		} else if (length == 1) {
			toBePrinted = toBePrinted.concat("    ").concat(Integer.toString(atomnum));
		}
		toBePrinted = toBePrinted.concat("  ");
		String atom = Character.toString(element);
		if (atom.length() ==  1) {
			toBePrinted = toBePrinted.concat(atom).concat("   ");
		} else if (atom.length() == 2) {
			toBePrinted = toBePrinted.concat(atom).concat("  ");
		} else if (atom.length() == 3) {
			toBePrinted = toBePrinted.concat(atom).concat(" ");
		}
		toBePrinted = toBePrinted.concat("  ");
		toBePrinted = toBePrinted.concat(" ").concat(Integer.toString(chainnum));
		String resnum = Integer.toString(this.resnum);
		if (resnum.length() == 4) {
			toBePrinted = toBePrinted.concat(resnum);
		} else if (resnum.length() == 3) {
			toBePrinted = toBePrinted.concat(" ").concat(resnum);
		} else if (resnum.length() == 2) {
			toBePrinted = toBePrinted.concat("  ").concat(resnum);
		} else if (resnum.length() == 1) {
			toBePrinted = toBePrinted.concat("   ").concat(resnum);
		}
		toBePrinted = toBePrinted.concat("    ");
		String xcoordraw = Double.toString(xcoord);
		String[] xcoords = xcoordraw.split("\\.");
		String xcoordstr = null;
		if (xcoords.length > 1) {
			if (xcoords[1].length() > 3) {
				xcoordstr = xcoords[0].concat(".").concat(xcoords[1].substring(0,3));
			} else {
				xcoordstr = xcoordraw;
			}
		} else {
			xcoordstr = xcoordraw;
		}
		if (xcoordstr.length() == 8) {
			toBePrinted = toBePrinted.concat(xcoordstr);
		} else if (xcoordstr.length() == 7) {
			toBePrinted = toBePrinted.concat(" ").concat(xcoordstr);
		} else if (xcoordstr.length() == 6) {
			toBePrinted = toBePrinted.concat("  ").concat(xcoordstr);
		} else if (xcoordstr.length() == 5) {
			toBePrinted = toBePrinted.concat("   ").concat(xcoordstr);
		} else if (xcoordstr.length() == 4) {
			toBePrinted = toBePrinted.concat("    ").concat(xcoordstr);
		} else if (xcoordstr.length() == 3) {
			toBePrinted = toBePrinted.concat("     ").concat(xcoordstr);
		} else if (xcoordstr.length() == 2) {
			toBePrinted = toBePrinted.concat("      ").concat(xcoordstr);
		} else if (xcoordstr.length() == 1) {
			toBePrinted = toBePrinted.concat("       ").concat(xcoordstr);
		}
		String ycoordraw = Double.toString(ycoord);
		String[] ycoords = ycoordraw.split("\\.");
		String ycoord = null;
		if (ycoords.length > 1) {
			if (ycoords[1].length() > 3) {
				ycoord = ycoords[0].concat(".").concat(ycoords[1].substring(0,3));
			} else {
				ycoord = ycoordraw;
			}
		} else {
			ycoord = ycoordraw;
		}
		if (ycoord.length() == 8) {
			toBePrinted = toBePrinted.concat(ycoord);
		} else if (ycoord.length() == 7) {
			toBePrinted = toBePrinted.concat(" ").concat(ycoord);
		} else if (ycoord.length() == 6) {
			toBePrinted = toBePrinted.concat("  ").concat(ycoord);
		} else if (ycoord.length() == 5) {
			toBePrinted = toBePrinted.concat("   ").concat(ycoord);
		} else if (ycoord.length() == 4) {
			toBePrinted = toBePrinted.concat("    ").concat(ycoord);
		} else if (ycoord.length() == 3) {
			toBePrinted = toBePrinted.concat("     ").concat(ycoord);
		} else if (ycoord.length() == 2) {
			toBePrinted = toBePrinted.concat("      ").concat(ycoord);
		} else if (ycoord.length() == 1) {
			toBePrinted = toBePrinted.concat("       ").concat(ycoord);
		} 
		String zcoordraw = Double.toString(zcoord);
		String[] zcoords = zcoordraw.split("\\.");
		String zcoord = null;
		if (zcoords.length > 1) {
			if (zcoords[1].length() > 3) {
				zcoord = zcoords[0].concat(".").concat(zcoords[1].substring(0,3));
			} else {
				zcoord = zcoordraw;
			}
		} else {
			zcoord = zcoordraw;
		}
		if (zcoord.length() == 8) {
			toBePrinted = toBePrinted.concat(zcoord);
		} else if (zcoord.length() == 7) {
			toBePrinted = toBePrinted.concat(" ").concat(zcoord);
		} else if (zcoord.length() == 6) {
			toBePrinted = toBePrinted.concat("  ").concat(zcoord);
		} else if (zcoord.length() == 5) {
			toBePrinted = toBePrinted.concat("   ").concat(zcoord);
		} else if (zcoord.length() == 4) {
			toBePrinted = toBePrinted.concat("    ").concat(zcoord);
		} else if (zcoord.length() == 3) {
			toBePrinted = toBePrinted.concat("     ").concat(zcoord);
		} else if (zcoord.length() == 2) {
			toBePrinted = toBePrinted.concat("      ").concat(zcoord);
		} else if (zcoord.length() == 1) {
			toBePrinted = toBePrinted.concat("       ").concat(zcoord);
		}
		toBePrinted = toBePrinted.concat("                       ");
		toBePrinted = toBePrinted.concat(atom.substring(0,1));
		System.out.println(toBePrinted);
	}
	public void setSpherical() {
		if (spherical == false) {
			spherical = true;
			double xcoordold = xcoord;
			double ycoordold = ycoord;
			double zcoordold = zcoord;
			xcoord = Math.sqrt((xcoordold * xcoordold) + (ycoordold * ycoordold) + (zcoordold * zcoordold));
			if (xcoord < 0.0d) {
				System.out.println(xcoord + " NEGATIVE " + xcoordold + " " + ycoordold + " " + zcoordold);
				System.out.println(Math.pow(xcoordold, 2) + " " + Math.pow(ycoordold, 2) + " " + Math.pow(zcoordold, 2));
				System.out.flush();
				System.exit(-1);
			}
			ycoord = Math.atan2(ycoordold, xcoordold);
			if (xcoord != 0.0) {
				zcoord = Math.acos(zcoordold/xcoord);
			} else {
				zcoord = 0.0;
			}    
			System.out.flush();
			if (Double.isNaN(xcoord) || Double.isNaN(ycoord) || Double.isNaN(zcoord)) {
				if (xcoord < 0.0d) {
					System.out.println(xcoord + " NEGATIVE " + xcoordold + " " + ycoordold + " " + zcoordold);
					System.out.println(Math.pow(xcoordold, 2) + " " + Math.pow(ycoordold, 2) + " " + Math.pow(zcoordold, 2));
					System.out.flush();
					System.exit(-1);
				}
				System.out.println("NaN created in spherical conversion " + xcoord + " " + ycoord + " " + zcoord + " " + " SPHERICAL " + xcoordold + " " + ycoordold + " " + zcoordold + " CARTESIAN");
				System.out.flush();
			}
		}
	}
	public void setCartesian() {
		if (spherical == true) {
			spherical = false;
			double xcoordold = xcoord;
			double ycoordold = ycoord;
			double zcoordold = zcoord;
			xcoord = xcoordold * Math.cos(ycoordold) * Math.sin(zcoordold);
			ycoord = xcoordold * Math.sin(ycoordold) * Math.sin(zcoordold);
			zcoord = xcoordold * Math.cos(zcoordold);
			System.out.flush();
			if (Double.isNaN(xcoord) || Double.isNaN(ycoord) || Double.isNaN(zcoord)) {
				System.out.println("NaN created in Cartesian conversion " + xcoord + " " + ycoord + " " + zcoord + " " + " CARTESIAN " + xcoordold + " " + ycoordold + " " + zcoordold + " SPHERICAL");
				System.out.flush();
			}
		}
	}
	public Atom transAtom(double xmov, double ymov, double zmov) {
		setCartesian();
		double xcoordnew = xmov + xcoord;
		double ycoordnew = ymov + ycoord;
		double zcoordnew = zmov + zcoord;
		Atom answer = new Atom(xcoordnew, ycoordnew, zcoordnew, element, resnum, atomnum, chainnum, eType, bonded);
		if (Double.isNaN(xcoordnew) || Double.isNaN(ycoordnew) || Double.isNaN(zcoordnew)) {
			System.err.println("NaN created in trans");
			answer.printAtomErr();
		}
		return answer;
	}
	public Atom rotateAtom(double xcent, double ycent, double zcent, double theta, double phi) {
		setCartesian();
		double xcoorda = xcoord - xcent;
		double ycoorda = ycoord - ycent;
		double zcoorda = zcoord - zcent;
		double rcoorda = Math.sqrt(Math.pow(xcoorda, 2) + Math.pow(ycoorda, 2) + Math.pow(zcoorda, 2));
		double thetaa = Math.atan2(ycoorda, xcoorda);
		double phia;
		if (rcoorda != 0.0) {
			phia = Math.acos(zcoorda/rcoorda);
		} else {
			phia = 0;
		}
		double thetab = thetaa + theta;
		if (thetab >= 2*Math.PI) {
			thetab -= 2*Math.PI;
		}
		double phib;
		if (xcoorda >= 0.0) {
			phib = phia + phi;
		} else {
			phib = phia - phi;
		}
		if (phib < 0.0) {
			phib = -phib;
		} if (phib > Math.PI) {
			phib = Math.PI - phib;
		}
		double xcoordb = rcoorda * Math.cos(thetab) * Math.sin(phib) + xcent;
		double ycoordb = rcoorda * Math.sin(thetab) * Math.sin(phib) + ycent;
		double zcoordb = rcoorda * Math.cos(phib) + zcent;
		Atom answer = new Atom(xcoordb, ycoordb, zcoordb, element, resnum, atomnum, chainnum, eType, bonded);
		if (Double.isNaN(xcoordb) || Double.isNaN(ycoordb) || Double.isNaN(zcoordb)) {
			System.err.println("NaN created in rot");
			answer.printAtomErr();
		}
		return answer;
	}
	public double distance(Atom other) {
		boolean changethis = false;
		if (spherical) {
			setCartesian();
			changethis = true;
		}
		boolean changeother = false;
		if (other.getSpherical()) {
			other.setCartesian();
			changeother = true;
		}
		double distance = Math.sqrt(Math.pow(xcoord - other.getXcoord(), 2) + Math.pow(ycoord - other.getYcoord(), 2) + Math.pow(zcoord - other.getZcoord(), 2));
		if (changethis) {
			setSpherical();
		}
		if (changeother) {
			other.setSpherical();
		}
		return distance;
	}
	public double angle(Atom second, Atom third) {
		double a = distance(second);
		double b = distance(third);
		double c = second.distance(third);
		return Math.acos((Math.pow(a,2) + Math.pow(b,2) - Math.pow(c,2))/(2*a*b));
	}
}
