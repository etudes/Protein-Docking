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

public class ImmutableArrayList<T> {
	private ArrayList<T> al;
	public ImmutableArrayList(ArrayList<T> al) {
		this.al = al;
	}
	public ImmutableArrayList(ImmutableArrayList<T> clone) {
		al = new ArrayList<T>();
		for (int i = 0; i < clone.size(); i++) {
			T next = clone.get(i);
			T clonea = null;
			if (next instanceof Atom) {
				clonea = (T)new Atom((Atom)next);
			}
			al.add(clonea);
		}
	}
	public T get(int index) {
		return al.get(index);
	}
	public int size() {
		return al.size();
	}
	public ArrayList<T> getAL() {
		return al;
	}
}
