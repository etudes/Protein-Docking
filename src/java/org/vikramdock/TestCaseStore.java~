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

public class TestCaseStore extends PriorityQueue {
	int maxcap;
	public TestCaseStore(int maxcap, int initialCapacity) {
		super(initialCapacity);
		this.maxcap = maxcap;
	}
	public TestCaseStore(int maxcap, int initialCapacity, Comparator comparator) {
		super(initialCapacity, comparator);
		this.maxcap = maxcap;
	}
	public boolean addCap(Object o) {
		if (super.size() < maxcap) {
			super.add(o);
		} else {
			reduceSize(maxcap);
			super.add(o);
			reduceSize(maxcap);
		}
		return true;
	}
	public void reduceSize(int cap) {
		if (super.size() > cap) {
			Object[] array = super.toArray();
			int size = super.size(); 
			for (int i = cap; i < size; i++) {
				super.remove(array[i]);
			}
		}
	}
}
	
