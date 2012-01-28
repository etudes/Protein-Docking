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

public class TestCaseComparator implements Comparator {
	public TestCaseComparator() {
	}
	public int compare(Object o1, Object o2) {
		TestCase tc1 = (TestCase)o1;
		TestCase tc2 = (TestCase)o2;
		if (tc1.getScore() > tc2.getScore()) {
			return 1;
		} else if (tc1.getScore() == tc2.getScore()) {
			return 0;
		} else {
			return -1;
		}
	}
}
