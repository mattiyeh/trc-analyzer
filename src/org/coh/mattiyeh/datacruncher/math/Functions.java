package org.coh.mattiyeh.datacruncher.math;

public class Functions {
	
	private Functions() {
	    throw new IllegalStateException("Utility class");
	  }
	
	public static Double safeLogTransform(double value) {
		return Math.log(value + 1);
	}

}
