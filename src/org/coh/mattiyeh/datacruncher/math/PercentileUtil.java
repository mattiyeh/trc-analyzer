package org.coh.mattiyeh.datacruncher.math;

import java.util.List;

import org.apache.commons.math3.stat.descriptive.rank.Percentile;

public class PercentileUtil {

	private PercentileUtil() {
		throw new IllegalStateException("Utility class");
	}

	public static double calculateNthPercentile(List<Double> data, int nthPercentile) {
		double[] levelsAsArray = data.stream().mapToDouble(Double::doubleValue).toArray();
		return new Percentile().evaluate(levelsAsArray, nthPercentile);
	}

}
