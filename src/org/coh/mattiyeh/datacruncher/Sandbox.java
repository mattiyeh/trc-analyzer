package org.coh.mattiyeh.datacruncher;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

public class Sandbox {

	public static void main(String[] args) {
	    List<Double> latencies = new ArrayList<>();
	    latencies.add(1d);
	    latencies.add(2d);
	    latencies.add(3d);
	    latencies.add(4d);
	    latencies.add(6d);
	    latencies.add(5d);
	    latencies.add(7d);
	    latencies.add(8d);
	    latencies.add(9d);
	    latencies.add(10d);
	    System.out.println(latencies);
	    //Collections.sort(latencies);
	    latencies.sort(Comparator.naturalOrder());
	    System.out.println(latencies);

	    System.out.println(percentile(latencies, 25));
	    System.out.println(percentile(latencies, 50));
	    System.out.println(percentile(latencies, 75));
	    System.out.println(percentile(latencies, 100));
	}

	public static double percentile(List<Double> latencies, double percentile) {
	    int index = (int) Math.ceil(percentile / 100.0 * latencies.size());
	    return latencies.get(index-1);
	}

}
