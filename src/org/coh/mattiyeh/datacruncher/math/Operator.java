package org.coh.mattiyeh.datacruncher.math;

public enum Operator {

	LESSTHANOREQUAL() {
		@Override
		public boolean apply(double x1, double x2) {
			return x1 <= x2;
		}
	},
	LESSTHAN {
		@Override
		public boolean apply(double x1, double x2) {
			return x1 < x2;
		}
	},
	GREATERTHANOREQUAL() {
		@Override
		public boolean apply(double x1, double x2) {
			return x1 >= x2;
		}
	},
	GREATERTHAN() {
		@Override
		public boolean apply(double x1, double x2) {
			return x1 > x2;
		}
	};

	public abstract boolean apply(double x1, double x2);

}
