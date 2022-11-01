package org.coh.mattiyeh.datacruncher.genome;

import java.io.Serializable;
import java.util.concurrent.ExecutionException;

public class TriSeq implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 3400328982142783019L;
	
	private String preBase;
	private String refBase;
	private String postBase;

	public TriSeq(String preBase, String refBase, String postBase) {
		this.preBase = preBase;
		this.refBase = refBase;
		this.postBase = postBase;
	}

	public TriSeq(String triSeq) throws ExecutionException {
		if (triSeq.length() != 3) {
			throw new ExecutionException("Ref_Tri does not contain correct number (3) of bases", null);
		}

		this.preBase = String.valueOf(triSeq.charAt(0)).toUpperCase();
		this.refBase = String.valueOf(triSeq.charAt(1)).toUpperCase();
		this.postBase = String.valueOf(triSeq.charAt(2)).toUpperCase();
	}

	public String getPreBase() {
		return preBase;
	}

	public String getRefBase() {
		return refBase;
	}

	public String getPostBase() {
		return postBase;
	}

	@Override
	public String toString() {
		return preBase + refBase + postBase;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((postBase == null) ? 0 : postBase.hashCode());
		result = prime * result + ((preBase == null) ? 0 : preBase.hashCode());
		result = prime * result + ((refBase == null) ? 0 : refBase.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		TriSeq other = (TriSeq) obj;
		if (postBase == null) {
			if (other.postBase != null)
				return false;
		} else if (!postBase.equals(other.postBase))
			return false;
		if (preBase == null) {
			if (other.preBase != null)
				return false;
		} else if (!preBase.equals(other.preBase))
			return false;
		if (refBase == null) {
			if (other.refBase != null)
				return false;
		} else if (!refBase.equals(other.refBase))
			return false;
		return true;
	}

}

