package org.coh.mattiyeh.datacruncher.genome;

import org.apache.commons.lang3.StringUtils;

public class GenomeUtils {
	
	private GenomeUtils() {		
	}

	/**
	 * @param str
	 * @return
	 */
	public static String getReverseComplement(String str) {
		return StringUtils.reverse(StringUtils.replaceChars(str, "acgtACGT", "TGCATGCA"));
	}
	
}
