package org.coh.mattiyeh.datacruncher.genome;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

public class GenomeExtractor {

	private FastaSequenceIndex fastAIndex;
	private IndexedFastaSequenceFile ifsf;

	/**
	 * @param filename
	 * @throws FileNotFoundException
	 */
	public GenomeExtractor(String filename) {
		this(new File(filename));
	}

	/**
	 * @param fastAFile
	 * @throws FileNotFoundException
	 */
	public GenomeExtractor(File fastAFile) {
		fastAIndex = new FastaSequenceIndex(new File(fastAFile.getAbsolutePath() + ".fai"));
		ifsf = new IndexedFastaSequenceFile(fastAFile, fastAIndex);
	}
	
	/**
	 * @param chr
	 * @param pos
	 * @return
	 */
	public String getBaseAt(String chr, int pos) {
		ReferenceSequence rs = ifsf.getSubsequenceAt(chr, pos, pos);
		return rs.getBaseString();
	}

	/**
	 * @param chr
	 * @param pos
	 * @return
	 * @throws IOException
	 */
	public TriSeq getTrinucleotideContext(String chr, int pos) {
		ReferenceSequence rs = ifsf.getSubsequenceAt("chr" + chr, pos - 1L, pos + 1L);
		String str = rs.getBaseString();
		return new TriSeq(str.substring(0, 1), str.substring(1, 2), str.substring(2, 3));
	}

	/**
	 * @param chr
	 * @return
	 */
	public boolean contigExistsInIndex(String chr) {
		return fastAIndex.hasIndexEntry(chr);
	}

	/**
	 * @throws IOException
	 */
	public void close() throws IOException {
		ifsf.close();
	}

}
