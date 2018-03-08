/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package org.embl.gbcs.je;

import java.text.NumberFormat;
import java.util.List;

import org.embl.cg.utilitytools.utils.StringUtil;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import htsjdk.samtools.fastq.FastqRecord;

public class HTSLibTest {
	private static Logger log = LoggerFactory.getLogger(HTSLibTest.class);

	
	
	
	/**
	 * Give a wrong layout and make sure we get an exception 
	 *  
	 */
	public static void main(String[] args) {
		String h = "@D3FCO8P1:178:C1WLBACXX:7:1101:1412:2021 1:N:0:";
		String s = "GGAGAATCAAAGGGCAGGGACGTAATCAATGCGAGTTAANNNNNNNNNNNNNNNNNNNNTNNNNAGTTCANGTGAACANNTNCNNTTNANNNNNNNNNGCN";
		String q = "CCCFFFFFGHHHHJJJJJJIIIHIJJIIJJJJJIJHIJI##############################################################";
		FastqRecord r = new FastqRecord(h, s, "+", q);
		
		for (byte b : r.getBaseQualities()) {
			NumberFormat nf = NumberFormat.getIntegerInstance();
			nf.setMinimumIntegerDigits(2);
			System.out.println(nf.format(b));
		}
	}
	
}
