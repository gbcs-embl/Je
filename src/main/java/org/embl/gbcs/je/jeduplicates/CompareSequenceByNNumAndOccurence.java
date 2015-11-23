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
package org.embl.gbcs.je.jeduplicates;

import java.util.Comparator;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;

public class CompareSequenceByNNumAndOccurence implements Comparator<String> {

	protected Map<String, Integer> seqCounts;
	
	/**
	 * @param seqCounts a map containing all sequence to be compared (keys) and their occurence number
	 */
	public CompareSequenceByNNumAndOccurence(Map<String, Integer> seqCounts){
		this.seqCounts = seqCounts;
	}

	/** 
	 * Sort sequences by the number of N they contain (with 0 Ns first) then by their occurence number (the highest occurence first)
	 */
	public int compare(String o1, String o2) {
		int numNs1 = StringUtils.countMatches(o1, "N");
		int numNs2 = StringUtils.countMatches(o2, "N");
		int d = numNs1 - numNs2;
		if(d == 0){
			d = seqCounts.get(o2) - seqCounts.get(o1); 
		}
		if(d == 0) return o1.compareTo(o2);
		else return d;
	}

	


}
