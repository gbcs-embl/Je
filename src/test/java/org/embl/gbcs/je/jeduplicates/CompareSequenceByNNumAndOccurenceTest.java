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


import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.embl.gbcs.je.jeduplicates.CompareSequenceByNNumAndOccurence;
import org.junit.Assert;
import org.junit.Test;

public class CompareSequenceByNNumAndOccurenceTest {

	
	@Test
	public void sequenceShouldGetSorted(){
		
		Map<String, Integer> counts = new HashMap<String, Integer>();
		counts.put("ATGC", 14);
		counts.put("CCCC", 4);
		counts.put("CCNC", 4);
		counts.put("GGNG", 3);
		counts.put("GGGG", 10);
		counts.put("NNAA", 8);
		counts.put("TNNT", 4);
		counts.put("AAAA", 8);
		counts.put("TTTT", 5);
		counts.put("ANGC", 5);
		
		
		List<String> expectedOrder = Arrays.asList(new String[]{
				"ATGC",
				"GGGG",
				"AAAA",
				"TTTT",
				"CCCC",
				"ANGC",
				"CCNC",
				"GGNG",
				"NNAA",
				"TNNT"
		});
		
		//
		List<String> seqs = new ArrayList<String>(counts.keySet());
		Collections.sort(seqs, new CompareSequenceByNNumAndOccurence(counts));
		
		for (int i = 0; i < seqs.size(); i++) {
			String seq = seqs.get(i);
			Assert.assertTrue("Position "+i+" : got "+seq+" while expected was "+expectedOrder.get(i), seq.equalsIgnoreCase(expectedOrder.get(i)));
		}
		
		
		
		
	}
	
	
	

}














