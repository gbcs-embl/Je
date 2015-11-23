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
package org.embl.gbcs.je.jemultiplexer;

/** Utility class to hang onto data about the best match for a given barcode */
public class BarcodeMatch {
	
		/**
		 * indicates if a barcode match has been found, in which case 'barcode' is not null 
		 */
		public boolean matched;
		
		/**
		 * sequence of the matched barcode 
		 */
		public String barcode;
		
		/**
		 * number of mismatches with 'barcode'
		 */
		public int mismatches;
		
		/**
		 * number of mismatches with second best matched barcode
		 */
		public int mismatchesToSecondBest;
		
		public String toString(){
			if(matched)
				return "matched :"+ barcode+" [MM="+mismatches+", MMD="+mismatchesToSecondBest+"]";
			return "no match";
		}
	
}
