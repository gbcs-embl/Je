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

import java.util.Map;

/** Utility class to hang onto data about the best match for a given barcode */
public class SampleMatch {
	
	
		/**
		 * The identified sample name or null 
		 * 
		 */
		protected String sample;
	
		/**
		 * The {@link BarcodeMatch} or null for each barcode slot
		 * A {@link BarcodeMatch} is always present if the sample has been identified
		 */
		protected Map<Integer, BarcodeMatch> barcodeMatches;
		
		/**
		 * A note to propagate for diagnostic file
		 */
		protected String diagnosticNote = "";
		
		public SampleMatch(String sample, Map<Integer, BarcodeMatch> barcodeMatches){
			this.sample= sample;
			this.barcodeMatches = barcodeMatches;
		}

		public SampleMatch(String sample, Map<Integer, BarcodeMatch> barcodeMatches, String diagnosticNote){
			this.sample= sample;
			this.barcodeMatches = barcodeMatches;
			this.diagnosticNote = diagnosticNote;
		}


		/**
		 * @return the sample
		 */
		public String getSample() {
			return sample;
		}


		/**
		 * @return the barcodeMatches
		 */
		public Map<Integer, BarcodeMatch> getBarcodeMatches() {
			return barcodeMatches;
		}


		/**
		 * @return the diagnosticNote
		 */
		public String getDiagnosticNote() {
			return diagnosticNote;
		}


		/**
		 * @param diagnosticNote the diagnosticNote to set
		 */
		public void setDiagnosticNote(String diagnosticNote) {
			this.diagnosticNote = diagnosticNote;
		}


		/**
		 * @param sample the sample to set
		 */
		public void setSample(String sample) {
			this.sample = sample;
		}


		/**
		 * @param barcodeMatches the barcodeMatches to set
		 */
		public void setBarcodeMatches(Map<Integer, BarcodeMatch> barcodeMatches) {
			this.barcodeMatches = barcodeMatches;
		}
		
	
}
