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
package org.embl.gbcs.je.demultiplexer;

import java.util.regex.Pattern;

import org.junit.Assert;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class BarcodeFileGeneralParserTest {

	private static Logger log = LoggerFactory.getLogger(BarcodeFileGeneralParserTest.class);

	private static final String TAB = "\t";
	
	@Test
	public final void testHeaderRegex() {
		
		String [] lines = new String []{
			"SAMPLE" + TAB + "BARCODE1",
			"SAMPLE" + TAB + "BARCODE1"+ TAB + "BARCODE2",
			"SAMPLE" + TAB + "BARCODE1"+ TAB + "BARCODE2" + TAB + "OUT1",
			"SAMPLE" + TAB + "BARCODE1"+ TAB + "BARCODE2" + TAB + "OUT1" + TAB + "OUT2" + TAB + "OUT3",
		};
		
		for(String line : lines){
			Assert.assertTrue(BarcodeFileGeneralParser.isValidHeaderLine(line));
		}
	}
	
	public final void testBarcodeContainsN(){
		Assert.assertTrue(BarcodeFileGeneralParser.containsN("ATGCANCAC"));
		Assert.assertFalse(BarcodeFileGeneralParser.containsN("ATGCATCAC"));
	}
	
	
}
