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

import java.util.List;

import org.embl.cg.utilitytools.utils.StringUtil;
import org.junit.Assert;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class ReadLayoutTest {
	private static Logger log = LoggerFactory.getLogger(ReadLayoutTest.class);

	
	
	
	/**
	 * Give a wrong layout and make sure we get an exception 
	 *  
	 */
	@Test
	public final void testIncorrectLayouts() {
		String [] layouts = {
				"<BARCODE:-6>NN<UMI:8>N<SAMPLE:20>",
				"<BARCODE:3>NN<UMI:x>N<SAMPLE:20>",
				"<BARCODE:>NN<UMI:8>N<SAMPLE:20>",
				"ATGCATGC",
				"<BARCOD:3>NN<UMI:8>N<SAMPLE:20>",
				"<BARCODE:3>NN<UMI:8>N<SAMPLE:20><WRONG:9>",
		};
		for (String layout : layouts) {
			try {
				ReadLayout l = new ReadLayout(layout);
				Assert.fail("Should have thrown an Jexception for layout "+layout);
			} catch (Jexception e) {
				// ok
			}
		}
	}
	
	/**
	 * Give a BC-UMI-SAMPLE layout with exact sample lenght and
	 *  make sure the correct pieces are extracted
	 */
	@Test
	public final void testMultiBlockLayout2() {
		String layout = "A<UMI:3><BARCODE:6><BARCODE:3><UMI:8>N<SAMPLE:x>";
		String bc1 = "ACACAC";
		String bc2 = "GGG";
		String umi1 = "ACT";
		String umi2 = "ACGTACGT";
		String smpl = "TACGACNACN"+"NACACACAGT"+"CG"; //22 long
		String read = "A"+umi1+bc1+bc2+umi2+"A"+smpl;
		
		ReadLayout l = new ReadLayout(layout);
		
		Assert.assertTrue("containsBarcode() should be true but was not for layout: "+layout, l.containsBarcode() );
		Assert.assertTrue("containsUMI() should be true but was not for layout: "+layout, l.containsUMI() );
		Assert.assertTrue("containsSampleSequence() should be true but was not for layout: "+layout, l.containsSampleSequence() );
		
		Assert.assertEquals(2, l.umiBlockNumber());
		Assert.assertEquals(2, l.barcodeBlockNumber());
		
		String[] arr  = l.extractBarcodes(read);
		log.debug(arr[0]);
		log.debug(arr[1]);
		Assert.assertEquals(bc1, arr[0]);
		Assert.assertEquals(bc2, arr[1]);
		
		Assert.assertEquals(umi1, l.extractUMIs(read)[0]);
		Assert.assertEquals(umi2, l.extractUMIs(read)[1]);
		
		Assert.assertEquals(bc1+bc2, l.extractBarcode(read));
		Assert.assertEquals(umi1+umi2, l.extractUMI(read));
		Assert.assertEquals(smpl, l.extractSampleSequence(read));	
	}
	
	/**
	 * Give a BC-UMI-SAMPLE layout with exact sample lenght and
	 *  make sure the correct pieces are extracted
	 */
	@Test
	public final void testMultiBlockLayout() {
		String layout = "A<UMI:3><BARCODE:6>NN<BARCODE:3><UMI:8>N<SAMPLE:x>";
		String bc1 = "ACACAC";
		String bc2 = "GGG";
		String umi1 = "ACT";
		String umi2 = "ACGTACGT";
		String smpl = "TACGACNACN"+"NACACACAGT"+"CG"; //22 long
		String read = "A"+umi1+bc1+"TT"+bc2+umi2+"A"+smpl;
		
		ReadLayout l = new ReadLayout(layout);
		
		Assert.assertTrue("containsBarcode() should be true but was not for layout: "+layout, l.containsBarcode() );
		Assert.assertTrue("containsUMI() should be true but was not for layout: "+layout, l.containsUMI() );
		Assert.assertTrue("containsSampleSequence() should be true but was not for layout: "+layout, l.containsSampleSequence() );
		
		Assert.assertEquals(2, l.umiBlockNumber());
		Assert.assertEquals(2, l.barcodeBlockNumber());
		
		String[] arr  = l.extractBarcodes(read);
		log.debug(arr[0]);
		log.debug(arr[1]);
		Assert.assertEquals(bc1, arr[0]);
		Assert.assertEquals(bc2, arr[1]);
		
		Assert.assertEquals(umi1, l.extractUMIs(read)[0]);
		Assert.assertEquals(umi2, l.extractUMIs(read)[1]);
		
		Assert.assertEquals(bc1+bc2, l.extractBarcode(read));
		Assert.assertEquals(umi1+umi2, l.extractUMI(read));
		Assert.assertEquals(smpl, l.extractSampleSequence(read));	
	}
	
	
	/**
	 * Give a BC-UMI-SAMPLE layout with exact sample lenght and
	 *  make sure the correct pieces are extracted
	 */
	@Test
	public final void testCorrectBcUmiSampleLayoutSimple2() {
		String layout = "<BARCODE:6>NN<UMI:8>N<SAMPLE:20>";
		String bc = "ACACAC";
		String umi = "ACGTACGT";
		String smpl = "TACGACNACN"+"NACACACAGT"+"CG"; //22 long
		String read = bc+"TT"+umi+"A"+smpl;
		
		ReadLayout l = new ReadLayout(layout);
		
		Assert.assertTrue("containsBarcode() should be true but was not for layout: "+layout, l.containsBarcode() );
		Assert.assertTrue("containsUMI() should be true but was not for layout: "+layout, l.containsUMI() );
		Assert.assertTrue("containsSampleSequence() should be true but was not for layout: "+layout, l.containsSampleSequence() );
		
		Assert.assertEquals(1, l.umiBlockNumber());
		Assert.assertEquals(1, l.barcodeBlockNumber());
		
		Assert.assertEquals(bc, l.extractBarcode(read));
		Assert.assertEquals(umi, l.extractUMI(read));
		Assert.assertEquals("TACGACNACN"+"NACACACAGT", l.extractSampleSequence(read));	
	}
	
	/**
	 * Give a BC-UMI-SAMPLE layout with relative sample lenght and
	 *  make sure the correct pieces are extracted
	 */
	@Test
	public final void testCorrectBcUmiSampleLayoutSimple3() {
		String layout = "<BARCODE:6>NN<UMI:8>N<SAMPLE:-2>";
		String bc = "ACACAC";
		String umi = "ACGTACGT";
		String smpl = "TACGACNACN"+"NACACACAGT"+"CG"; //22 long
		String read = bc+"TT"+umi+"A"+smpl;
		
		ReadLayout l = new ReadLayout(layout);
		
		Assert.assertTrue("containsBarcode() should be true but was not for layout: "+layout, l.containsBarcode() );
		Assert.assertTrue("containsUMI() should be true but was not for layout: "+layout, l.containsUMI() );
		Assert.assertTrue("containsSampleSequence() should be true but was not for layout: "+layout, l.containsSampleSequence() );
		
		Assert.assertEquals(1, l.umiBlockNumber());
		Assert.assertEquals(1, l.barcodeBlockNumber());
		
		Assert.assertEquals(bc, l.extractBarcode(read));
		Assert.assertEquals(umi, l.extractUMI(read));
		Assert.assertEquals("TACGACNACN"+"NACACACAGT", l.extractSampleSequence(read));	
	}
	
	
	/**
	 * Give a BC-UMI layout and make sure the correct pieces are extracted
	 */
	@Test
	public final void testCorrectBcUmiLayoutSimple() {
		String layout = "<BARCODE:6>NN<UMI:8>NN";
		String bc = "ACACAC";
		String umi = "ACGTACGT";
		
		String read = bc+"TT"+umi+"AA";
		
		ReadLayout l = new ReadLayout(layout);
		
		Assert.assertTrue("containsBarcode() should be true but was not for layout: "+layout, l.containsBarcode() );
		Assert.assertTrue("containsUMI() should be true but was not for layout: "+layout, l.containsUMI() );
		Assert.assertFalse("containsSampleSequence() should be false but was not for layout: "+layout, l.containsSampleSequence() );
		
		Assert.assertEquals(1, l.umiBlockNumber());
		Assert.assertEquals(1, l.barcodeBlockNumber());
		
		Assert.assertEquals(bc, l.extractBarcode(read));
		Assert.assertEquals(umi, l.extractUMI(read));
		Assert.assertNull(l.extractSampleSequence(read));	
	}
	
	
	/**
	 * Give a BC-UMI-SAMPLE layout and make sure the correct pieces are extracted
	 */
	@Test
	public final void testCorrectBcUmiSampleLayoutSimple() {
		String layout = "<BARCODE:6>NN<UMI:8>N<SAMPLE:x>";
		String bc = "ACACAC";
		String umi = "ACGTACGT";
		String smpl = "TACGACNACN"+"NACACACAGT"+"CG"; //22 long
		String read = bc+"TT"+umi+"A"+smpl;
		
		ReadLayout l = new ReadLayout(layout);
		
		Assert.assertTrue("containsBarcode() should be true but was not for layout: "+layout, l.containsBarcode() );
		Assert.assertTrue("containsUMI() should be true but was not for layout: "+layout, l.containsUMI() );
		Assert.assertTrue("containsSampleSequence() should be true but was not for layout: "+layout, l.containsSampleSequence() );
		
		Assert.assertEquals(1, l.umiBlockNumber());
		Assert.assertEquals(1, l.barcodeBlockNumber());
		
		Assert.assertEquals(bc, l.extractBarcode(read));
		Assert.assertEquals(umi, l.extractUMI(read));
		Assert.assertEquals(smpl, l.extractSampleSequence(read));	
	}
	
	
	/**
	 * Give a BC-SAMPLE layout and make sure the correct pieces are extracted
	 */
	@Test
	public final void testCorrectBcSampleLayoutSimple() {
		String layout = "<BARCODE:6>N<SAMPLE:x>";
		String bc = "ACACAC";
		String smpl = "TACGACNACN"+"NACACACAGT"+"CG"; //22 long
		String read = bc+"A"+smpl;
		
		ReadLayout l = new ReadLayout(layout);
		
		Assert.assertTrue("containsBarcode() should be true but was not for layout: "+layout, l.containsBarcode() );
		Assert.assertFalse("containsUMI() should be false but was not for layout: "+layout, l.containsUMI() );
		Assert.assertTrue("containsSampleSequence() should be true but was not for layout: "+layout, l.containsSampleSequence() );
		
		Assert.assertEquals(0, l.umiBlockNumber());
		Assert.assertEquals(1, l.barcodeBlockNumber());
		
		Assert.assertEquals(bc, l.extractBarcode(read));
		Assert.assertNull(l.extractUMI(read));
		Assert.assertEquals(smpl, l.extractSampleSequence(read));	
	}
	
	/**
	 * Give a UMI-SAMPLE layout and make sure the correct pieces are extracted
	 */
	@Test
	public final void testCorrectUmiSampleLayoutSimple() {
		String layout = "NN<UMI:8>N<SAMPLE:x>";
		String umi = "ACGTACGT";
		String smpl = "TACGACNACN"+"NACACACAGT"+"CG"; //22 long
		String read = "TT"+umi+"A"+smpl;
		
		ReadLayout l = new ReadLayout(layout);
		
		Assert.assertFalse("containsBarcode() should be false but was not for layout: "+layout, l.containsBarcode() );
		Assert.assertTrue("containsUMI() should be true but was not for layout: "+layout, l.containsUMI() );
		Assert.assertTrue("containsSampleSequence() should be true but was not for layout: "+layout, l.containsSampleSequence() );
		
		Assert.assertEquals(1, l.umiBlockNumber());
		Assert.assertEquals(0, l.barcodeBlockNumber());
		
		Assert.assertNull( l.extractBarcode(read));
		Assert.assertEquals(umi, l.extractUMI(read));
		Assert.assertEquals(smpl, l.extractSampleSequence(read));	
	}
	
	
	@Test
	public final void testCorrectUmiOnlyLayoutSimple() {
		String layout = "<UMI:8>";
		String umi = "ACGTACGT";
		String read = umi;
		
		ReadLayout l = new ReadLayout(layout);
		
		Assert.assertFalse("containsBarcode() should be false but was not for layout: "+layout, l.containsBarcode() );
		Assert.assertTrue("containsUMI() should be true but was not for layout: "+layout, l.containsUMI() );
		Assert.assertFalse("containsSampleSequence() should be false but was not for layout: "+layout, l.containsSampleSequence() );
		
		Assert.assertEquals(1, l.umiBlockNumber());
		Assert.assertEquals(0, l.barcodeBlockNumber());
		
		Assert.assertNull( l.extractBarcode(read));
		Assert.assertEquals(umi, l.extractUMI(read));
		Assert.assertNull(l.extractSampleSequence(read));	
	}
	
	@Test
	public final void testCorrectSampleOnlyLayout() {
		String layout = "<SAMPLE:x>";
		String spl = "ACGTACGT";
		String read = spl;
		
		ReadLayout l = new ReadLayout(layout);
		
		Assert.assertFalse("containsBarcode() should be false but was not for layout: "+layout, l.containsBarcode() );
		Assert.assertFalse("containsUMI() should be false but was not for layout: "+layout, l.containsUMI() );
		Assert.assertTrue("containsSampleSequence() should be true but was not for layout: "+layout, l.containsSampleSequence() );
		
		Assert.assertEquals(0, l.umiBlockNumber());
		Assert.assertEquals(0, l.barcodeBlockNumber());
		
		Assert.assertNull( l.extractBarcode(read));
		Assert.assertEquals(spl, l.extractSampleSequence(read));
		Assert.assertNull(l.extractUMI(read));	
	}
	
	@Test
	public final void testCorrectSampleOnlyLayout2() {
		String layout = "<SAMPLE:-2>";
		String spl = "ACGTAC";
		String read = spl+"GT";
		
		ReadLayout l = new ReadLayout(layout);
		
		Assert.assertFalse("containsBarcode() should be false but was not for layout: "+layout, l.containsBarcode() );
		Assert.assertFalse("containsUMI() should be false but was not for layout: "+layout, l.containsUMI() );
		Assert.assertTrue("containsSampleSequence() should be true but was not for layout: "+layout, l.containsSampleSequence() );
		
		Assert.assertEquals(0, l.umiBlockNumber());
		Assert.assertEquals(0, l.barcodeBlockNumber());
		
		Assert.assertNull( l.extractBarcode(read));
		Assert.assertEquals(spl, l.extractSampleSequence(read));
		Assert.assertNull(l.extractUMI(read));	
	}
	
	@Test
	public final void testCorrectBarcodeOnlyLayout() {
		String layout = "<BARCODE:6>";
		String bc = "ACGTAC";
		String read = bc + "GT";
		
		ReadLayout l = new ReadLayout(layout);
		
		Assert.assertTrue("containsBarcode() should be true but was not for layout: "+layout, l.containsBarcode() );
		Assert.assertFalse("containsUMI() should be false but was not for layout: "+layout, l.containsUMI() );
		Assert.assertFalse("containsSampleSequence() should be false but was not for layout: "+layout, l.containsSampleSequence() );
		
		Assert.assertEquals(0, l.umiBlockNumber());
		Assert.assertEquals(1, l.barcodeBlockNumber());
		
		Assert.assertEquals( bc, l.extractBarcode(read));
		Assert.assertNull(l.extractSampleSequence(read));
		Assert.assertNull(l.extractUMI(read));	
	}
	
	


	private String printErrors(List<String> errors) {
		return StringUtil.mergeList(errors, "\n");
	}

	
	
	
}
