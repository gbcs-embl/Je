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


import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.StringUtil;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.embl.gbcs.je.jemultiplexer.BarcodeMatch;
import org.junit.Assert;
import org.junit.Test;

public class JemultiplexerTest {

	
	//this is a method to play with quality int / byte / strings
	public final void testQuality(){
		
		System.out.println("0xff is " + 0xff);
		
		String s = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^";
		byte [] arr = StringUtil.stringToBytes(s);
		for (int i = 0; i < s.length(); i++) {
			System.out.print(s.charAt(i)+"\t");
		}
		System.out.println();
		for (byte b : arr) {
			System.out.print(b+"\t");
		}
		System.out.println();
		for (byte b : arr) {
			System.out.print((b & 0xff)+"\t");
		}
		SAMUtils.fastqToPhred(arr);
		System.out.println();
		for (byte b : arr) {
			System.out.print(b+"\t");
		}
		System.out.println();
		for (final byte qual : arr) {
            final int uQual = qual & 0xff;
            System.out.print(uQual+"\t");
		}
		System.out.println();
		 
	}
	
	@Test
	public final void testTrim(){
		int bcl=7;
		int ztrim=2;
		
		String seq  = "CCCCCCATTTTTTTTTTTTTTTTTTTTTTTTTTGG";
		Assert.assertEquals("TTTTTTTTTTTTTTTTTTTTTTTTTT", seq.subSequence(bcl, seq.length()-ztrim));
	}
	
	@Test
	public final void testWrite(){
		Jemultiplexer j = new Jemultiplexer();
		
		/*
		 * first check that if CLIP_BARCODE is false, nothing happens to the read no matter method options
		 */
		j.CLIP_BARCODE = false;
		
		String h    = "@D3FCO8P1:178:C1WLBACXX:7:1101:4883:26830 1:N:0:";
		String seq  = "GATGCTTCGCGCGTGCAGCCCAGGACATCTAAGGGCATCACAGACCTGTTATTGCTCAATCTCATTATTGCTAGACGCAATTTGTCCATTTAAGAAGCTAG";
		String qual = "BCCFDFFFHHHHGBHIGHIJIJJIIJIIIIIIJGIJIJJGIIGFHIGGHGHHHHEHHEDFFFFEEDEEEEDA>>CCBDDDDCDDECDDEEEEEDDCDDDDA";
		final FastqRecord in = new FastqRecord(h, seq, "+", qual);
		
		Assert.assertTrue( areReadsEqual(in, j.write(in, 6, 1, 0, null, "TCTCTC", false, false)) );
		Assert.assertTrue( areReadsEqual(in, j.write(in, 8, 1, 0,null, "TCTCTC", false, false)) );
		Assert.assertTrue( areReadsEqual(in, j.write(in, 8, 10, 0,null, "TCTC", false, false)) );
		
		FastqRecord expected = null;
		expected = new FastqRecord(h+"TCTCTC", seq, "+", qual);
		Assert.assertTrue( areReadsEqual(expected, j.write(in, 6, 1, 0,null, "TCTCTC", true, false)) );
		
		
		//clip BC now
		j.CLIP_BARCODE = true;
		expected = new FastqRecord(
				h+"TCTCTC", 
				"TCGCGCGTGCAGCCCAGGACATCTAAGGGCATCACAGACCTGTTATTGCTCAATCTCATTATTGCTAGACGCAATTTGTCCATTTAAGAAGCTAG", 
				"+", 
				"FFHHHHGBHIGHIJIJJIIJIIIIIIJGIJIJJGIIGFHIGGHGHHHHEHHEDFFFFEEDEEEEDA>>CCBDDDDCDDECDDEEEEEDDCDDDDA");
		Assert.assertTrue( areReadsEqual(expected, j.write(in, 6, 0, 0,null, "TCTCTC", true, false)) );
		expected = new FastqRecord(
				h, 
				"TCGCGCGTGCAGCCCAGGACATCTAAGGGCATCACAGACCTGTTATTGCTCAATCTCATTATTGCTAGACGCAATTTGTCCATTTAAGAAGCTAG", 
				"+", 
				"FFHHHHGBHIGHIJIJJIIJIIIIIIJGIJIJJGIIGFHIGGHGHHHHEHHEDFFFFEEDEEEEDA>>CCBDDDDCDDECDDEEEEEDDCDDDDA");
		Assert.assertTrue( areReadsEqual(expected, j.write(in, 6, 0, 0,null, "TCTCTC", false, false)) );
		
		expected = new FastqRecord(
				h+"TCTCTC", 
				"CGCGCGTGCAGCCCAGGACATCTAAGGGCATCACAGACCTGTTATTGCTCAATCTCATTATTGCTAGACGCAATTTGTCCATTTAAGAAGCTAG", 
				"+", 
				"FHHHHGBHIGHIJIJJIIJIIIIIIJGIJIJJGIIGFHIGGHGHHHHEHHEDFFFFEEDEEEEDA>>CCBDDDDCDDECDDEEEEEDDCDDDDA");
		Assert.assertTrue( areReadsEqual(expected, j.write(in, 6, 1, 0,null, "TCTCTC", true, false)) );
		
		expected = new FastqRecord(
				h, 
				"CGCGCGTGCAGCCCAGGACATCTAAGGGCATCACAGACCTGTTATTGCTCAATCTCATTATTGCTAGACGCAATTTGTCCATTTAAGAAGCTAG", 
				"+", 
				"FHHHHGBHIGHIJIJJIIJIIIIIIJGIJIJJGIIGFHIGGHGHHHHEHHEDFFFFEEDEEEEDA>>CCBDDDDCDDECDDEEEEEDDCDDDDA");
		Assert.assertTrue( areReadsEqual(expected, j.write(in, 6, 1, 0,null, "TCTCTC", false, false)) );
		
		
		
		
	}
	
	
	
	private boolean areReadsEqual(FastqRecord in, FastqRecord out) {
		return in.getReadHeader().equals(out.getReadHeader()) 
				&& in.getBaseQualityHeader().equals(out.getBaseQualityHeader())
				&& in.getBaseQualityHeader().equals(out.getBaseQualityHeader())
				&& in.getReadString().equals(out.getReadString());
	}


	@Test
	public final void testWriteRCHAR(){
		Jemultiplexer j = new Jemultiplexer();
		j.CLIP_BARCODE = true;
		j.READ_NAME_REPLACE_CHAR = ":";
		
		String h    = "@D3FCO8P1:178:C1WLBACXX:7:1101:4883:26830:1:N:0:";
		String seq  = "GATGCTTCGCGCGTGCAGCCCAGGACATCTAAGGGCATCACAGACCTGTTATTGCTCAATCTCATTATTGCTAGACGCAATTTGTCCATTTAAGAAGCTAG";
		String qual = "BCCFDFFFHHHHGBHIGHIJIJJIIJIIIIIIJGIJIJJGIIGFHIGGHGHHHHEHHEDFFFFEEDEEEEDA>>CCBDDDDCDDECDDEEEEEDDCDDDDA";
		final FastqRecord in = new FastqRecord(h, seq, "+", qual);
		
		//clip BC now
		FastqRecord expected = new FastqRecord(
				h+"TCTCTC", 
				"TCGCGCGTGCAGCCCAGGACATCTAAGGGCATCACAGACCTGTTATTGCTCAATCTCATTATTGCTAGACGCAATTTGTCCATTTAAGAAGCTAG", 
				"+", 
				"FFHHHHGBHIGHIJIJJIIJIIIIIIJGIJIJJGIIGFHIGGHGHHHHEHHEDFFFFEEDEEEEDA>>CCBDDDDCDDECDDEEEEEDDCDDDDA");
		Assert.assertTrue( areReadsEqual(expected, j.write(in, 6, 0, 0,null, "TCTCTC", true, false)) );


		j.CLIP_BARCODE = true;
		j.READ_NAME_REPLACE_CHAR = "_";
		h    = "@D3FCO8P1:178:C1WLBACXX:7:1101:4883:26830_1:N:0:";
		seq  = "GATGCTTCGCGCGTGCAGCCCAGGACATCTAAGGGCATCACAGACCTGTTATTGCTCAATCTCATTATTGCTAGACGCAATTTGTCCATTTAAGAAGCTAG";
		qual = "BCCFDFFFHHHHGBHIGHIJIJJIIJIIIIIIJGIJIJJGIIGFHIGGHGHHHHEHHEDFFFFEEDEEEEDA>>CCBDDDDCDDECDDEEEEEDDCDDDDA";
		final FastqRecord in2 = new FastqRecord(h, seq, "+", qual);

		expected = new FastqRecord(
				h+"TCTCTC", 
				"TCGCGCGTGCAGCCCAGGACATCTAAGGGCATCACAGACCTGTTATTGCTCAATCTCATTATTGCTAGACGCAATTTGTCCATTTAAGAAGCTAG", 
				"+", 
				"FFHHHHGBHIGHIJIJJIIJIIIIIIJGIJIJJGIIGFHIGGHGHHHHEHHEDFFFFEEDEEEEDA>>CCBDDDDCDDECDDEEEEEDDCDDDDA");
		Assert.assertTrue( areReadsEqual(expected, j.write(in2, 6, 0, 0,null, "TCTCTC", true, false)) );
		
	}
	
	
	@Test
	public final void testfindBarcode(){
		Jemultiplexer j = new Jemultiplexer();
		j.MAX_MISMATCHES_1 = 1;
		j.MAX_MISMATCHES_2 = 1;
		j.MIN_BASE_QUALITY_1 = 10;
		j.MIN_BASE_QUALITY_2 = 10;
		j.MIN_MISMATCH_DELTA_1 = 1;
		j.MIN_MISMATCH_DELTA_2 = 1;
		
		String [] bcArr = {"AAAA", "AAAC", "AACC", "ACCC"};
		j.barcodeValidator = new BarcodeValidator(null, j);
		j.barcodeValidator.barcodeSetRead1 = new TreeSet<String>(Arrays.asList(bcArr));
		j.barcodeValidator.barcodeSetRead2 = new TreeSet<String>(Arrays.asList(bcArr));
		j.barcodeBytes1 = new byte [bcArr.length][];
		j.barcodeBytes2 = new byte [bcArr.length][];
		for (int i = 0; i < bcArr.length; i++) {
			j.barcodeBytes1[i] = StringUtil.stringToBytes(bcArr[i]);
			j.barcodeBytes2[i] = StringUtil.stringToBytes(bcArr[i]);
		}
		
		Map<String, String> bc2sample = new HashMap<String, String>();
		Map<String, Set<String>> sample2bcSet = new HashMap<String, Set<String>>();
		for (int i = 0; i < bcArr.length; i++) {
			bc2sample.put(bcArr[i], "S"+i);
			TreeSet<String> s = new TreeSet<String>();
			s.add(bcArr[i]);
			sample2bcSet.put("S"+i, s);
		}
		
		/*
		 * The difference with the findBestBarcode() is that here the quality is not considered if the subseq matches a barcode
		 */
		
		BarcodeMatch bm = null;
		
		//give exact seq first , result should be always the same no matter the quality
		bm = j.findBarcode(BarcodePosition.READ_1, "AAAA", "HHHH", bc2sample, sample2bcSet);
		Assert.assertTrue(bm.matched);
		Assert.assertEquals("AAAA", bm.barcode);
		Assert.assertEquals(0, bm.mismatches);
		Assert.assertEquals(4, bm.mismatchesToSecondBest);
		
		bm = j.findBarcode(BarcodePosition.READ_1,"AAAA", "#HHH", bc2sample, sample2bcSet);
		Assert.assertTrue(bm.matched);
		Assert.assertEquals("AAAA", bm.barcode);
		Assert.assertEquals(0, bm.mismatches);
		Assert.assertEquals(4, bm.mismatchesToSecondBest);
		
		bm = j.findBarcode(BarcodePosition.READ_1,"AAAA", "####", bc2sample, sample2bcSet);
		Assert.assertTrue(bm.matched);
		Assert.assertEquals("AAAA", bm.barcode);
		Assert.assertEquals(0, bm.mismatches);
		Assert.assertEquals(4, bm.mismatchesToSecondBest);
		
		//start having mismatches
		bm = j.findBarcode(BarcodePosition.READ_1,"CCCC", "HHHH", bc2sample, sample2bcSet);
		Assert.assertTrue(bm.matched);
		Assert.assertEquals("ACCC", bm.barcode);
		Assert.assertEquals(1, bm.mismatches);
		Assert.assertEquals(2, bm.mismatchesToSecondBest);
		
		// mismatch is 2 with both AACC and ACCC due to quality => no match => should be null
		bm = j.findBarcode(BarcodePosition.READ_1,"CCCC", "H#HH", bc2sample, sample2bcSet);
		Assert.assertNull(bm);
		
		
		
		
	}
	
	
	@Test
	public final void testfindBestBarcode(){
		Jemultiplexer j = new Jemultiplexer();
		j.MAX_MISMATCHES_1 = 1;
		j.MAX_MISMATCHES_2 = 1;
		j.MIN_BASE_QUALITY_1 = 10;
		j.MIN_BASE_QUALITY_2 = 10;
		j.MIN_MISMATCH_DELTA_1 = 1;
		j.MIN_MISMATCH_DELTA_2 = 1;
		
		String [] bcArr = {"AAAA", "AAAC", "AACC", "ACCC"};
		j.barcodeBytes1 = new byte [bcArr.length][];
		j.barcodeBytes2 = new byte [bcArr.length][];
		for (int i = 0; i < bcArr.length; i++) {
			j.barcodeBytes1[i] = StringUtil.stringToBytes(bcArr[i]);
			j.barcodeBytes2[i] = StringUtil.stringToBytes(bcArr[i]);
		}
		
		BarcodeMatch bm = null;
		
		bm = j.findBestBarcode(BarcodePosition.READ_1, "AAAA", "HHHH");
		Assert.assertTrue(bm.matched);
		Assert.assertEquals("AAAA", bm.barcode);
		Assert.assertEquals(0, bm.mismatches);
		Assert.assertEquals(1, bm.mismatchesToSecondBest);
		
		bm = j.findBestBarcode(BarcodePosition.READ_1, "AAAA", "#HHH");
		Assert.assertTrue(bm.matched);
		Assert.assertEquals("AAAA", bm.barcode);
		Assert.assertEquals(1, bm.mismatches);
		Assert.assertEquals(2, bm.mismatchesToSecondBest);
		
		bm = j.findBestBarcode(BarcodePosition.READ_1, "CCCC", "HHHH");
		Assert.assertTrue(bm.matched);
		Assert.assertEquals("ACCC", bm.barcode);
		Assert.assertEquals(1, bm.mismatches);
		Assert.assertEquals(2, bm.mismatchesToSecondBest);
		
		bm = j.findBestBarcode(BarcodePosition.READ_1, "CCCC", "H#HH");
		Assert.assertFalse(bm.matched); // mismatch is 2 with both AACC and ACCC due to quality
		Assert.assertNull(bm.barcode);
		
		
		bm = j.findBestBarcode(BarcodePosition.READ_1, "AAAA", "####");
		Assert.assertFalse(bm.matched);
		Assert.assertNull(bm.barcode);
		
		bm = j.findBestBarcode(BarcodePosition.READ_1, "AGGC", "HHHH");
		Assert.assertTrue(!bm.matched);
		Assert.assertNull(bm.barcode);
		
		//2 MM
		j.MAX_MISMATCHES_1 = 2;
		j.MAX_MISMATCHES_2 = 2;
		j.MIN_MISMATCH_DELTA_1 = 1;
		j.MIN_MISMATCH_DELTA_2 = 1;
		bm = j.findBestBarcode(BarcodePosition.READ_1, "GGAA", "HHHH");
		Assert.assertTrue(bm.matched);
		Assert.assertEquals("AAAA", bm.barcode);
		Assert.assertEquals(2, bm.mismatches);
		Assert.assertEquals(3, bm.mismatchesToSecondBest);
		
		j.MAX_MISMATCHES_1 = 2;
		j.MAX_MISMATCHES_2 = 2;
		j.MIN_MISMATCH_DELTA_1 = 2;
		j.MIN_MISMATCH_DELTA_2 = 2;
		bm = j.findBestBarcode(BarcodePosition.READ_1, "GGAA", "HHHH");
		Assert.assertTrue(!bm.matched);
		Assert.assertNull(bm.barcode);
		
		
	}
	
	@Test
	public final void testQualityConversion() {
		
		//first understand how this works...
		int i = 0;
		for(byte b : StringUtil.stringToBytes("!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^")){
			Assert.assertEquals((byte)(i+33), b);
			i++;
		}
		
		//
		Jemultiplexer j = new Jemultiplexer();
		byte[] barr = j.convertQualityStringToByteArray("#1=ADFH+*");
		i=0;
		Assert.assertEquals((byte)2, barr[i++]); //# is 2
		Assert.assertEquals((byte)16, barr[i++]); //1 is 16
		Assert.assertEquals((byte)28, barr[i++]); //= is 28
		Assert.assertEquals((byte)32, barr[i++]); // A is 32
		Assert.assertEquals((byte)35, barr[i++]); // D is 35
		Assert.assertEquals((byte)37, barr[i++]); // F is 37
		Assert.assertEquals((byte)39, barr[i++]); // H is 39
		Assert.assertEquals((byte)10, barr[i++]); // + is 10
		Assert.assertEquals((byte)9, barr[i++]); // * is 9
	}
	
	@Test
	public final void testcountMismatches() {
		
		//
		Jemultiplexer j = new Jemultiplexer();

		//no quality given
		Assert.assertEquals("No Mismatch", 0, j.countMismatches(
				StringUtil.stringToBytes("AAAA"), 
				StringUtil.stringToBytes("AAAA"), null, 10) );
		
		Assert.assertEquals("2 Mismatches", 2, j.countMismatches(
				StringUtil.stringToBytes("AATT"), 
				StringUtil.stringToBytes("AAAA"), null, 10) );
		
		Assert.assertEquals("3 Mismatches", 3, j.countMismatches(
				StringUtil.stringToBytes("CCAC"), 
				StringUtil.stringToBytes("AAAA"), null, 10) );
		
		//all good qualities : only test mismatch counting
		byte[] qualities = j.convertQualityStringToByteArray("FFFF"); //all 37 quality
		Assert.assertEquals("No Mismatch", 0, j.countMismatches(
				StringUtil.stringToBytes("AAAA"), 
				StringUtil.stringToBytes("AAAA"), qualities, 10) );
		
		Assert.assertEquals("2 Mismatches", 2, j.countMismatches(
				StringUtil.stringToBytes("AATT"), 
				StringUtil.stringToBytes("AAAA"), qualities, 10) );
		
		Assert.assertEquals("3 Mismatches", 3, j.countMismatches(
				StringUtil.stringToBytes("CCAC"), 
				StringUtil.stringToBytes("AAAA"), qualities, 10) );
		
		
		//get one bad quality base 
		qualities = j.convertQualityStringToByteArray("FF#F"); //all 37 quality but 1 '#' ie '2' 
		Assert.assertEquals("0 Mismatch but 1 quality issue => 1 mismatch", 1, j.countMismatches(
				StringUtil.stringToBytes("AAAA"), 
				StringUtil.stringToBytes("AAAA"), qualities, 10) );
		
		Assert.assertEquals("2 Mismatches", 2, j.countMismatches(
				StringUtil.stringToBytes("AATT"), 
				StringUtil.stringToBytes("AAAA"), qualities, 10) );
		
		Assert.assertEquals("3 Mismatches and 1 quality mismatch => 3 mismatches", 4, j.countMismatches(
				StringUtil.stringToBytes("CCAC"), 
				StringUtil.stringToBytes("AAAA"), qualities, 10) );
		
		//get one bad quality base just below threshold
		
		qualities = j.convertQualityStringToByteArray("FF*F");  
		Assert.assertEquals("0 Mismatch but 1 quality issue => 1 mismatch", 1, j.countMismatches(
				StringUtil.stringToBytes("AAAA"), 
				StringUtil.stringToBytes("AAAA"), qualities, 10) );

		Assert.assertEquals("2 Mismatches", 2, j.countMismatches(
				StringUtil.stringToBytes("AATT"), 
				StringUtil.stringToBytes("AAAA"), qualities, 10) );

		Assert.assertEquals("3 Mismatches and 1 quality mismatch => 3 mismatches", 4, j.countMismatches(
				StringUtil.stringToBytes("CCAC"), 
				StringUtil.stringToBytes("AAAA"), qualities, 10) );
				
		
		qualities = j.convertQualityStringToByteArray("FF+F");  
		Assert.assertEquals("0 Mismatch", 0, j.countMismatches(
				StringUtil.stringToBytes("AAAA"), 
				StringUtil.stringToBytes("AAAA"), qualities, 10) );

		Assert.assertEquals("2 Mismatches", 2, j.countMismatches(
				StringUtil.stringToBytes("AATT"), 
				StringUtil.stringToBytes("AAAA"), qualities, 10) );

		Assert.assertEquals("3 Mismatches ", 3, j.countMismatches(
				StringUtil.stringToBytes("CCAC"), 
				StringUtil.stringToBytes("AAAA"), qualities, 10) );
				
		
		
		//all crapy qualities
		qualities = j.convertQualityStringToByteArray("####");  
		Assert.assertEquals(4, j.countMismatches(
				StringUtil.stringToBytes("AAAA"), 
				StringUtil.stringToBytes("AAAA"), qualities, 10) );
		Assert.assertEquals(4, j.countMismatches(
				StringUtil.stringToBytes("AATT"), 
				StringUtil.stringToBytes("AAAA"), qualities, 10) );
		Assert.assertEquals(4, j.countMismatches(
				StringUtil.stringToBytes("CCAC"), 
				StringUtil.stringToBytes("AAAA"), qualities, 10) );
		
		
	}
	
	

	

}














