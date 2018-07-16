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

import org.embl.cg.utilitytools.utils.ExceptionUtil;
import org.embl.cg.utilitytools.utils.StringUtil;
import org.junit.Assert;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.FastqQualityFormat;

public class FastqWriterLayoutTest {
	private static Logger log = LoggerFactory.getLogger(FastqWriterLayoutTest.class);

	
	
	
	/**
	 * Give a wrong layout and make sure we get an exception 
	 *  
	 */
	@Test
	public final void testIncorrectLayouts() {
		
		ReadLayout[] rls = new ReadLayout[]{
				new ReadLayout("<BARCODE:6><SAMPLE:-2>"),
				new ReadLayout("<UMI:10>")
		};
		
		String [] layouts = {
				"<BARCODE><UMI1><SAMPLE2>", // missign idx on BARCODE
				"B1G1S1H1", //wrong letters
				"<BARCODE1><UMI13><SAMPLE24>"//wrong ids
		};
		for (String layout : layouts) {
			try {
				FastqWriterLayout l = new FastqWriterLayout(layout, null, rls, false, FastqWriterLayout.DEFAULT_READNAME_DELIMITOR, FastqQualityFormat.Standard);
				Assert.fail("Should have thrown an Jexception for layout "+layout);
			} catch (Jexception e) {
				log.debug(e.getMessage());
				// ok
			}
		}
	}
	
	/**
	 *  
	 */
	@Test
	public final void testStandardSE() {
		
		String readname = "fakename";
		String qualityheader = "+";
		String bcSeq =  "AAATTT";
		String bcQual = "AAAAAE";
		String samplSeq =  "CGATGACTGCTAGCTGCTAGCTAGCAT";
		String samplQual = "EE#EEEEEEEEEEEEEEEEEEEEEEEE";
		
		String seqlayout = "<SAMPLE1>" ;
		String headlayout = "<READBAR1>" ;
		FastqRecord ori = new FastqRecord(
				readname, 
				bcSeq+samplSeq, 
				qualityheader, 
				bcQual+samplQual);
		try {
			FastqWriterLayout fwl = new FastqWriterLayout(seqlayout, headlayout, new ReadLayout[] {new ReadLayout("<BARCODE:6><SAMPLE:x>")}, false, FastqWriterLayout.DEFAULT_READNAME_DELIMITOR, FastqQualityFormat.Standard);		
			FastqRecord r = fwl.assembleRecord(ori, null);
			
			Assert.assertEquals(readname+":"+bcSeq, r.getReadName());
			Assert.assertEquals(qualityheader, r.getBaseQualityHeader());
			Assert.assertEquals(samplSeq, r.getReadString());
			Assert.assertEquals(samplQual, r.getBaseQualityString());
			
			seqlayout = "<READBAR1><SAMPLE1>" ; //write the real barcode sequence
			fwl = new FastqWriterLayout(seqlayout, headlayout, new ReadLayout[] {new ReadLayout("<BARCODE:6><SAMPLE:x>")}, false, FastqWriterLayout.DEFAULT_READNAME_DELIMITOR, FastqQualityFormat.Standard);		
			r = fwl.assembleRecord(ori, null);
			Assert.assertEquals(bcSeq+samplSeq, r.getReadString());
			Assert.assertEquals(bcQual+samplQual, r.getBaseQualityString());
			Assert.assertEquals(ori.getReadString(), r.getReadString());
			Assert.assertEquals(ori.getBaseQualityString(), r.getBaseQualityString());
			
		} catch (Jexception e) {
			log.error(ExceptionUtil.getStackTrace(e));
			Assert.fail("Got error with message : "+e.getMessage());
		}
		
	}
	
	
	/**
	 *  
	 */
	@Test
	public final void testStandardSEWithQualityInReadName() {
		
		String readname = "fakename";
		String qualityheader = "+";
		String bcSeq =  "AAATTT";
		String bcQual = "AAAAAE"; //A is 32 + 33 and E is 36 + 33 
		String bcQualNum = "323232323236";
		String samplSeq =  "CGATGACTGCTAGCTGCTAGCTAGCAT";
		String samplQual = "EE#EEEEEEEEEEEEEEEEEEEEEEEE";
		
		String seqlayout = "<SAMPLE1>" ;
		String headlayout = "<READBAR1>" ;
		FastqRecord ori = new FastqRecord(
				readname, 
				bcSeq+samplSeq, 
				qualityheader, 
				bcQual+samplQual);
		try {
			FastqWriterLayout fwl = new FastqWriterLayout(seqlayout, headlayout, new ReadLayout[] {new ReadLayout("<BARCODE:6><SAMPLE:x>")}, true, FastqWriterLayout.DEFAULT_READNAME_DELIMITOR, FastqQualityFormat.Standard);		
			FastqRecord r = fwl.assembleRecord(ori, null);
			
			Assert.assertEquals(readname+":"+bcSeq+bcQualNum, r.getReadName());
			Assert.assertEquals(qualityheader, r.getBaseQualityHeader());
			Assert.assertEquals(samplSeq, r.getReadString());
			Assert.assertEquals(samplQual, r.getBaseQualityString());
			
			seqlayout = "<READBAR1><SAMPLE1>" ; //write the real barcode sequence
			fwl = new FastqWriterLayout(seqlayout, headlayout, new ReadLayout[] {new ReadLayout("<BARCODE:6><SAMPLE:x>")}, true, FastqWriterLayout.DEFAULT_READNAME_DELIMITOR, FastqQualityFormat.Standard);		
			r = fwl.assembleRecord(ori, null);
			Assert.assertEquals(readname+":"+bcSeq+bcQualNum, r.getReadName());
			Assert.assertEquals(bcSeq+samplSeq, r.getReadString());
			Assert.assertEquals(bcQual+samplQual, r.getBaseQualityString());
			Assert.assertEquals(ori.getReadString(), r.getReadString());
			Assert.assertEquals(ori.getBaseQualityString(), r.getBaseQualityString());
			
		} catch (Jexception e) {
			log.error(ExceptionUtil.getStackTrace(e));
			Assert.fail("Got error with message : "+e.getMessage());
		}
		
	}
	
	
	
	/**
	 *  
	 */
	@Test
	public final void testStandardPE() {
		
		String readname = "NB501764:229:HV253BGX2:1:11101:6294:1041";
		String qualityheader = "+";
		String r1 = "CTGAGTACGTCAACACGACCAGTCAAACTTGACTGAAGGGATCTTCTTACATTTCCCTTCTTCAACATATTCTCTAATAT";
		String q1 = "AAAAA#EE//EEEEEEE/AEEAAEEEEEEAEEEA/EEAEEEEAEEE/EEEAEEEEEEEEEEAEEAEEEEAEEAEEAE/EA";
		String r2 = "CTGAGTGTTCGANCAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
		String q2 = "AAAAAEEA/AEE#EEE################################################################";
		FastqRecord fwd = new FastqRecord(readname, r1, qualityheader, q1);
		FastqRecord rev = new FastqRecord(readname, r2, qualityheader, q2);
		
		//consider similar barcodes
		ReadLayout rl1 = new ReadLayout("<BARCODE:6><SAMPLE1:x>");
		ReadLayout rl2 = new ReadLayout("<BARCODE:6><SAMPLE2:x>");
		
		
		try {
			FastqWriterLayout fwdl = new FastqWriterLayout("<SAMPLE1>", "<READBAR1>", new ReadLayout[]{rl1, rl2}, false, FastqWriterLayout.DEFAULT_READNAME_DELIMITOR, FastqQualityFormat.Standard);
			FastqWriterLayout revl = new FastqWriterLayout("<SAMPLE2>", "<READBAR1>", new ReadLayout[]{rl1, rl2}, false, FastqWriterLayout.DEFAULT_READNAME_DELIMITOR, FastqQualityFormat.Standard);
			
			FastqRecord nfwd = fwdl.assembleRecord(new FastqRecord[]{fwd, rev},  null);
			FastqRecord nrev = revl.assembleRecord(new FastqRecord[]{fwd, rev},  null);
			
			Assert.assertEquals(readname+":CTGAGT", nfwd.getReadName());
			Assert.assertEquals(readname+":CTGAGT", nrev.getReadName());
			
			
			Assert.assertEquals(r1.substring(6), nfwd.getReadString());
			Assert.assertEquals(r2.substring(6), nrev.getReadString());
			
			Assert.assertEquals(q1.substring(6), nfwd.getBaseQualityString());
			Assert.assertEquals(q2.substring(6), nrev.getBaseQualityString());
			
		} catch (Jexception e) {
			log.error(ExceptionUtil.getStackTrace(e));
			Assert.fail("Got error with message : "+e.getMessage());
		}
		
		//consider different barcodes
		rl1 = new ReadLayout("<BARCODE1:6><SAMPLE1:x>");
		rl2 = new ReadLayout("<BARCODE2:6><SAMPLE2:x>");
				
		try {
			FastqWriterLayout fwdl = new FastqWriterLayout("<SAMPLE1>", "<READBAR1><READBAR2>", new ReadLayout[]{rl1, rl2}, false, FastqWriterLayout.DEFAULT_READNAME_DELIMITOR, FastqQualityFormat.Standard);
			FastqWriterLayout revl = new FastqWriterLayout("<SAMPLE2>", "<READBAR1><READBAR2>", new ReadLayout[]{rl1, rl2}, false, FastqWriterLayout.DEFAULT_READNAME_DELIMITOR, FastqQualityFormat.Standard);

			FastqRecord nfwd = fwdl.assembleRecord(new FastqRecord[]{fwd, rev},  null);
			FastqRecord nrev = revl.assembleRecord(new FastqRecord[]{fwd, rev},  null);

			Assert.assertEquals(readname+":CTGAGT:CTGAGT", nfwd.getReadName());
			Assert.assertEquals(readname+":CTGAGT:CTGAGT", nrev.getReadName());


			Assert.assertEquals(r1.substring(6), nfwd.getReadString());
			Assert.assertEquals(r2.substring(6), nrev.getReadString());

			Assert.assertEquals(q1.substring(6), nfwd.getBaseQualityString());
			Assert.assertEquals(q2.substring(6), nrev.getBaseQualityString());



		} catch (Jexception e) {
			log.error(ExceptionUtil.getStackTrace(e));
			Assert.fail("Got error with message : "+e.getMessage());
		}
	}

	@Test
	public final void testStandardPEWithUMI() {
		
		String readname = "NB501764:229:HV253BGX2:1:11101:6294:1041";
		String qualityheader = "+";
		String r1 = "CTGAGTACGTCAACACGACCAGTCAAACTTGACTGAAGGGATCTTCTTACATTTCCCTTCTTCAACATATTCTCTAATAT";
		String q1 = "AAAAA#EE//EEEEEEE/AEEAAEEEEEEAEEEA/EEAEEEEAEEE/EEEAEEEEEEEEEEAEEAEEEEAEEAEEAE/EA";
		String r2 = "CTGAGTGTTCGANCAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
		String q2 = "AAAAAEEA/AEE#EEE################################################################";
		FastqRecord fwd = new FastqRecord(readname, r1, qualityheader, q1);
		FastqRecord rev = new FastqRecord(readname, r2, qualityheader, q2);
		
		
		//consider 1 BC 1 UMI
		ReadLayout rl1 = new ReadLayout("<BARCODE1:6><SAMPLE1:x>");
		ReadLayout rl2 = new ReadLayout("<UMI1:10><SAMPLE2:x>");

		try {
			FastqWriterLayout fwdl = new FastqWriterLayout("<SAMPLE1>", "<READBAR1><UMI1>", new ReadLayout[]{rl1, rl2}, false, FastqWriterLayout.DEFAULT_READNAME_DELIMITOR, FastqQualityFormat.Standard);
			FastqWriterLayout revl = new FastqWriterLayout("<SAMPLE2>", "<UMI1><READBAR1>", new ReadLayout[]{rl1, rl2}, false, FastqWriterLayout.DEFAULT_READNAME_DELIMITOR, FastqQualityFormat.Standard);

			FastqRecord nfwd = fwdl.assembleRecord(new FastqRecord[]{fwd, rev},  null);
			FastqRecord nrev = revl.assembleRecord(new FastqRecord[]{fwd, rev},  null);

			Assert.assertEquals(readname+":CTGAGT:CTGAGTGTTC", nfwd.getReadName());
			Assert.assertEquals(readname+":CTGAGTGTTC:CTGAGT", nrev.getReadName());


			Assert.assertEquals(r1.substring(6), nfwd.getReadString());
			Assert.assertEquals(r2.substring(10), nrev.getReadString());

			Assert.assertEquals(q1.substring(6), nfwd.getBaseQualityString());
			Assert.assertEquals(q2.substring(10), nrev.getBaseQualityString());



		} catch (Jexception e) {
			log.error(ExceptionUtil.getStackTrace(e));
			Assert.fail("Got error with message : "+e.getMessage());
		}


	}
	
	@Test
	public final void testCrazyLayout() {
		
		String readname = "NB501764:229:HV253BGX2:1:11101:6294:1041";
		String qualityheader = "+";
		String r1 = "CTGAGTACGTCAACACGACCAGTCAAACTTGACTGAAGGGATCTTCTTACATTTCCCTTCTTCAACATATTCTCTAATAT";
		String q1 = "AAAAA#EE//EEEEEEE/AEEAAEEEEEEAEEEA/EEAEEEEAEEE/EEEAEEEEEEEEEEAEEAEEEEAEEAEEAE/EA";
		String r2 = "CTGAGTGTTCGANCAAACTGAAGGGNNGACCAGTCAAACTTNNNNACTGAAGGGNTTCNNACTGAAGGGACTGAAGGGAA";
		String q2 = "AAAAAEEA/AEE#EEE################################################################";
		final FastqRecord fwd = new FastqRecord(readname, r1, qualityheader, q1);
		final FastqRecord rev = new FastqRecord(readname, r2, qualityheader, q2);
		
		
		//consider 1 BC 1 UMI
		ReadLayout rl1 = new ReadLayout("<BARCODE1:6>NN<SAMPLE1:4>NN<UMI2:4>NNNN<SAMPLE2:x>");
		ReadLayout rl2 = new ReadLayout("<BARCODE2:6><SAMPLE3:24>NN<UMI1:8><SAMPLE4:38><BARCODE3:2>");

		try {
//			log.debug("testCrazyLayout => layout1");
//			FastqWriterLayout l1 = new FastqWriterLayout("<SAMPLE1><SAMPLE4>", "<READBAR1><READBAR2><READBAR3><UMI1><UMI2>", new ReadLayout[]{rl1, rl2}, false, FastqWriterLayout.DEFAULT_READNAME_DELIMITOR);
//			FastqRecord rec1 = l1.assembleRecord(new FastqRecord[]{fwd, rev}, null);
//			Assert.assertEquals(readname+":CTGAGT:CTGAGT:AA:GTCAAACT:ACGA", rec1.getReadName());
//			Assert.assertEquals("GTCA"+"TNNNNACTGAAGGGNTTCNNACTGAAGGGACTGAAGGG", rec1.getReadString());
			
			log.debug("testCrazyLayout => layout2");
			FastqWriterLayout l2 = new FastqWriterLayout("<BARCODE3><SAMPLE2><SAMPLE3>", "<UMI1><UMI2>", new ReadLayout[]{rl1, rl2}, false, FastqWriterLayout.DEFAULT_READNAME_DELIMITOR, FastqQualityFormat.Standard);
			FastqRecord rec2 = l2.assembleRecord(new FastqRecord[]{fwd, rev},  null);
			Assert.assertEquals(readname+":GTCAAACT:ACGA", rec2.getReadName());
			Assert.assertEquals("AA"+"TCAAACTTGACTGAAGGGATCTTCTTACATTTCCCTTCTTCAACATATTCTCTAATAT"+"GTTCGANCAAACTGAAGGGNNGAC", rec2.getReadString());
			
			log.debug("testCrazyLayout => layout3");
			FastqWriterLayout l3 = new FastqWriterLayout("<BARCODE1><BARCODE2><BARCODE3>", "",new ReadLayout[]{rl1, rl2}, false, FastqWriterLayout.DEFAULT_READNAME_DELIMITOR, FastqQualityFormat.Standard);
			FastqRecord rec3 = l3.assembleRecord(new FastqRecord[]{fwd, rev}, null);
			Assert.assertEquals(readname, rec3.getReadName());
			Assert.assertEquals("CTGAGT"+"CTGAGT"+"AA", rec3.getReadString());


		} catch (Jexception e) {
			log.error(ExceptionUtil.getStackTrace(e));
			Assert.fail("Got error with message : "+e.getMessage());
		}


	}
	
	public String printErrors(List<String> errors) {
		return StringUtil.mergeList(errors, "\n");
	}

	
	
	
}
