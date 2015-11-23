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

import java.io.File;
import java.net.URISyntaxException;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.embl.cg.utilitytools.utils.StringUtil;
import org.junit.Assert;
import org.junit.Test;

public class CommandLineParsingTest {


	private static Logger log = LoggerFactory.getLogger(CommandLineParsingTest	.class);

	
	
	/**
	 *  
	 */
	@Test
	public final void testRCHAR() {
		try {
			File bcFile = new File(CommandLineParsingTest.class.getResource("barcodes_dualindexing_C49NHACXX.txt").toURI()) ;
			File f1 = new File(CommandLineParsingTest.class.getResource("C49NHACXX_testCharles-lane8NoIndex_1_sequence_head_100.txt").toURI());
			File f2 = new File(CommandLineParsingTest.class.getResource("C49NHACXX_testCharles-lane8NoIndex_2_sequence_head_100.txt").toURI());
			File i1 = new File(CommandLineParsingTest.class.getResource("C49NHACXX_testCharles-lane8NoIndex_IDX1_sequence_head_100.txt").toURI());
			File i2 = new File(CommandLineParsingTest.class.getResource("C49NHACXX_testCharles-lane8NoIndex_IDX2_sequence_head_100.txt").toURI());

			
			//test RCHAR set to null
			String[] argv = new String[] {
					"F1="+f1.getAbsolutePath(), 
					"F2="+f2.getAbsolutePath(), 
					"I1="+i1.getAbsolutePath(), 
					"I2="+i2.getAbsolutePath(), 
					"BF="+bcFile.getAbsolutePath(),
					"BPOS=BOTH", //if not given must be turned to NONE 
					"LEN=8",
					"RCHAR=NULL",
					"SAME_HEADERS=true"
			};

			JemultiplexerIllumina j = new JemultiplexerIllumina();
			j.TEST_MODE_STOP_AFTER_PARSING = true;
			//parse
			Assert.assertEquals(0, j.instanceMain(argv));
			//
			Assert.assertTrue("RCHAR is not blank: '"+j.READ_NAME_REPLACE_CHAR+"'", j.READ_NAME_REPLACE_CHAR.equals(" "));
			
			//test not setting it
			argv = new String[] {
					"F1="+f1.getAbsolutePath(), 
					"F2="+f2.getAbsolutePath(), 
					"I1="+i1.getAbsolutePath(), 
					"I2="+i2.getAbsolutePath(), 
					"BF="+bcFile.getAbsolutePath(),
					"BPOS=BOTH", //if not given must be turned to NONE 
					"LEN=8"
			};

			j = new JemultiplexerIllumina();
			j.TEST_MODE_STOP_AFTER_PARSING = true;
			//parse
			Assert.assertEquals(0, j.instanceMain(argv));
			//sine same_header is true by default, we expect rchar to be set to :
			Assert.assertEquals(":", j.READ_NAME_REPLACE_CHAR ); 
			
			
			

		} catch (URISyntaxException e) {
			throw new RuntimeException(e);
		}
	}
	
	
	/**
	 * barcode at both ends, non-redundant 
	 */
	@Test
	public final void testPEIllumina() {

		try {
			File bcFile = new File(CommandLineParsingTest.class.getResource("barcodes_dualindexing_C49NHACXX.txt").toURI()) ;
			File f1 = new File(CommandLineParsingTest.class.getResource("C49NHACXX_testCharles-lane8NoIndex_1_sequence_head_100.txt").toURI());
			File f2 = new File(CommandLineParsingTest.class.getResource("C49NHACXX_testCharles-lane8NoIndex_2_sequence_head_100.txt").toURI());

			File i1 = new File(CommandLineParsingTest.class.getResource("C49NHACXX_testCharles-lane8NoIndex_IDX1_sequence_head_100.txt").toURI());
			File i2 = new File(CommandLineParsingTest.class.getResource("C49NHACXX_testCharles-lane8NoIndex_IDX2_sequence_head_100.txt").toURI());


			String[] argv = new String[] {
					"F1="+f1.getAbsolutePath(), 
					"F2="+f2.getAbsolutePath(), 
					"I1="+i1.getAbsolutePath(), 
					"I2="+i2.getAbsolutePath(), 
					"BF="+bcFile.getAbsolutePath(),
					//"XT=1:1", //XT must be turned to 0 
					"ZT=5:6",
					//"BPOS=BOTH", //if not given must be turned to NONE 
					"BRED=false", 
					"S=true", //irrelevant
					//"BCLEN=6", // must be forced to 0 
					"MM=3",
					"MMD=2",
					"Q=20"
			};

			JemultiplexerIllumina j = new JemultiplexerIllumina();
			j.TEST_MODE_STOP_AFTER_PARSING = true;
			//parse
			Assert.assertEquals(0, j.instanceMain(argv));

			//check j values are correctly positioned
			Assert.assertTrue(j.RUNNING_PAIRED_END);
			Assert.assertEquals( 0, j.XTRIMLEN_1.intValue());
			Assert.assertEquals( 0, j.XTRIMLEN_2.intValue());
			Assert.assertEquals( 5, j.ZTRIMLEN_1.intValue());
			Assert.assertEquals( 6, j.ZTRIMLEN_2.intValue());
			Assert.assertTrue(j.BARCODE_READ_POS == BarcodePosition.NONE);
			Assert.assertTrue(j.getBarcodeForSampleMatchingOptionValue() == BarcodePosition.NONE);
			Assert.assertFalse(j.REDUNDANT_BARCODES);
			Assert.assertEquals(0, j.BCLEN_1.intValue());
			Assert.assertEquals(0, j.BCLEN_2.intValue());
			Assert.assertEquals( 3, j.MAX_MISMATCHES_1.intValue());
			Assert.assertEquals( 3, j.MAX_MISMATCHES_2.intValue());
			Assert.assertEquals(2, j.MIN_MISMATCH_DELTA_1.intValue());
			Assert.assertEquals(2, j.MIN_MISMATCH_DELTA_2.intValue());
			Assert.assertEquals( 20, j.MIN_BASE_QUALITY_1.intValue());
			Assert.assertEquals( 20, j.MIN_BASE_QUALITY_2.intValue());
			Assert.assertEquals( Jemultiplexer.DEFAULT_CLIP_BARCODE, j.CLIP_BARCODE);
			Assert.assertEquals( Jemultiplexer.DEFAULT_ADD_BARCODE_TO_HEADER , j.ADD_BARCODE_TO_HEADER);
			Assert.assertEquals(Jemultiplexer.DEFAULT_GZIP_OUTPUTS, j.GZIP_OUTPUTS);
			Assert.assertEquals(Jemultiplexer.DEFAULT_KEEP_UNASSIGNED_READ,  j.KEEP_UNASSIGNED_READ);
			Assert.assertEquals( Jemultiplexer.DEFAULT_QUALITY_FORMAT, j.QUALITY_FORMAT);

			
			//single index file
			bcFile = new File(CommandLineParsingTest.class.getResource("barcodes_Illumina_simple_indexing_C49NHACXX.txt").toURI()) ;
			f1 = new File(CommandLineParsingTest.class.getResource("C49NHACXX_testCharles-lane8NoIndex_1_sequence_head_100.txt").toURI());
			f2 = new File(CommandLineParsingTest.class.getResource("C49NHACXX_testCharles-lane8NoIndex_2_sequence_head_100.txt").toURI());

			i1 = new File(CommandLineParsingTest.class.getResource("C49NHACXX_testCharles-lane8NoIndex_IDX1_sequence_head_100.txt").toURI());
			


			argv = new String[] {
					"F1="+f1.getAbsolutePath(), 
					"F2="+f2.getAbsolutePath(), 
					"I1="+i1.getAbsolutePath(), 
					"BF="+bcFile.getAbsolutePath(),
					//"XT=1:1", //XT must be turned to 0 
					"ZT=5:6",
					//"BPOS=BOTH", //if not given must be turned to NONE 
					"BRED=false", 
					"S=true", //irrelevant
					//"BCLEN=6", // must be forced to 0 
					"MM=3",
					"MMD=2",
					"Q=20"
			};

			j = new JemultiplexerIllumina();
			j.TEST_MODE_STOP_AFTER_PARSING = true;
			//parse
			Assert.assertEquals(0, j.instanceMain(argv));

			//check j values are correctly positioned
			Assert.assertTrue(j.RUNNING_PAIRED_END);
			Assert.assertEquals( 0, j.XTRIMLEN_1.intValue());
			Assert.assertEquals( 0, j.XTRIMLEN_2.intValue());
			Assert.assertEquals( 5, j.ZTRIMLEN_1.intValue());
			Assert.assertEquals( 6, j.ZTRIMLEN_2.intValue());
			Assert.assertTrue(j.BARCODE_READ_POS == BarcodePosition.NONE);
			Assert.assertTrue(j.getBarcodeForSampleMatchingOptionValue() == BarcodePosition.NONE);
			Assert.assertFalse(j.REDUNDANT_BARCODES);
			Assert.assertEquals(0, j.BCLEN_1.intValue());
			Assert.assertEquals(0, j.BCLEN_2.intValue());
			Assert.assertEquals( 3, j.MAX_MISMATCHES_1.intValue());
			Assert.assertEquals( 3, j.MAX_MISMATCHES_2.intValue());
			Assert.assertEquals(2, j.MIN_MISMATCH_DELTA_1.intValue());
			Assert.assertEquals(2, j.MIN_MISMATCH_DELTA_2.intValue());
			Assert.assertEquals( 20, j.MIN_BASE_QUALITY_1.intValue());
			Assert.assertEquals( 20, j.MIN_BASE_QUALITY_2.intValue());
			Assert.assertEquals( Jemultiplexer.DEFAULT_CLIP_BARCODE, j.CLIP_BARCODE);
			Assert.assertEquals( Jemultiplexer.DEFAULT_ADD_BARCODE_TO_HEADER , j.ADD_BARCODE_TO_HEADER);
			Assert.assertEquals(Jemultiplexer.DEFAULT_GZIP_OUTPUTS, j.GZIP_OUTPUTS);
			Assert.assertEquals(Jemultiplexer.DEFAULT_KEEP_UNASSIGNED_READ,  j.KEEP_UNASSIGNED_READ);
			Assert.assertEquals( Jemultiplexer.DEFAULT_QUALITY_FORMAT, j.QUALITY_FORMAT);
			
			//try define additional barcodes
			
		} catch (URISyntaxException e) {
			throw new RuntimeException(e);
		}
	}
	
	
	
	/**
	 * barcode at both ends, non-redundant 
	 */
	@Test
	public final void testPEComplexCase1() {

		try {
			File bcFile = new File(CommandLineParsingTest.class.getResource("barcodefiles/correct_barcodes_PE_both-ends_with_fnames.txt").toURI()) ;
			File f1 = new File(CommandLineParsingTest.class.getResource("C1WLBACXX_lane7_1_sequence.txt").toURI());
			File f2 = new File(CommandLineParsingTest.class.getResource("C1WLBACXX_lane7_2_sequence.txt").toURI());


			String[] argv = new String[] {
					"F1="+f1.getAbsolutePath(), 
					"F2="+f2.getAbsolutePath(), 
					"BF="+bcFile.getAbsolutePath(),
					"XT=1:1", 
					"ZT=5:5", 
					"BPOS=BOTH", 
					"BM=BOTH", 
					"BRED=false", 
					"S=true", //ignored
					//"BCLEN=6", // let discover
					"MM=3",
					"MMD=2",
					"Q=20"
			};

			Jemultiplexer j = new Jemultiplexer();
			j.TEST_MODE_STOP_AFTER_PARSING = true;
			//parse
			Assert.assertEquals(0, j.instanceMain(argv));

			//check j values are correctly positioned
			Assert.assertTrue(j.RUNNING_PAIRED_END);
			Assert.assertEquals( 1, j.XTRIMLEN_1.intValue());
			Assert.assertEquals( 1, j.XTRIMLEN_2.intValue());
			Assert.assertEquals( 5, j.ZTRIMLEN_1.intValue());
			Assert.assertEquals( 5, j.ZTRIMLEN_2.intValue());
			Assert.assertTrue(j.BARCODE_READ_POS == BarcodePosition.BOTH);
			Assert.assertTrue(j.BARCODE_FOR_SAMPLE_MATCHING == BarcodePosition.BOTH);
			Assert.assertFalse(j.REDUNDANT_BARCODES);
			Assert.assertEquals(6, j.BCLEN_1.intValue());
			Assert.assertEquals(7, j.BCLEN_2.intValue());
			Assert.assertEquals( 3, j.MAX_MISMATCHES_1.intValue());
			Assert.assertEquals( 3, j.MAX_MISMATCHES_2.intValue());
			Assert.assertEquals(2, j.MIN_MISMATCH_DELTA_1.intValue());
			Assert.assertEquals(2, j.MIN_MISMATCH_DELTA_2.intValue());
			Assert.assertEquals( 20, j.MIN_BASE_QUALITY_1.intValue());
			Assert.assertEquals( 20, j.MIN_BASE_QUALITY_2.intValue());
			Assert.assertEquals( Jemultiplexer.DEFAULT_CLIP_BARCODE, j.CLIP_BARCODE);
			Assert.assertEquals( Jemultiplexer.DEFAULT_ADD_BARCODE_TO_HEADER , j.ADD_BARCODE_TO_HEADER);
			Assert.assertEquals(Jemultiplexer.DEFAULT_GZIP_OUTPUTS, j.GZIP_OUTPUTS);
			Assert.assertEquals(Jemultiplexer.DEFAULT_KEEP_UNASSIGNED_READ,  j.KEEP_UNASSIGNED_READ);
			Assert.assertEquals( Jemultiplexer.DEFAULT_QUALITY_FORMAT, j.QUALITY_FORMAT);

			//should be just the same...
			argv = new String[] {
					"F1="+f1.getAbsolutePath(), 
					"F2="+f2.getAbsolutePath(), 
					"BF="+bcFile.getAbsolutePath(),
					"XT=1", 
					"ZT=5", 
					"BPOS=BOTH", 
					"BM=BOTH", 
					"BRED=false", 
					"BCLEN=6:7",
					"MM=3:4",
					"MMD=2:3",
					"Q=20:30"
			};

			j = new Jemultiplexer();
			j.TEST_MODE_STOP_AFTER_PARSING = true;
			//parse
			Assert.assertEquals(0, j.instanceMain(argv));

			//check j values are correctly positioned
			Assert.assertTrue(j.RUNNING_PAIRED_END);
			Assert.assertEquals( 1, j.XTRIMLEN_1.intValue());
			Assert.assertEquals( 1, j.XTRIMLEN_2.intValue());
			Assert.assertEquals( 5, j.ZTRIMLEN_1.intValue());
			Assert.assertEquals( 5, j.ZTRIMLEN_2.intValue());
			Assert.assertTrue(j.BARCODE_READ_POS == BarcodePosition.BOTH);
			Assert.assertTrue(j.BARCODE_FOR_SAMPLE_MATCHING == BarcodePosition.BOTH);
			Assert.assertFalse(j.REDUNDANT_BARCODES);
			Assert.assertEquals(6, j.BCLEN_1.intValue());
			Assert.assertEquals(7, j.BCLEN_2.intValue());
			Assert.assertEquals( 3, j.MAX_MISMATCHES_1.intValue());
			Assert.assertEquals( 4, j.MAX_MISMATCHES_2.intValue());
			Assert.assertEquals(2, j.MIN_MISMATCH_DELTA_1.intValue());
			Assert.assertEquals(3, j.MIN_MISMATCH_DELTA_2.intValue());
			Assert.assertEquals( 20, j.MIN_BASE_QUALITY_1.intValue());
			Assert.assertEquals( 30, j.MIN_BASE_QUALITY_2.intValue());
			Assert.assertEquals( Jemultiplexer.DEFAULT_CLIP_BARCODE, j.CLIP_BARCODE);
			Assert.assertEquals( Jemultiplexer.DEFAULT_ADD_BARCODE_TO_HEADER , j.ADD_BARCODE_TO_HEADER);
			Assert.assertEquals(Jemultiplexer.DEFAULT_GZIP_OUTPUTS, j.GZIP_OUTPUTS);
			Assert.assertEquals(Jemultiplexer.DEFAULT_KEEP_UNASSIGNED_READ,  j.KEEP_UNASSIGNED_READ);
			Assert.assertEquals( Jemultiplexer.DEFAULT_QUALITY_FORMAT, j.QUALITY_FORMAT);

			
			//error in BCLEN => crash
			argv = new String[] {
					"F1="+f1.getAbsolutePath(), 
					"F2="+f2.getAbsolutePath(), 
					"BF="+bcFile.getAbsolutePath(),
					"BPOS=BOTH", 
					"BM=BOTH", 
					"BRED=false", 
					"BCLEN=6:6", //wrong second length
			};

			j = new Jemultiplexer();
			j.TEST_MODE_STOP_AFTER_PARSING = true;
			//parse
			Assert.assertFalse(0 == j.instanceMain(argv));

			
		} catch (URISyntaxException e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * barcode at both ends, non-redundant, only one used for look up, the other one is rdm 
	 */
	@Test
	public final void testPEComplexCase2() {

		try {
			File bcFile = new File(CommandLineParsingTest.class.getResource("barcodefiles/correct_barcodes_PE_multi_with_fnames.txt").toURI()) ;
			File f1 = new File(CommandLineParsingTest.class.getResource("C1WLBACXX_lane7_1_sequence.txt").toURI());
			File f2 = new File(CommandLineParsingTest.class.getResource("C1WLBACXX_lane7_2_sequence.txt").toURI());

			for(BarcodePosition bp : new BarcodePosition[]{BarcodePosition.READ_1, BarcodePosition.READ_2}){
				String bclenOPT = bp == BarcodePosition.READ_1 ? "BCLEN=6:8" : "BCLEN=8:6";
				
				
				String[] argv = new String[] {
						"F1="+f1.getAbsolutePath(), 
						"F2="+f2.getAbsolutePath(), 
						"BF="+bcFile.getAbsolutePath(),
						"XT=2:1", 
						"ZT=5:3", 
						"BPOS=BOTH", 
						"BM="+bp, 
						"BRED=false", 
						bclenOPT, //BLEN must be give if we want to CLIP barcode ie how do we know BLEN of the random side otherwise?? 
						"MM=3",
						"MMD=2",
						"Q=20"
				};

				Jemultiplexer j = new Jemultiplexer();
				j.TEST_MODE_STOP_AFTER_PARSING = true;
				//parse
				log.warn("exec with args : "+StringUtil.mergeArray(argv, " ; "));
				int retcode = j.instanceMain(argv);
				log.warn("AAA RET code  : "+retcode);
				Assert.assertEquals(0, retcode);

				//check j values are correctly positioned
				Assert.assertTrue(j.RUNNING_PAIRED_END);
				Assert.assertEquals( 2, j.XTRIMLEN_1.intValue());
				Assert.assertEquals( 1, j.XTRIMLEN_2.intValue());
				Assert.assertEquals( 5, j.ZTRIMLEN_1.intValue());
				Assert.assertEquals( 3, j.ZTRIMLEN_2.intValue());
				Assert.assertTrue(j.BARCODE_READ_POS == BarcodePosition.BOTH);
				Assert.assertTrue(j.BARCODE_FOR_SAMPLE_MATCHING == bp);
				Assert.assertFalse(j.REDUNDANT_BARCODES);
				if(bp == BarcodePosition.READ_1){
					Assert.assertEquals(6, j.BCLEN_1.intValue());
					Assert.assertEquals(8, j.BCLEN_2.intValue());
				}else{
					Assert.assertEquals(8, j.BCLEN_1.intValue());
					Assert.assertEquals(6, j.BCLEN_2.intValue());
				}
				Assert.assertEquals( 3, j.MAX_MISMATCHES_1.intValue());
				Assert.assertEquals( 3, j.MAX_MISMATCHES_2.intValue());
				Assert.assertEquals(2, j.MIN_MISMATCH_DELTA_1.intValue());
				Assert.assertEquals(2, j.MIN_MISMATCH_DELTA_2.intValue());
				Assert.assertEquals( 20, j.MIN_BASE_QUALITY_1.intValue());
				Assert.assertEquals( 20, j.MIN_BASE_QUALITY_2.intValue());
				
				//
				
				argv = new String[] {
						"F1="+f1.getAbsolutePath(), 
						"F2="+f2.getAbsolutePath(), 
						"BF="+bcFile.getAbsolutePath(),
						"XT=2", //since BPOS=BOTH ,will apply to both ends
						"ZT=5",  //since BPOS=BOTH ,will apply to both ends
						"BPOS=BOTH", 
						"BM="+bp, 
						"BRED=false", 
						bclenOPT, //BLEN must be give if we want to CLIP barcode ie how do we know BLEN of the random side otherwise?? 
						"MM=3",
						"MMD=2",
						"Q=20"
				};

				
				j = new Jemultiplexer();
				j.TEST_MODE_STOP_AFTER_PARSING = true;
				//parse
				log.warn("exec with args : "+StringUtil.mergeArray(argv, " ; "));
				retcode = j.instanceMain(argv);
				log.warn("BBB RET code  : "+retcode);
				Assert.assertEquals(0, retcode);

				//check j values are correctly positioned
				Assert.assertTrue(j.RUNNING_PAIRED_END);
				Assert.assertEquals( 2, j.XTRIMLEN_1.intValue());
				Assert.assertEquals( 2, j.XTRIMLEN_2.intValue());
				Assert.assertEquals( 5, j.ZTRIMLEN_1.intValue());
				Assert.assertEquals( 5, j.ZTRIMLEN_2.intValue());
				Assert.assertTrue(j.BARCODE_READ_POS == BarcodePosition.BOTH);
				Assert.assertTrue(j.BARCODE_FOR_SAMPLE_MATCHING == bp);
				Assert.assertFalse(j.REDUNDANT_BARCODES);
				if(bp == BarcodePosition.READ_1){
					Assert.assertEquals(6, j.BCLEN_1.intValue());
					Assert.assertEquals(8, j.BCLEN_2.intValue());
				}else{
					Assert.assertEquals(8, j.BCLEN_1.intValue());
					Assert.assertEquals(6, j.BCLEN_2.intValue());
				}Assert.assertEquals( 3, j.MAX_MISMATCHES_1.intValue());
				Assert.assertEquals( 3, j.MAX_MISMATCHES_2.intValue());
				Assert.assertEquals(2, j.MIN_MISMATCH_DELTA_1.intValue());
				Assert.assertEquals(2, j.MIN_MISMATCH_DELTA_2.intValue());
				Assert.assertEquals( 20, j.MIN_BASE_QUALITY_1.intValue());
				Assert.assertEquals( 20, j.MIN_BASE_QUALITY_2.intValue());
				

				
				// test it fails when on bclen is not passed with CLIP=true
				
				argv = new String[] {
						"F1="+f1.getAbsolutePath(), 
						"F2="+f2.getAbsolutePath(), 
						"BF="+bcFile.getAbsolutePath(),
						"XT=2", //since BPOS=BOTH ,will apply to both ends
						"ZT=5",  //since BPOS=BOTH ,will apply to both ends
						"BPOS=BOTH", 
						"BM="+bp, 
						"BRED=false", 
						//"BCLEN=6:8", //BLEN must be give if we want to CLIP barcode ie how do we know BLEN of the random side otherwise?? 
						"C=true"
				};

				
				j = new Jemultiplexer();
				j.TEST_MODE_STOP_AFTER_PARSING = true;
				//parse
				log.warn("exec with args : "+StringUtil.mergeArray(argv, " ; "));
				retcode = j.instanceMain(argv);
				log.warn("CCC  RET code  : "+retcode);
				Assert.assertTrue(0 != retcode);

				
				
			}
			
		} catch (URISyntaxException e) {
			throw new RuntimeException(e);
		}
	}

	
	/**
	 * barcode at both ends, non-redundant, only one used for look up, the other one is rdm 
	 * check how BCLEN behaves
	 */
	@Test
	public final void testPEComplexBCLenTest() {

		try {
			File bcFile = new File(CommandLineParsingTest.class.getResource("barcodefiles/correct_barcodes_PE_multi_with_fnames.txt").toURI()) ;
			File f1 = new File(CommandLineParsingTest.class.getResource("C1WLBACXX_lane7_1_sequence.txt").toURI());
			File f2 = new File(CommandLineParsingTest.class.getResource("C1WLBACXX_lane7_2_sequence.txt").toURI());


			String[] argv = new String[] {
					"F1="+f1.getAbsolutePath(), 
					"F2="+f2.getAbsolutePath(), 
					"BF="+bcFile.getAbsolutePath(),
					"BPOS=BOTH", 
					"BM=READ_1", 
					"BRED=false", 
					"BCLEN=6" //only give one length => it must match both the length in BC file and will be the same used for the 2d BC  
			};

			Jemultiplexer j = new Jemultiplexer();
			j.TEST_MODE_STOP_AFTER_PARSING = true;
			//parse
			Assert.assertEquals(0, j.instanceMain(argv));

			//check j values are correctly positioned
			Assert.assertTrue(j.RUNNING_PAIRED_END);
			Assert.assertTrue(j.BARCODE_READ_POS == BarcodePosition.BOTH);
			Assert.assertTrue(j.BARCODE_FOR_SAMPLE_MATCHING == BarcodePosition.READ_1);
			Assert.assertFalse(j.REDUNDANT_BARCODES);
			Assert.assertEquals(6, j.BCLEN_1.intValue());
			Assert.assertEquals(6, j.BCLEN_2.intValue());

			
			//should fail with another length
			argv = new String[] {
					"F1="+f1.getAbsolutePath(), 
					"F2="+f2.getAbsolutePath(), 
					"BF="+bcFile.getAbsolutePath(),
					"BPOS=BOTH", 
					"BM=READ_1", 
					"BRED=false", 
					"BCLEN=8" //BC in file are 6 bp  
			};

			j = new Jemultiplexer();
			j.TEST_MODE_STOP_AFTER_PARSING = true;
			//parse
			Assert.assertFalse(0 == j.instanceMain(argv));

			
			
			//should fail with no length
			argv = new String[] {
					"F1="+f1.getAbsolutePath(), 
					"F2="+f2.getAbsolutePath(), 
					"BF="+bcFile.getAbsolutePath(),
					"BPOS=BOTH", 
					"BM=READ_1", 
					"BRED=false" 
					  
			};

			j = new Jemultiplexer();
			j.TEST_MODE_STOP_AFTER_PARSING = true;
			//parse
			Assert.assertFalse(0 == j.instanceMain(argv));
			
		} catch (URISyntaxException e) {
			throw new RuntimeException(e);
		}
	}

	
	

	/**
	 * barcode at both ends, redundant ; or at one end only
	 */
	@Test
	public final void testPESimpleCase() {

		try {
			File bcFile = new File(CommandLineParsingTest.class.getResource("C1WLBACXX_lane7.bs").toURI()) ;
			File f1 = new File(CommandLineParsingTest.class.getResource("C1WLBACXX_lane7_1_sequence.txt").toURI());
			File f2 = new File(CommandLineParsingTest.class.getResource("C1WLBACXX_lane7_2_sequence.txt").toURI());


			String[] argv = new String[] {
					"F1="+f1.getAbsolutePath(), 
					"F2="+f2.getAbsolutePath(), 
					"BF="+bcFile.getAbsolutePath(),
					"XT=2", 
					"ZT=1", 
					"BPOS=BOTH", 
					"BM=BOTH", 
					"BRED=true", 
					"S=true", 
					"BCLEN=6",
					"MM=3",
					"MMD=2",
					"Q=30",
					"ADD=false",
					"GZ=false",
					"UN=false",
			};

			Jemultiplexer j = new Jemultiplexer();
			j.TEST_MODE_STOP_AFTER_PARSING = true;
			//parse
			Assert.assertEquals(0, j.instanceMain(argv));

			//check j values are correctly positioned
			Assert.assertTrue(j.RUNNING_PAIRED_END);
			Assert.assertEquals( 2, j.XTRIMLEN_1.intValue());
			Assert.assertEquals( 2, j.XTRIMLEN_2.intValue());
			Assert.assertEquals( 1, j.ZTRIMLEN_1.intValue());
			Assert.assertEquals( 1, j.ZTRIMLEN_2.intValue());
			Assert.assertTrue(j.BARCODE_READ_POS == BarcodePosition.BOTH);
			Assert.assertTrue(j.BARCODE_FOR_SAMPLE_MATCHING == BarcodePosition.BOTH);
			Assert.assertTrue(j.REDUNDANT_BARCODES);
			Assert.assertEquals(6, j.BCLEN_1.intValue());
			Assert.assertEquals(6, j.BCLEN_2.intValue());
			Assert.assertEquals( 3, j.MAX_MISMATCHES_1.intValue());
			Assert.assertEquals( 3, j.MAX_MISMATCHES_2.intValue());
			Assert.assertEquals(2, j.MIN_MISMATCH_DELTA_1.intValue());
			Assert.assertEquals(2, j.MIN_MISMATCH_DELTA_2.intValue());
			Assert.assertEquals( 30, j.MIN_BASE_QUALITY_1.intValue());
			Assert.assertEquals( 30, j.MIN_BASE_QUALITY_2.intValue());
			Assert.assertTrue(j.STRICT);
			Assert.assertEquals( Jemultiplexer.DEFAULT_CLIP_BARCODE, j.CLIP_BARCODE);
			Assert.assertFalse( j.ADD_BARCODE_TO_HEADER);
			Assert.assertFalse(j.GZIP_OUTPUTS);
			Assert.assertFalse( j.KEEP_UNASSIGNED_READ);
			Assert.assertEquals( Jemultiplexer.DEFAULT_QUALITY_FORMAT, j.QUALITY_FORMAT);


			/*
			 * with different options that gives just a warning
			 */

			argv = new String[] {
					"F1="+f1.getAbsolutePath(), 
					"F2="+f2.getAbsolutePath(), 
					"BF="+bcFile.getAbsolutePath(),
					"XT=2:3", 
					"ZT=1:2", 
					"BPOS=BOTH", 
					"BM=BOTH", 
					"BRED=true", 
					"S=true", 
					"BCLEN=6",
					"MM=3",
					"MMD=2",
					"Q=30",
					"ADD=false",
					"GZ=false",
					"UN=false",
			};

			j = new Jemultiplexer();
			j.TEST_MODE_STOP_AFTER_PARSING = true;
			//parse
			Assert.assertEquals(0, j.instanceMain(argv));

			//check j values are correctly positioned
			Assert.assertTrue(j.RUNNING_PAIRED_END);
			Assert.assertEquals( 2, j.XTRIMLEN_1.intValue());
			Assert.assertEquals( 3, j.XTRIMLEN_2.intValue());
			Assert.assertEquals( 1, j.ZTRIMLEN_1.intValue());
			Assert.assertEquals( 2, j.ZTRIMLEN_2.intValue());
			Assert.assertTrue(j.BARCODE_READ_POS == BarcodePosition.BOTH);
			Assert.assertTrue(j.BARCODE_FOR_SAMPLE_MATCHING == BarcodePosition.BOTH);
			Assert.assertTrue(j.REDUNDANT_BARCODES);
			Assert.assertEquals(6, j.BCLEN_1.intValue());
			Assert.assertEquals(6, j.BCLEN_2.intValue());
			Assert.assertEquals( 3, j.MAX_MISMATCHES_1.intValue());
			Assert.assertEquals( 3, j.MAX_MISMATCHES_2.intValue());
			Assert.assertEquals(2, j.MIN_MISMATCH_DELTA_1.intValue());
			Assert.assertEquals(2, j.MIN_MISMATCH_DELTA_2.intValue());
			Assert.assertEquals( 30, j.MIN_BASE_QUALITY_1.intValue());
			Assert.assertEquals( 30, j.MIN_BASE_QUALITY_2.intValue());
			Assert.assertTrue(j.STRICT);
			Assert.assertEquals( Jemultiplexer.DEFAULT_CLIP_BARCODE, j.CLIP_BARCODE);
			Assert.assertFalse( j.ADD_BARCODE_TO_HEADER);
			Assert.assertFalse(j.GZIP_OUTPUTS);
			Assert.assertFalse( j.KEEP_UNASSIGNED_READ);
			Assert.assertEquals( Jemultiplexer.DEFAULT_QUALITY_FORMAT, j.QUALITY_FORMAT);


			for(String s : new String[]{"MM=3:2", "MMD=1:1", "Q=30:45", "XT=2:3", "ZT=1:2" }){
				argv = new String[]{
						"F1="+f1.getAbsolutePath(), 
						"F2="+f2.getAbsolutePath(), 
						"BF="+bcFile.getAbsolutePath(),
						"BPOS=BOTH", 
						"BM=BOTH", 
						"BRED=true", 
						"S=true", 
						s

				};

				j = new Jemultiplexer();
				j.TEST_MODE_STOP_AFTER_PARSING = true;
				//parse
				Assert.assertTrue("Should have validated with "+s, 0 == j.instanceMain(argv));
			}
			
			
			//with options that crash
			for(String s : new String[]{"BCLEN=6:7" }){
				argv = new String[]{
						"F1="+f1.getAbsolutePath(), 
						"F2="+f2.getAbsolutePath(), 
						"BF="+bcFile.getAbsolutePath(),
						"BPOS=BOTH", 
						"BM=BOTH", 
						"BRED=true", 
						"S=true", 
						s

				};

				j = new Jemultiplexer();
				j.TEST_MODE_STOP_AFTER_PARSING = true;
				//parse
				Assert.assertFalse("Should not have validated with "+s, 0 == j.instanceMain(argv));
			}

			
			/*
			 * now have a barcode at only one side
			 */
			for(BarcodePosition bp : new BarcodePosition[]{BarcodePosition.READ_1, BarcodePosition.READ_2}){
				argv = new String[] {
						"F1="+f1.getAbsolutePath(), 
						"F2="+f2.getAbsolutePath(), 
						"BF="+bcFile.getAbsolutePath(),
						"XT=2", 
						"ZT=0:8", 
						"BPOS="+bp.toString(),  
						"BRED=true", //should be ignored, leave it to check that no failure is observed
						"S=true", //should be ignored, leave it to check that no failure is observed
						"BCLEN=6", 
						"MM=3",
						"MMD=2",
						"Q=30",
						"ADD=true",
						"GZ=true",
						"UN=true",
				};

				j = new Jemultiplexer();
				j.TEST_MODE_STOP_AFTER_PARSING = true;
				//parse
				Assert.assertEquals(0, j.instanceMain(argv));

				//check j values are correctly positioned
				Assert.assertTrue(j.RUNNING_PAIRED_END);
				Assert.assertEquals( bp==BarcodePosition.READ_1? 2:0, j.XTRIMLEN_1.intValue());
				Assert.assertEquals( bp==BarcodePosition.READ_1? 0:2, j.XTRIMLEN_2.intValue());
				Assert.assertEquals( 0, j.ZTRIMLEN_1.intValue()); //this won't change with bp as we gave a 0:8 syntax
				Assert.assertEquals(  8, j.ZTRIMLEN_2.intValue()); //this won't change with bp as we gave a 0:8 syntax
				Assert.assertTrue(j.BARCODE_READ_POS == bp);
				Assert.assertTrue(j.BARCODE_FOR_SAMPLE_MATCHING == bp);
				Assert.assertEquals(bp==BarcodePosition.READ_1? 6:0, j.BCLEN_1.intValue());
				Assert.assertEquals(bp==BarcodePosition.READ_1? 0:6, j.BCLEN_2.intValue());
				Assert.assertEquals( 3, j.MAX_MISMATCHES_1.intValue());
				Assert.assertEquals( 3, j.MAX_MISMATCHES_2.intValue());
				Assert.assertEquals(2, j.MIN_MISMATCH_DELTA_1.intValue());
				Assert.assertEquals(2, j.MIN_MISMATCH_DELTA_2.intValue());
				Assert.assertEquals( 30, j.MIN_BASE_QUALITY_1.intValue());
				Assert.assertEquals( 30, j.MIN_BASE_QUALITY_2.intValue());
				Assert.assertEquals( Jemultiplexer.DEFAULT_CLIP_BARCODE, j.CLIP_BARCODE);
				Assert.assertTrue( j.ADD_BARCODE_TO_HEADER);
				Assert.assertTrue(j.GZIP_OUTPUTS);
				Assert.assertTrue( j.KEEP_UNASSIGNED_READ);
				Assert.assertEquals( Jemultiplexer.DEFAULT_QUALITY_FORMAT, j.QUALITY_FORMAT);

				//	now options that are forbidden in such a situation
				for(String s : new String[]{"BCLEN=6:7", "MM=3:4", "MMD=3:1", "Q=30:43" }){
					argv = new String[] {
							"F1="+f1.getAbsolutePath(), 
							"F2="+f2.getAbsolutePath(), 
							"BF="+bcFile.getAbsolutePath(),
							"XT=2:4", 
							"ZT=0:8", 
							"BPOS="+bp.toString(),  
							"BRED=true", //should be ignored, leave it to check that no failure is observed
							"S=true", //should be ignored, leave it to check that no failure is observed
							s
					};

					j = new Jemultiplexer();
					j.TEST_MODE_STOP_AFTER_PARSING = true;
					//parse
					Assert.assertFalse("Should not have validated with "+s, 0 == j.instanceMain(argv));

				}
			}
		} catch (URISyntaxException e) {
			throw new RuntimeException(e);
		}
	}


	@Test
	public final void testMinimalSE() {
		try {
			File bcFile = new File(CommandLineParsingTest.class.getResource("C1WLBACXX_lane7.bs").toURI()) ;
			File f1 = new File(CommandLineParsingTest.class.getResource("C1WLBACXX_lane7_1_sequence.txt").toURI());
			String[] argv = {
					"F1="+f1.getAbsolutePath(), 
					"BF="+bcFile.getAbsolutePath() 
			};

			Jemultiplexer j = new Jemultiplexer();
			j.TEST_MODE_STOP_AFTER_PARSING = true;
			//parse
			Assert.assertEquals(0, j.instanceMain(argv));

			//check j values are correctly positioned
			Assert.assertTrue(!j.RUNNING_PAIRED_END);
			Assert.assertTrue("Wrong BARCODE_READ_POS value", j.BARCODE_READ_POS == BarcodePosition.READ_1);
			Assert.assertTrue(j.BARCODE_FOR_SAMPLE_MATCHING == BarcodePosition.READ_1);
			Assert.assertEquals("Wrong BCLEN_1 in SE mode", 6, j.BCLEN_1.intValue());
			Assert.assertEquals( Jemultiplexer.DEFAULT_MAX_MISMATCHES, j.MAX_MISMATCHES_1);
			Assert.assertEquals( Jemultiplexer.DEFAULT_MIN_MISMATCH_DELTA, j.MIN_MISMATCH_DELTA_1);
			Assert.assertEquals( Jemultiplexer.DEFAULT_MIN_BASE_QUALITY, j.MIN_BASE_QUALITY_1);
			Assert.assertEquals( Jemultiplexer.DEFAULT_XTRIMLEN, j.XTRIMLEN_1);
			Assert.assertEquals( Jemultiplexer.DEFAULT_ZTRIMLEN, j.ZTRIMLEN_1);
			Assert.assertEquals( Jemultiplexer.DEFAULT_CLIP_BARCODE, j.CLIP_BARCODE);
			Assert.assertEquals( Jemultiplexer.DEFAULT_ADD_BARCODE_TO_HEADER, j.ADD_BARCODE_TO_HEADER);
			Assert.assertEquals( Jemultiplexer.DEFAULT_QUALITY_FORMAT, j.QUALITY_FORMAT);
			Assert.assertEquals( Jemultiplexer.DEFAULT_GZIP_OUTPUTS, j.GZIP_OUTPUTS);
			Assert.assertEquals( Jemultiplexer.DEFAULT_KEEP_UNASSIGNED_READ, j.KEEP_UNASSIGNED_READ);


		} catch (URISyntaxException e) {
			throw new RuntimeException(e);
		}
	}

	@Test
	public final void testMinimalSEWithRCHAR() {
		try {
			File bcFile = new File(CommandLineParsingTest.class.getResource("C1WLBACXX_lane7.bs").toURI()) ;
			File f1 = new File(CommandLineParsingTest.class.getResource("C1WLBACXX_lane7_1_sequence.txt").toURI());
			for(String c : new String[]{":", "_", "-"}){
				String[] argv = {
						"F1="+f1.getAbsolutePath(), 
						"BF="+bcFile.getAbsolutePath(),
						"RCHAR="+c 
				};

				Jemultiplexer j = new Jemultiplexer();
				j.TEST_MODE_STOP_AFTER_PARSING = true;
				//parse
				Assert.assertEquals(0, j.instanceMain(argv));

				//check RCHAR
				Assert.assertEquals("Wrong READ_NAME_REPLACE_CHAR", c, j.READ_NAME_REPLACE_CHAR);
			}
		} catch (URISyntaxException e) {
			throw new RuntimeException(e);
		}
	}

	
	@Test
	public final void testMinimalSEWithDifferentOptions() {
		try {
			File bcFile = new File(CommandLineParsingTest.class.getResource("C1WLBACXX_lane7.bs").toURI()) ;
			File f1 = new File(CommandLineParsingTest.class.getResource("C1WLBACXX_lane7_1_sequence.txt").toURI());
			String[] argv = {
					"F1="+f1.getAbsolutePath(), 
					"BF="+bcFile.getAbsolutePath(),
					"XT=2", 
					"ZT=1", 
					"BPOS=READ_2", //should be ignored and reset correctly
					"BRED=false", //should be ignored 
					"BM=BOTH", //should be ignored and reset correctly
					"S=true", //should be ignored
					"BCLEN=6",
					"MM=3",
					"MMD=2",
					"Q=30",
					"ADD=false",
					"GZ=false",
					"UN=false"
			};

			Jemultiplexer j = new Jemultiplexer();
			j.TEST_MODE_STOP_AFTER_PARSING = true;
			//parse
			Assert.assertEquals(0, j.instanceMain(argv));

			//check j values are correctly positioned
			Assert.assertTrue(!j.RUNNING_PAIRED_END);
			Assert.assertEquals( 2, j.XTRIMLEN_1.intValue());
			Assert.assertEquals( 1, j.ZTRIMLEN_1.intValue());
			Assert.assertTrue(j.BARCODE_READ_POS == BarcodePosition.READ_1);
			Assert.assertTrue(j.BARCODE_FOR_SAMPLE_MATCHING == BarcodePosition.READ_1);
			Assert.assertEquals("Wrong BCLEN_1 in SE mode", 6, j.BCLEN_1.intValue());
			Assert.assertEquals( 3, j.MAX_MISMATCHES_1.intValue());
			Assert.assertEquals(2, j.MIN_MISMATCH_DELTA_1.intValue());
			Assert.assertEquals( 30, j.MIN_BASE_QUALITY_1.intValue());
			Assert.assertEquals( Jemultiplexer.DEFAULT_CLIP_BARCODE, j.CLIP_BARCODE);
			Assert.assertFalse( j.ADD_BARCODE_TO_HEADER);
			Assert.assertFalse(j.GZIP_OUTPUTS);
			Assert.assertFalse( j.KEEP_UNASSIGNED_READ);
			Assert.assertEquals( Jemultiplexer.DEFAULT_QUALITY_FORMAT, j.QUALITY_FORMAT);


		} catch (URISyntaxException e) {
			throw new RuntimeException(e);
		}
	}

	
	
	
	@Test
	public final void testMinimalSEWithWrongOptions() {
		try {
			File bcFile = new File(CommandLineParsingTest.class.getResource("C1WLBACXX_lane7.bs").toURI()) ;
			File f1 = new File(CommandLineParsingTest.class.getResource("C1WLBACXX_lane7_1_sequence.txt").toURI());
			String[] argv = {
					"F1="+f1.getAbsolutePath(), 
					"BF="+bcFile.getAbsolutePath(),
					"XT=2:2", 
					"ZT=1:1", 
					"BCLEN=6",
					"MM=3",
					"MMD=2",
					"Q=30",
			};

			Jemultiplexer j = new Jemultiplexer();
			j.TEST_MODE_STOP_AFTER_PARSING = true;
			//parse
			Assert.assertFalse(0 == j.instanceMain(argv));



		} catch (URISyntaxException e) {
			throw new RuntimeException(e);
		}
	}

}
