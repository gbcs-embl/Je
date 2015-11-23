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
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeSet;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.embl.cg.utilitytools.utils.ExceptionUtil;
import org.embl.cg.utilitytools.utils.StringUtil;
import org.embl.cg.utilitytools.utils.parser.csv.ParsingException;
import org.junit.Assert;
import org.junit.Test;

public class BarcodeFileTest {
	private static Logger log = LoggerFactory.getLogger(BarcodeFileTest.class);

	
	@Test
	public final void testShortIncorrectFormats() {
		Jemultiplexer jse = null;
		BarcodeValidator v = null;
		List<String> errors = null;

		try {

			//both these files shoudl work in SE and PE (simple cases)
			File f1 = new File(BarcodeFileTest.class.getResource("barcodefiles/incorrect_barcodes_barcodetooshort.txt").toURI()) ;
			File f2 = new File(BarcodeFileTest.class.getResource("barcodefiles/incorrect_barcodes_simple_redundantbarcode.txt").toURI()) ;
			File f3 = new File(BarcodeFileTest.class.getResource("barcodefiles/incorrect_barcodes_simple_barcodetooshort.txt").toURI()) ;
			File f4 = new File(BarcodeFileTest.class.getResource("barcodefiles/incorrect_barcodes_simple_colnumbers.txt").toURI()) ;
			File [] fileArr = new File[]{f1, f2, f3, f4}; 
			for (int i = 0; i < fileArr.length; i++) {
				File f = fileArr[i];
				log.info("using "+f.getName());
				/*
				 * single end
				 */
				for (int l: new int[]{0,6} ){
					jse = new Jemultiplexer();
					jse.RUNNING_PAIRED_END = false;
					jse.BARCODE_FILE = f;
					jse.BARCODE_READ_POS = BarcodePosition.READ_1;
					jse.BARCODE_FOR_SAMPLE_MATCHING = BarcodePosition.READ_1;
					jse.BCLEN_1 = l;

					v = new BarcodeValidator(f, jse);
					try{
						errors = v.validate();
						Assert.assertTrue("Should have got an error here!", errors.size() > 0); 
					}catch (RuntimeException e) {
						// this is also ok , ie some malformed files do launch exceptions
						Assert.assertTrue(e.getCause() instanceof ParsingException);
					}
				}
				/*
				 * this format should also work for simple PAIRED END but not for BARCODE_FOR_SAMPLE_MATCHING = BOTH and 
				 * REDUNDANT_BARCODES=false as a bc1:bc2 should be given in such situation 
				 */
				for(BarcodePosition bp : new BarcodePosition[]{BarcodePosition.READ_1, BarcodePosition.READ_2, BarcodePosition.BOTH}){
					log.debug(bp.toString());
					for(int l : new int []{0,6}){

						jse = new Jemultiplexer();
						jse.RUNNING_PAIRED_END = true;
						jse.BARCODE_FILE = f;
						jse.BARCODE_READ_POS =bp;
						jse.BARCODE_FOR_SAMPLE_MATCHING = bp;


						log.debug("preset BCLEN to "+l);
						if(bp != BarcodePosition.READ_2)
							jse.BCLEN_1 = l;
						if(bp != BarcodePosition.READ_1)
							jse.BCLEN_2 = l;

						v = new BarcodeValidator(f, jse);
						try{
							log.info("VALIDATE file "+f.getName()+" and BPOS="+bp);
							errors = v.validate();
							Assert.assertTrue("Should have got an error here with bc file "+f.getName()+" and BPOS="+bp, errors.size() > 0); 
						}catch (RuntimeException e) {
							// this is also ok , ie some malformed files do launch exceptions
							Assert.assertTrue(e.getCause() instanceof ParsingException);
						}
					}
				}

			}
		} catch (URISyntaxException e) {
			//this is not a test failure!
			throw new RuntimeException(e);
		}

	}
	
	
	
	@Test
	public final void testShortCorrectFormats() {
		Jemultiplexer jse = null;
		BarcodeValidator v = null;
		List<String> errors = null;
		
		
		try {
			
			//both these files shoudl work in SE and PE (simple cases)
			File f1 = new File(BarcodeFileTest.class.getResource("barcodefiles/correct_barcodes_simple.txt").toURI()) ;
			File f2 = new File(BarcodeFileTest.class.getResource("barcodefiles/correct_barcodes_multi.txt").toURI()) ;
			File f3 = new File(BarcodeFileTest.class.getResource("barcodefiles/correct_barcodes_simple_iCLIP.txt").toURI()) ;
			File [] fileArr = new File[]{f1, f2, f3}; //NEVER CHANGE ORDER AS VALIDATION MAPS ARE BUILD ASSUMING THIS ORDER
			for (int i = 0; i < fileArr.length; i++) {
				File f = fileArr[i];
				
				int realblen = i==2 ? 11 : 6;
				
				//init validation maps
				HashMap<String, Set<String>> validation_sample2barcode = new HashMap<String, Set<String>>();
				if(i == 0){
					validation_sample2barcode.put("sample1", new TreeSet<String>( Arrays.asList(new String[]{"AAAAAA"}) ) );
					validation_sample2barcode.put("sample2", new TreeSet<String>( Arrays.asList(new String[]{"TTTTTT"}) ) );
					validation_sample2barcode.put("sample3", new TreeSet<String>( Arrays.asList(new String[]{"CCCCCC"}) ) );
					validation_sample2barcode.put("sample4", new TreeSet<String>( Arrays.asList(new String[]{"GGGGGG"}) ) );
				}else if(i==1){
					validation_sample2barcode.put("sample1", new TreeSet<String>( Arrays.asList(new String[]{"AAAAAA", "TTTTTT"}) ) );
					validation_sample2barcode.put("sample2", new TreeSet<String>( Arrays.asList(new String[]{"CCCCCC","GGGGGG","GTGTGT"}) ) );
					validation_sample2barcode.put("sample3", new TreeSet<String>( Arrays.asList(new String[]{"ATGCAG","TATATA","CCGGTT"}) ) );
					
				}else if(i==2){
					validation_sample2barcode.put("sample1", new TreeSet<String>( Arrays.asList(new String[]{"NNNAAAAAANN"}) ) );
					validation_sample2barcode.put("sample2", new TreeSet<String>( Arrays.asList(new String[]{"NNNTTTTTTNN"}) ) );
					validation_sample2barcode.put("sample3", new TreeSet<String>( Arrays.asList(new String[]{"NNNCCCCCCNN"}) ) );
					validation_sample2barcode.put("sample4", new TreeSet<String>( Arrays.asList(new String[]{"NNNGGGGGGNN"}) ) );
				}else{
					throw new RuntimeException("IMPLEMENT ME NOW!!");
				}

				/*
				 * single end
				 */
				jse = new Jemultiplexer();
				jse.RUNNING_PAIRED_END = false;
				jse.BARCODE_FILE = f;
				jse.BARCODE_READ_POS = BarcodePosition.READ_1;
				jse.BARCODE_FOR_SAMPLE_MATCHING = BarcodePosition.READ_1;

				//try with no lenght initialized
				jse.BCLEN_1 = 0;

				v = new BarcodeValidator(f, jse);
				errors = v.validate();
				Assert.assertEquals("Unexpected errors during validations : \n"+printErrors(errors), 0, errors.size()); 
				Assert.assertNotNull( v.bclen1);
				Assert.assertEquals("Incorrect barcode length", realblen, v.bclen1.intValue()); 
				Assert.assertFalse(v.isLongFormat);
				try {
					v.parse();
					validateBarcodeMaps(jse, v, validation_sample2barcode);
				} catch (Exception e) {
					Assert.fail(ExceptionUtil.getStackTrace(e));
				}
				
				//try with length initialized
				jse.BCLEN_1 = realblen;
				errors = v.validate();
				Assert.assertEquals("Unexpected errors during validations : \n"+printErrors(errors), 0, errors.size()); 
				Assert.assertNotNull( v.bclen1);
				Assert.assertEquals("Incorrect barcode length", realblen, v.bclen1.intValue()); 
				Assert.assertFalse(v.isLongFormat);
				try {
					v.parse();
					validateBarcodeMaps(jse, v, validation_sample2barcode);
				} catch (Exception e) {
					Assert.fail(ExceptionUtil.getStackTrace(e));
				}
				
				/*
				 * this format should also work for simple PAIRED END but not for BARCODE_FOR_SAMPLE_MATCHING = BOTH and 
				 * jse.REDUNDANT_BARCODES=false as a bc1:bc2 should be given in such situation 
				 */
				for(BarcodePosition bp : new BarcodePosition[]{BarcodePosition.READ_1, BarcodePosition.READ_2, BarcodePosition.BOTH}){
					log.debug(bp.toString());
					for(int l : new int []{0,realblen}){

						jse = new Jemultiplexer();
						jse.RUNNING_PAIRED_END = true;
						jse.BARCODE_FILE = f;
						jse.BARCODE_READ_POS =bp;
						jse.BARCODE_FOR_SAMPLE_MATCHING = bp;


						log.debug("preset BCLEN to "+l);
						if(bp != BarcodePosition.READ_2)
							jse.BCLEN_1 = l;
						if(bp != BarcodePosition.READ_1)
							jse.BCLEN_2 = l;

						v = new BarcodeValidator(f, jse);
						errors = v.validate();
						Assert.assertEquals("Unexpected errors during validations : \n"+printErrors(errors), 0, errors.size()); 
						Assert.assertFalse(v.isLongFormat); 

						if(bp != BarcodePosition.READ_2){
							Assert.assertNotNull( v.bclen1);
							Assert.assertEquals("Incorrect barcode length", realblen, v.bclen1.intValue());
						}
						if(bp != BarcodePosition.READ_1){
							Assert.assertNotNull( v.bclen2);
							Assert.assertEquals("Incorrect barcode length", realblen, v.bclen2.intValue());
						}
						
						
						try {
							v.parse();
							validateBarcodeMaps(jse, v, validation_sample2barcode);
						} catch (Exception e) {
							Assert.fail(ExceptionUtil.getStackTrace(e));
						}
					}
				}
				
				/*
				 * this should fail validation 
				 */
				jse = new Jemultiplexer();
				jse.RUNNING_PAIRED_END = true;
				jse.BARCODE_FILE = f;
				jse.BARCODE_READ_POS =BarcodePosition.BOTH;
				jse.BARCODE_FOR_SAMPLE_MATCHING = BarcodePosition.BOTH;
				jse.REDUNDANT_BARCODES= false;
				v = new BarcodeValidator(f, jse);
				try{
					errors = v.validate();
					Assert.assertTrue("Barcode syntax with ':' is a must when using both ends for sample lookup in a non redundant fashion. Should have got an error here!", errors.size() > 0); 
				}catch (RuntimeException e) {
					// this is also ok , ie some malformed files do launch exceptions
					Assert.assertTrue(e.getCause() instanceof ParsingException);
				}
				
				
			
			}
		} catch (URISyntaxException e) {
			//this is not a test failure!
			throw new RuntimeException(e);
		}


		

	}

	@Test
	public final void testExtendedCorrectFormats() {
		Jemultiplexer jse = null;
		BarcodeValidator v = null;
		List<String> errors = null;


		/*
		 * SE extended format
		 */

		try {
			File f = new File(BarcodeFileTest.class.getResource("barcodefiles/correct_barcodes_SE_multi_with_fnames.txt").toURI()) ;
			jse = new Jemultiplexer();
			jse.RUNNING_PAIRED_END = false;
			jse.BARCODE_FILE = f;
			for(int l : new int []{0,6}){
				jse.BCLEN_1 = l;
				v = new BarcodeValidator(f, jse);
				errors = v.validate();
				Assert.assertEquals("Unexpected errors during validations : \n"+printErrors(errors), 0, errors.size()); 
				Assert.assertNotNull( v.bclen1);
				Assert.assertEquals("Incorrect barcode length", 6, v.bclen1.intValue()); 
				Assert.assertTrue(v.isLongFormat);
				try {
					v.parse();
					Map<String, String> s2f = v.sample2file1;
					//sample barcode mapping
					Assert.assertEquals(2, v.sample2barcodes.get("sample1").size());
					Assert.assertTrue(v.sample2barcodes.get("sample1").contains("AAAAAA"));
					Assert.assertTrue(v.sample2barcodes.get("sample1").contains("TTTTTT"));
					Assert.assertEquals(2, v.sample2barcodes.get("sample2").size());
					Assert.assertTrue(v.sample2barcodes.get("sample2").contains("CCCCCC"));
					Assert.assertTrue(v.sample2barcodes.get("sample2").contains("GGGGGG"));
					//sample file mapping
					Assert.assertEquals(2, s2f.size());
					Assert.assertEquals("file1.txt", s2f.get("sample1"));
					Assert.assertEquals("file2.txt", s2f.get("sample2"));
				} catch (Exception e) {
					Assert.fail(ExceptionUtil.getStackTrace(e));
				}
			}

			/*
			 * but this should NOT work for any PE situation 
			 */
			jse = new Jemultiplexer();
			jse.RUNNING_PAIRED_END = true;
			jse.BARCODE_FILE = f;
			v = new BarcodeValidator(f, jse);
			errors = v.validate(); 
			Assert.assertTrue("Long format for pair end should have 4 columns, not 3", errors.size() > 0); 

		} catch (URISyntaxException e) {
			//this is not a test failure!
			throw new RuntimeException(e);
		}
		
		/*
		 * PE extended format
		 */

		try {
			File f = new File(BarcodeFileTest.class.getResource("barcodefiles/correct_barcodes_PE_multi_with_fnames.txt").toURI()) ;
			jse = new Jemultiplexer();
			jse.RUNNING_PAIRED_END = true;
			jse.BARCODE_READ_POS = BarcodePosition.BOTH;
			jse.BARCODE_FOR_SAMPLE_MATCHING = BarcodePosition.BOTH;
			jse.REDUNDANT_BARCODES = true;
			jse.BARCODE_FILE = f;
			v = new BarcodeValidator(f, jse);
			errors = v.validate();
			Assert.assertEquals("Unexpected errors during validations : \n"+printErrors(errors), 0, errors.size()); 
			Assert.assertNotNull( v.bclen1);
			Assert.assertEquals("Incorrect barcode length", 6, v.bclen1.intValue()); 
			Assert.assertNotNull( v.bclen2);
			Assert.assertEquals("Incorrect barcode length", 6, v.bclen2.intValue()); 
			Assert.assertTrue(v.isLongFormat);
			try {
				v.parse();
				Map<String, String> s2f1 = v.sample2file1;
				Map<String, String> s2f2 = v.sample2file2;
				//sample barcode mapping
				Assert.assertEquals(2, v.sample2barcodes.get("sample1").size());
				Assert.assertTrue(v.sample2barcodes.get("sample1").contains("AAAAAA"));
				Assert.assertTrue(v.sample2barcodes.get("sample1").contains("TTTTTT"));
				Assert.assertEquals(2, v.sample2barcodes.get("sample2").size());
				Assert.assertTrue(v.sample2barcodes.get("sample2").contains("CCCCCC"));
				Assert.assertTrue(v.sample2barcodes.get("sample2").contains("GGGGGG"));
				//sample file mapping
				Assert.assertEquals(2, s2f1.size());
				Assert.assertEquals(2, s2f2.size());
				Assert.assertEquals("file1_1.txt", s2f1.get("sample1"));
				Assert.assertEquals("file1_2.txt", s2f2.get("sample1"));
				Assert.assertEquals("file2_1.txt", s2f1.get("sample2"));
				Assert.assertEquals("file2_2.txt", s2f2.get("sample2"));
			} catch (Exception e) {
				Assert.fail(ExceptionUtil.getStackTrace(e));
			}
			

			/*
			 * but this should NOT work for SE
			 */
			jse = new Jemultiplexer();
			jse.RUNNING_PAIRED_END = false;
			jse.BARCODE_FILE = f;
			v = new BarcodeValidator(f, jse);
			errors = v.validate();
			Assert.assertTrue("Long format for single end should have 3 columns, not 4", errors.size() > 0); 

		} catch (URISyntaxException e) {
			//this is not a test failure!
			throw new RuntimeException(e);
		}
		
	}
	
	@Test
	public final void testBothEndsNonRedundantCase(){
		Jemultiplexer jse = null;
		BarcodeValidator v = null;
		List<String> errors = null;

		String fname = "barcodefiles/correct_barcodes_PE_both-ends_with_fnames.txt";
		log.debug("testing with BF file : "+fname);
		try {
			File f = new File(BarcodeFileTest.class.getResource(fname).toURI()) ;
			jse = new Jemultiplexer();
			jse.RUNNING_PAIRED_END = true;
			jse.BARCODE_READ_POS = BarcodePosition.BOTH;
			jse.BARCODE_FOR_SAMPLE_MATCHING = BarcodePosition.BOTH;
			jse.REDUNDANT_BARCODES = false;
			jse.BARCODE_FILE = f;
			v = new BarcodeValidator(f, jse);
			log.debug("validating...");			
			errors = v.validate();
			log.debug("validating..done..checking values now");
			Assert.assertEquals("Unexpected errors during validations : \n"+printErrors(errors), 0, errors.size()); 
			Assert.assertNotNull( v.bclen1);
			Assert.assertEquals("Incorrect barcode length", 6, v.bclen1.intValue()); 
			Assert.assertNotNull( v.bclen2);
			Assert.assertEquals("Incorrect barcode length", 7, v.bclen2.intValue()); 
			Assert.assertTrue(v.isLongFormat);
			try {
				v.parse();
				Map<String, String> s2f1 = v.sample2file1;
				Map<String, String> s2f2 = v.sample2file2;
				
				//sample barcode mapping
				Assert.assertEquals(2, v.sample2barcodes.get("sample1").size());
				Assert.assertTrue(v.sample2barcodes.get("sample1").contains("AAAAAAATATATT"));
				Assert.assertTrue(v.sample2barcodes.get("sample1").contains("TTTTTTCGCGCGG"));
				Assert.assertEquals("sample1", v.barcode2sample.get("AAAAAAATATATT"));
				Assert.assertEquals("sample1", v.barcode2sample.get("TTTTTTCGCGCGG"));
				
				Assert.assertEquals(4, v.sample2barcodes.get("sample2").size());
				Assert.assertTrue(v.sample2barcodes.get("sample2").contains("CCCCCCAAATTTC"));
				Assert.assertTrue(v.sample2barcodes.get("sample2").contains("GGGGGGAAATTTC"));
				Assert.assertTrue(v.sample2barcodes.get("sample2").contains("CCCCCCTTTAAAG"));
				Assert.assertTrue(v.sample2barcodes.get("sample2").contains("GGGGGGTTTAAAG"));
				Assert.assertEquals("sample2", v.barcode2sample.get("CCCCCCAAATTTC"));
				Assert.assertEquals("sample2", v.barcode2sample.get("GGGGGGAAATTTC"));
				Assert.assertEquals("sample2", v.barcode2sample.get("CCCCCCTTTAAAG"));
				Assert.assertEquals("sample2", v.barcode2sample.get("GGGGGGTTTAAAG"));
				
				
				//sample file mapping
				Assert.assertEquals(2, s2f1.size());
				Assert.assertEquals(2, s2f2.size());
				Assert.assertEquals("file1_1.txt", s2f1.get("sample1"));
				Assert.assertEquals("file1_2.txt", s2f2.get("sample1"));
				Assert.assertEquals("file2_1.txt", s2f1.get("sample2"));
				Assert.assertEquals("file2_2.txt", s2f2.get("sample2"));
			} catch (Exception e) {
				Assert.fail(ExceptionUtil.getStackTrace(e));
			}
			

		} catch (URISyntaxException e) {
			//this is not a test failure!
			throw new RuntimeException(e);
		}
	}
	
	/**
	 * Checks if v has its map content 'sample2barcodes' correctly set according to the given expected mapping 
	 * 
	 * @param jse the {@link Jemultiplexer} given to the {@link BarcodeValidator} v 
	 * @param v the {@link BarcodeValidator} after parse() was called
	 * @param validation_sample2barcode the expected mapping
	 */
	private void validateBarcodeMaps(Jemultiplexer jse, BarcodeValidator v,
			HashMap<String, Set<String>> validation_sample2barcode) {
		TreeSet<String> allbc1 = new TreeSet<String>();
		TreeSet<String> allbc2 = new TreeSet<String>();
		for (Entry<String, Set<String>> e : validation_sample2barcode.entrySet()) {
			String spl = e.getKey();
			Set<String> bcs = e.getValue();
			
			if(!jse.RUNNING_PAIRED_END || jse.BARCODE_FOR_SAMPLE_MATCHING!=BarcodePosition.READ_2)
				allbc1.addAll(bcs);
			if(jse.RUNNING_PAIRED_END && jse.BARCODE_FOR_SAMPLE_MATCHING!=BarcodePosition.READ_1)
				allbc2.addAll(bcs);
				
			Assert.assertEquals("Wrong barcode number for "+spl, bcs.size(), v.sample2barcodes.get(spl).size());
			for (String bc : bcs) {
				Assert.assertTrue("Sample "+spl+" misses the barcode "+bc, v.sample2barcodes.get(spl).contains(bc));
			}
		}
		
		if(allbc1.size()>0){
			Assert.assertNotNull("barcodeSetRead1 is null while RUNNING_PAIRED_END="+jse.RUNNING_PAIRED_END
					+" , BARCODE_READ_POS="+jse.BARCODE_READ_POS
					+" , BARCODE_FOR_SAMPLE_MATCHING="+jse.BARCODE_FOR_SAMPLE_MATCHING
					+" , REDUNDANT_BARCODES="+jse.REDUNDANT_BARCODES,
					v.barcodeSetRead1);
			Assert.assertTrue("Wrong bc set size for READ1", allbc1.size() == v.barcodeSetRead1.size());
		}

		if(allbc2.size()>0){
			Assert.assertNotNull("barcodeSetRead2 is null while RUNNING_PAIRED_END="+jse.RUNNING_PAIRED_END
					+" , BARCODE_READ_POS="+jse.BARCODE_READ_POS
					+" , BARCODE_FOR_SAMPLE_MATCHING="+jse.BARCODE_FOR_SAMPLE_MATCHING
					+" , REDUNDANT_BARCODES="+jse.REDUNDANT_BARCODES,
					v.barcodeSetRead2);
			Assert.assertTrue("Wrong bc set size for READ2", allbc2.size() == v.barcodeSetRead2.size());
		}
	}

	private String printErrors(List<String> errors) {
		return StringUtil.mergeList(errors, "\n");
	}

	
	
	
}
