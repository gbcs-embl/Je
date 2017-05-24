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
package org.embl.gbcs.je.jemultiplexer.functest;

import java.io.File;
import java.net.URISyntaxException;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.embl.cg.utilitytools.utils.StringUtil;
import org.embl.gbcs.je.jemultiplexer.BarcodePosition;
import org.embl.gbcs.je.jemultiplexer.JemultiplexerIllumina;
import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;





/**
 * This class launches Jemultiplexer in a number of conditions (SE, PE) 
 * and checks that read have be distributed as they should 
 * 
 * @author girardot
 *
 */
public class JemultiplexerIlluminaFunctionalTest extends BaseTesterForJemultiplexer{

	private static Logger log = LoggerFactory.getLogger(JemultiplexerIlluminaFunctionalTest.class);


	@After
	public void cleanUpResultFiles(){
		try {
			File f1 = new File(JemultiplexerIlluminaFunctionalTest.class.getResource("file_1_sequence.txt").toURI());

			File outdir = f1.getParentFile();
			String [] fnames = new String []{
					"sample1_SE.txt", "sample2_SE.txt", 
					"sample3_SE.txt", "sample4_SE.txt",
					"sample1_1_PE.txt","sample1_2_PE.txt",
					"sample2_1_PE.txt","sample2_2_PE.txt",
					"sample3_1_PE.txt","sample3_2_PE.txt",
					"sample4_1_PE.txt","sample4_2_PE.txt"
			};

			for (String fname : fnames) {
				
				File f = new File(outdir, fname) ;
				
				if(f.exists()){
					f.delete();
					//log.info("DELETING left over result files "+f.getAbsolutePath());
				}
				//else{ log.info("Ignoring "+f.getAbsolutePath());}
			}
			log.info("Test Env. Cleanup done");
		} catch (URISyntaxException e) {
			throw new RuntimeException("Test Env. Tear Down : Cleanup failed");
		}
	}
	
	@Before
	public void cleanUpLeftOverResultFiles(){
		cleanUpResultFiles();
	}


	@Test
	public void testWithIndexFilesAndUMISE(){

		try {
			File bcFile = new File(JemultiplexerIlluminaFunctionalTest.class.getResource("barcodes_SE.txt").toURI()) ;
			File f1 = new File(JemultiplexerIlluminaFunctionalTest.class.getResource("file_1_sequence.txt").toURI());
			File i1 = new File(JemultiplexerIlluminaFunctionalTest.class.getResource("IDX_1_sequence.txt").toURI());


			File outdir = f1.getParentFile();

			String[] argv = new String[] {
					"F1="+f1.getAbsolutePath(), 
					"I1="+i1.getAbsolutePath(),
					"BF="+bcFile.getAbsolutePath(),
					"LEN=8",
					"BPOS=READ_1", 
					"XT=0", 
					"Q=0", //make sure quality is not considered
					"O="+outdir.getAbsolutePath(), 
					"GZ=false",
					"UN=false",
					"ASYNC=false",
					"READ_NAME_REPLACE_CHAR=:"
			};

			JemultiplexerIllumina j = new JemultiplexerIllumina();

			//parse
			Assert.assertEquals(0, j.instanceMain(argv));
			Assert.assertTrue(j.ADD_BARCODE_TO_HEADER);
			Assert.assertTrue(j.CLIP_BARCODE);
			Assert.assertTrue(j.BARCODE_READ_POS == BarcodePosition.READ_1);
			Assert.assertNotNull(j.READ_NAME_REPLACE_CHAR);

			
			//check header format
			Map<String, String> m = fetchBarcodesInHeader(new File(outdir, "sample1_SE.txt")); //CACTGT is the BC
			for (String h : m.values()) {
				log.debug(h);
				Assert.assertTrue("Wrong sample barcode", h.contains("CACTGT"));
				Assert.assertTrue("Wrong UMI length", h.replace("CACTGT:","").length() == 8);
			}
			
			
		} catch (URISyntaxException e) {
			throw new RuntimeException(e);
		}


	}

	
	@Test
	public void testWithIndexFilesSE(){

		try {
			File bcFile = new File(JemultiplexerIlluminaFunctionalTest.class.getResource("barcodes_SE.txt").toURI()) ;
			File f1 = new File(JemultiplexerIlluminaFunctionalTest.class.getResource("file_1_sequence.txt").toURI());
			File i1 = new File(JemultiplexerIlluminaFunctionalTest.class.getResource("IDX_1_sequence.txt").toURI());


			File outdir = f1.getParentFile();

			String[] argv = new String[] {
					"F1="+f1.getAbsolutePath(), 
					"I1="+i1.getAbsolutePath(),
					"BF="+bcFile.getAbsolutePath(),
					"XT=0", 
					"Q=0", //make sure quality is not considered
					"O="+outdir.getAbsolutePath(), 
					"GZ=false",
					"UN=false",
					"ASYNC=false",
					"ENSURE_IDENTICAL_HEADER_NAMES=false"
			};

			JemultiplexerIllumina j = new JemultiplexerIllumina();

			//parse
			Assert.assertEquals(0, j.instanceMain(argv));

			//check result file content
			Assert.assertTrue("Dmltplexd file sample1_SE.txt does not contain expected reads", 
					checkReadsInFile(outdir, "sample1_SE.txt", "expout/CACTGT_1_noclip_top6reads.txt") );
			Assert.assertTrue("Dmltplexd file sample1_SE.txt does not contain expected reads", 
					checkReadsInFile(outdir,  "sample2_SE.txt", "expout/ATTCCG_1_noclip_7reads.txt") );
			Assert.assertTrue("Dmltplexd file sample1_SE.txt does not contain expected reads", 
					checkReadsInFile(outdir, "sample3_SE.txt", "expout/GCTACC_1_noclip_6reads.txt") );
			Assert.assertTrue("Dmltplexd file sample1_SE.txt does not contain expected reads", 
					checkReadsInFile(outdir, "sample4_SE.txt", "expout/CGAAAC_1_noclip_last6reads.txt") );

			//check header format
			Map<String, String> m = fetchBarcodesInHeader(new File(outdir, "sample1_SE.txt")); //CACTGT is the BC
			for (String h : m.values()) {
				Assert.assertTrue("Wrong barcode", h.equals("CACTGT"));
			}
			
			
		} catch (URISyntaxException e) {
			throw new RuntimeException(e);
		}


	}

	
	@Test
	public void testWithIndexFilesPEIDX1Only(){

		//bc in read 1 only

		try {
			File bcFile = new File(JemultiplexerIlluminaFunctionalTest.class.getResource("barcodes_PE_read1only.txt").toURI()) ;
			File f1 = new File(JemultiplexerIlluminaFunctionalTest.class.getResource("file_1_sequence.txt").toURI());
			File f2 = new File(JemultiplexerIlluminaFunctionalTest.class.getResource("file_2_sequence.txt").toURI());
			File i1 = new File(JemultiplexerIlluminaFunctionalTest.class.getResource("IDX_1_sequence.txt").toURI());

			File outdir = f1.getParentFile();

			String[] argv = new String[] {
					"F1="+f1.getAbsolutePath(), 
					"F2="+f2.getAbsolutePath(), 
					"I1="+i1.getAbsolutePath(),
					"BF="+bcFile.getAbsolutePath(),
					"XT=0", //XT must be turned to 0 
					"Q=0", //make sure quality is not considered
					"O="+outdir.getAbsolutePath(), 
					"GZ=false",
					"UN=false"
					,"ASYNC=false"
					,"CLIP_BARCODE=false"
					,"SAME_HEADERS=false"
					//,"CREATE_MD5_FILE=true"
			};

			JemultiplexerIllumina j = new JemultiplexerIllumina();

			//parse
			Assert.assertEquals(0, j.instanceMain(argv));


			//check result file content
			Assert.assertTrue("Dmltplexd file sample1_1_PE.txt does not contain expected reads", 
					checkReadsInFile(outdir, "sample1_1_PE.txt", "expout/CACTGT_1_noclip_top6reads.txt") );
			Assert.assertTrue("Dmltplexd file sample1_2_PE.txt does not contain expected reads", 
					checkReadsInFile(outdir, "sample1_2_PE.txt", "expout/GTATAG_2_noclip_top6-CACTGT-reads.txt") );


			Assert.assertTrue("Dmltplexd file sample2_1_PE.txt does not contain expected reads", 
					checkReadsInFile(outdir,  "sample2_1_PE.txt", "expout/ATTCCG_1_noclip_7reads.txt") );
			Assert.assertTrue("Dmltplexd file sample2_2_PE.txt does not contain expected reads", 
					checkReadsInFile(outdir, "sample2_2_PE.txt", "expout/TCCGTC_2_noclip_next7-ATTCCG-reads.txt") );


			Assert.assertTrue("Dmltplexd file sample3_1_PE.txt does not contain expected reads", 
					checkReadsInFile(outdir, "sample3_1_PE.txt", "expout/GCTACC_1_noclip_6reads.txt") );
			Assert.assertTrue("Dmltplexd file sample3_2_PE.txt does not contain expected reads", 
					checkReadsInFile(outdir, "sample3_2_PE.txt", "expout/TGGTCA_2_noclip_next6-GCTACC-reads.txt") );


			Assert.assertTrue("Dmltplexd file sample4_1_PE.txt does not contain expected reads", 
					checkReadsInFile(outdir, "sample4_1_PE.txt", "expout/CGAAAC_1_noclip_last6reads.txt") );
			Assert.assertTrue("Dmltplexd file sample4_2_PE.txt does not contain expected reads", 
					checkReadsInFile(outdir, "sample4_2_PE.txt", "expout/CACTGT_2_noclip_last6-CGAAAC-reads.txt") );


			//check header format
			Map<String, String> m = fetchBarcodesInHeader(new File(outdir, "sample1_1_PE.txt")); //CACTGT is the BC
			for (String h : m.values()) {
				Assert.assertTrue("Wrong barcode", h.equals("CACTGT"));
			}
			//since a single BC was used, the same shoudl be found in PE files
			m = fetchBarcodesInHeader(new File(outdir, "sample1_2_PE.txt")); //CACTGT is the BC
			for (String h : m.values()) {
				Assert.assertTrue("Wrong barcode", h.equals("CACTGT"));
			}

		} catch (URISyntaxException e) {
			throw new RuntimeException(e);
		}

	}

	@Test
	public void testWithIndexFilesPEBothIDX(){

		//bc in read 1 only

		try {
			File bcFile = new File(JemultiplexerIlluminaFunctionalTest.class.getResource("barcodes_PE.txt").toURI()) ;
			File f1 = new File(JemultiplexerIlluminaFunctionalTest.class.getResource("file_1_sequence.txt").toURI());
			File f2 = new File(JemultiplexerIlluminaFunctionalTest.class.getResource("file_2_sequence.txt").toURI());
			File i1 = new File(JemultiplexerIlluminaFunctionalTest.class.getResource("IDX_1_sequence.txt").toURI());
			File i2 = new File(JemultiplexerIlluminaFunctionalTest.class.getResource("IDX_2_sequence.txt").toURI());

			File outdir = f1.getParentFile();

			String[] argv = new String[] {
					"F1="+f1.getAbsolutePath(), 
					"F2="+f2.getAbsolutePath(), 
					"I1="+i1.getAbsolutePath(),
					"I2="+i2.getAbsolutePath(),
					"BF="+bcFile.getAbsolutePath(),
					"XT=0", //XT must be turned to 0 
					"Q=0", //make sure quality is not considered
					"O="+outdir.getAbsolutePath(), 
					"BRED=false",
					"GZ=false",
					"UN=false",
					"ASYNC=false"
					,"SAME_HEADERS=false"
			};

			JemultiplexerIllumina j = new JemultiplexerIllumina();

			//parse
			Assert.assertEquals(0, j.instanceMain(argv));


			//check result file content
			Assert.assertTrue("Dmltplexd file sample1_1_PE.txt does not contain expected reads", 
					checkReadsInFile(outdir, "sample1_1_PE.txt", "expout/CACTGT_1_noclip_top6reads.txt") );
			Assert.assertTrue("Dmltplexd file sample2_1_PE.txt does not contain expected reads", 
					checkReadsInFile(outdir, "sample1_2_PE.txt", "expout/GTATAG_2_noclip_top6-CACTGT-reads.txt") );


			Assert.assertTrue("Dmltplexd file sample2_1_PE.txt does not contain expected reads", 
					checkReadsInFile(outdir,  "sample2_1_PE.txt", "expout/ATTCCG_1_noclip_7reads.txt") );
			Assert.assertTrue("Dmltplexd file sample2_2_PE.txt does not contain expected reads", 
					checkReadsInFile(outdir, "sample2_2_PE.txt", "expout/TCCGTC_2_noclip_next7-ATTCCG-reads.txt") );


			Assert.assertTrue("Dmltplexd file sample3_1_PE.txt does not contain expected reads", 
					checkReadsInFile(outdir, "sample3_1_PE.txt", "expout/GCTACC_1_noclip_6reads.txt") );
			Assert.assertTrue("Dmltplexd file sample3_2_PE.txt does not contain expected reads", 
					checkReadsInFile(outdir, "sample3_2_PE.txt", "expout/TGGTCA_2_noclip_next6-GCTACC-reads.txt") );


			Assert.assertTrue("Dmltplexd file sample4_1_PE.txt does not contain expected reads", 
					checkReadsInFile(outdir, "sample4_1_PE.txt", "expout/CGAAAC_1_noclip_last6reads.txt") );
			Assert.assertTrue("Dmltplexd file sample4_2_PE.txt does not contain expected reads", 
					checkReadsInFile(outdir, "sample4_2_PE.txt", "expout/CACTGT_2_noclip_last6-CGAAAC-reads.txt") );


			//check header format
			Map<String, String> m = fetchBarcodesInHeader(new File(outdir, "sample1_1_PE.txt")); //CACTGT:GTATAG is the combined BC
			for (String h : m.values()) {
				Assert.assertTrue("Wrong barcode", h.equals("CACTGT:GTATAG"));
			}
			//the same shoudl be found in PE files
			m = fetchBarcodesInHeader(new File(outdir, "sample1_2_PE.txt")); 
			for (String h : m.values()) {
				Assert.assertTrue("Wrong barcode", h.equals("CACTGT:GTATAG"));
			}

		} catch (URISyntaxException e) {
			throw new RuntimeException(e);
		}

	}
	
	@Test
	public void testWithIndexFilesPEBothIDXAndInternalBCs(){

		//bc in read 1 only

		try {
			File bcFile = new File(JemultiplexerIlluminaFunctionalTest.class.getResource("barcodes_PE.txt").toURI()) ;
			File f1 = new File(JemultiplexerIlluminaFunctionalTest.class.getResource("file_1_sequence.txt").toURI());
			File f2 = new File(JemultiplexerIlluminaFunctionalTest.class.getResource("file_2_sequence.txt").toURI());
			File i1 = new File(JemultiplexerIlluminaFunctionalTest.class.getResource("IDX_1_sequence.txt").toURI());
			File i2 = new File(JemultiplexerIlluminaFunctionalTest.class.getResource("IDX_2_sequence.txt").toURI());

			File outdir = f1.getParentFile();

			String[] argv = new String[] {
					"F1="+f1.getAbsolutePath(), 
					"F2="+f2.getAbsolutePath(), 
					"I1="+i1.getAbsolutePath(),
					"I2="+i2.getAbsolutePath(),
					"BF="+bcFile.getAbsolutePath(),
					"XT=0:5", 
					"ZT=0:5", 
					"Q=0", //make sure quality is not considered
					"O="+outdir.getAbsolutePath(), 
					"BPOS=READ_1",
					"BCLEN=10",
					"GZ=false",
					"UN=false",
					"ASYNC=false"
					,"SAME_HEADERS=false"
			};

			JemultiplexerIllumina j = new JemultiplexerIllumina();

			//parse
			Assert.assertEquals(0, j.instanceMain(argv));


			//check result file content
			Assert.assertTrue("Dmltplexd file sample1_1_PE.txt does not contain expected reads", 
					checkReadsInFile(outdir, "sample1_1_PE.txt", "expout/CACTGT_1_noclip_top6reads.txt") );
			Assert.assertTrue("Dmltplexd file sample1_2_PE.txt does not contain expected reads", 
					checkReadsInFile(outdir, "sample1_2_PE.txt", "expout/GTATAG_2_noclip_top6-CACTGT-reads.txt") );

			//read 1 should be 10 base less due to internal BC
			Assert.assertTrue("Dmltplexd file sample1_1_PE.txt does not contain reads of expected length", 
					checkReadLength( new File(outdir, "sample1_1_PE.txt"), 90) );
			
			
			//read 2 should be 10 base less due XT and ZT option
			Assert.assertTrue("Dmltplexd file sample1_2_PE.txt does not contain reads of expected length", 
					checkReadLength( new File(outdir, "sample1_2_PE.txt"), 90) );
			
			//no need to check the other demultiplexed files

			/*
			 * check header format :
			 * CACTGT:GTATAG is the combined BC
			 * but we declared an extra 10-base long BC in the READ_1 => file should hold an extra 10 base long sequence ; but not file 2
			 */
			Map<String, String> m = fetchBarcodesInHeader(new File(outdir, "sample1_1_PE.txt")); 
			Pattern p = Pattern.compile("^CACTGT:GTATAG:[ACTGN]{10,10}$");
			for (String h : m.values()) {
				Matcher matcher = p.matcher(h);
				Assert.assertTrue("Wrong barcode end (expected 'CACTGT:GTATAG:[ACTGN]{10,10}' ): "+h, matcher.matches());
			}
			//only combined sample BC in read 2 file
			m = fetchBarcodesInHeader(new File(outdir, "sample1_2_PE.txt")); 
			for (String h : m.values()) {
				Assert.assertTrue("Wrong barcode", h.equals("CACTGT:GTATAG"));
			}
			
			
			argv = new String[] {
					"F1="+f1.getAbsolutePath(), 
					"F2="+f2.getAbsolutePath(), 
					"I1="+i1.getAbsolutePath(),
					"I2="+i2.getAbsolutePath(),
					"BF="+bcFile.getAbsolutePath(),
					"XT=0", 
					"ZT=0", 
					"Q=0", //make sure quality is not considered
					"O="+outdir.getAbsolutePath(), 
					"BPOS=BOTH",
					"BCLEN=10",
					"GZ=false",
					"UN=false"
					,"ASYNC=false"
					,"FORCE=true"
					,"SAME_HEADERS=false"
			};

			j = new JemultiplexerIllumina();

			//parse
			Assert.assertEquals(0, j.instanceMain(argv));


			//read 1 should be 10 base less due to internal BC
			Assert.assertTrue("Dmltplexd file sample1_1_PE.txt does not contain reads of expected length", 
					checkReadLength( new File(outdir, "sample1_1_PE.txt"), 90) );
			
			
			//read 2 should be 10 base less due to internal BC
			Assert.assertTrue("Dmltplexd file sample1_2_PE.txt does not contain reads of expected length", 
					checkReadLength( new File(outdir, "sample1_2_PE.txt"), 90) );
			
			//no need to check the other demulitplexed files
			
			
			/*
			 * check header format :
			 * CACTGT:GTATAG is the combined BC
			 * but we declared an extra 10-base long BC in the both reads => files should hold an extra 10 base long sequence 
			 */
			Map<String, String> m1 = fetchBarcodesInHeader(new File(outdir, "sample1_1_PE.txt")); 
			Map<String, String> m2 = fetchBarcodesInHeader(new File(outdir, "sample1_2_PE.txt"));
			p = Pattern.compile("^CACTGT:GTATAG:([ACTGN]{10,10})$");
			for (String k : m1.keySet()) {
				String h1 = m1.get(k);
				String h2 = m2.get(k);
				Matcher matcher = p.matcher(h1);
				Assert.assertTrue("Wrong barcode end (expected 'CACTGT:GTATAG:[ACTGN]{10,10}' ): "+h1, matcher.matches());
				String xtra1 = matcher.group(1);
				
				matcher = p.matcher(h2);
				Assert.assertTrue("Wrong barcode end (expected 'CACTGT:GTATAG:[ACTGN]{10,10}' ): "+h2, matcher.matches());
				String xtra2 = matcher.group(1);
				
				Assert.assertFalse("Should not be the same barcoding sequence", xtra1.equals(xtra2));
				
			}
			
			
			
		} catch (URISyntaxException e) {
			throw new RuntimeException(e);
		}

	}
	
	
	
	@Test
	public void testWithIndexFilesPEBothIDXAndInternalBCsSameHeaders(){

		//bc in read 1 only

		try {
			File bcFile = new File(JemultiplexerIlluminaFunctionalTest.class.getResource("barcodes_PE.txt").toURI()) ;
			File f1 = new File(JemultiplexerIlluminaFunctionalTest.class.getResource("file_1_sequence.txt").toURI());
			File f2 = new File(JemultiplexerIlluminaFunctionalTest.class.getResource("file_2_sequence.txt").toURI());
			File i1 = new File(JemultiplexerIlluminaFunctionalTest.class.getResource("IDX_1_sequence.txt").toURI());
			File i2 = new File(JemultiplexerIlluminaFunctionalTest.class.getResource("IDX_2_sequence.txt").toURI());

			File outdir = f1.getParentFile();

			String[]  argv = new String[] {
					"F1="+f1.getAbsolutePath(), 
					"F2="+f2.getAbsolutePath(), 
					"I1="+i1.getAbsolutePath(),
					"I2="+i2.getAbsolutePath(),
					"BF="+bcFile.getAbsolutePath(),
					"XT=0", 
					"ZT=0", 
					"Q=0", //make sure quality is not considered
					"O="+outdir.getAbsolutePath(), 
					"BPOS=BOTH",
					"BCLEN=10",
					"GZ=false",
					"UN=false"
					,"ASYNC=false"
					,"FORCE=true"
					,"SAME_HEADERS=true"
			};

			JemultiplexerIllumina j = new JemultiplexerIllumina();

			//parse
			Assert.assertEquals(0, j.instanceMain(argv));


			//read 1 should be 10 base less due to internal BC
			Assert.assertTrue("Dmltplexd file sample1_1_PE.txt does not contain reads of expected length", 
					checkReadLength( new File(outdir, "sample1_1_PE.txt"), 90) );
			
			
			//read 2 should be 10 base less due to internal BC
			Assert.assertTrue("Dmltplexd file sample1_2_PE.txt does not contain reads of expected length", 
					checkReadLength( new File(outdir, "sample1_2_PE.txt"), 90) );
			
			//no need to check the other demulitplexed files
			
			
			/*
			 * check header format :
			 * CACTGT:GTATAG is the combined BC
			 * but we declared an extra 10-base long BC in the both reads => files should hold an extra 10 base long sequence 
			 */
			Map<String, String[]> m1 = fetchBarcodesInSameHeader(new File(outdir, "sample1_1_PE.txt")); 
			Map<String, String[]> m2 = fetchBarcodesInSameHeader(new File(outdir, "sample1_2_PE.txt"));
			Pattern p = Pattern.compile("^CACTGT:GTATAG:[ACTGN]{10,10}:[ACTGN]{10,10}$");
			Matcher matcher = p.matcher("");
			for (String k : m1.keySet()) {
				String[] h1arr = m1.get(k);
				String h1 = StringUtil.mergeArray(h1arr, ":");
				String[] h2arr = m2.get(k);
				String h2 = StringUtil.mergeArray(h2arr, ":");
				Assert.assertTrue("Note same header !: read1="+h1+" ; read2="+h2, h1.equals(h2));
				Assert.assertTrue("expected 4 entries in barcode section :"+h1, h1arr.length == 4);
				//log.debug("'"+h1+"'");
				matcher.reset(h1);
				Assert.assertTrue("Wrong barcode string (expected 'CACTGT:GTATAG:[ACTGN]{10,10}:[ACTGN]{10,10}' ): "+h1, matcher.matches());
				
			}
			
			
			//only one internal bc
			argv = new String[] {
					"F1="+f1.getAbsolutePath(), 
					"F2="+f2.getAbsolutePath(), 
					"I1="+i1.getAbsolutePath(),
					"I2="+i2.getAbsolutePath(),
					"BF="+bcFile.getAbsolutePath(),
					"XT=0", 
					"ZT=0:10", 
					"Q=0", //make sure quality is not considered
					"O="+outdir.getAbsolutePath(), 
					"BPOS=READ_1",
					"BCLEN=10",
					"GZ=false",
					"UN=false"
					,"ASYNC=false"
					,"FORCE=true"
					,"SAME_HEADERS=true"
			};

			j = new JemultiplexerIllumina();

			//parse
			Assert.assertEquals(0, j.instanceMain(argv));


			//read 1 should be 10 base less due to internal BC
			Assert.assertTrue("Dmltplexd file sample1_1_PE.txt does not contain reads of expected length", 
					checkReadLength( new File(outdir, "sample1_1_PE.txt"), 90) );
			
			
			//read 2 should be 10 base less due to internal BC
			Assert.assertTrue("Dmltplexd file sample1_2_PE.txt does not contain reads of expected length", 
					checkReadLength( new File(outdir, "sample1_2_PE.txt"), 90) );
			
			//no need to check the other demulitplexed files
			
			
			/*
			 * check header format :
			 * CACTGT:GTATAG is the combined BC
			 * but we declared an extra 10-base long BC in the both reads => files should hold an extra 10 base long sequence 
			 */
			m1 = fetchBarcodesInSameHeader(new File(outdir, "sample1_1_PE.txt")); 
			m2 = fetchBarcodesInSameHeader(new File(outdir, "sample1_2_PE.txt"));
			p = Pattern.compile("^CACTGT:GTATAG:[ACTGN]{10,10}$");
			matcher = p.matcher("");
			for (String k : m1.keySet()) {
				String[] h1arr = m1.get(k);
				String h1 = StringUtil.mergeArray(h1arr, ":");
				String[] h2arr = m2.get(k);
				String h2 = StringUtil.mergeArray(h2arr, ":");
				Assert.assertTrue("Note same header !: read1="+h1+" ; read2="+h2, h1.equals(h2));
				Assert.assertTrue("expected 3 entries in barcode section :"+h1, h1arr.length == 3);
				//log.debug("'"+h1+"'");
				matcher.reset(h1);
				Assert.assertTrue("Wrong barcode string (expected 'CACTGT:GTATAG:[ACTGN]{10,10}' ): "+h1, matcher.matches());
				
			}
			
			
			
			
		} catch (URISyntaxException e) {
			throw new RuntimeException(e);
		}

	}
	




}














