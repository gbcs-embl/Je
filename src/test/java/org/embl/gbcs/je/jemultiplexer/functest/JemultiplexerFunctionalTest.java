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

import org.embl.gbcs.je.jemultiplexer.Jemultiplexer;
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
public class JemultiplexerFunctionalTest extends BaseTesterForJemultiplexer{

	private static Logger log = LoggerFactory.getLogger(JemultiplexerFunctionalTest.class);


	@After
	public void cleanUpResultFiles(){
		try {
			File f1 = new File(JemultiplexerFunctionalTest.class.getResource("file_1_sequence.txt").toURI());

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
	public void testNoIndexFilesSE(){

		try {
			File bcFile = new File(JemultiplexerFunctionalTest.class.getResource("barcodes_SE.txt").toURI()) ;
			File f1 = new File(JemultiplexerFunctionalTest.class.getResource("file_1_sequence.txt").toURI());


			File outdir = f1.getParentFile();

			String[] argv = new String[] {
					"F1="+f1.getAbsolutePath(), 
					"BF="+bcFile.getAbsolutePath(),
					"XT=0", //XT must be turned to 0 
					"Q=0", //make sure quality is not considered
					"O="+outdir.getAbsolutePath(), 
					"GZ=false",
					"UN=false"
					,"SAME_HEADERS=false"
					,"ASYNC=false"
			};

			Jemultiplexer j = new Jemultiplexer();

			//parse
			Assert.assertEquals(0, j.instanceMain(argv));

			//System.out.println(j.OUTPUT_DIR.getAbsolutePath());

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
	public void testNoIndexFilesPERead1Only(){

		//bc in read 1 only

		try {
			File bcFile = new File(JemultiplexerFunctionalTest.class.getResource("barcodes_PE_read1only.txt").toURI()) ;
			File f1 = new File(JemultiplexerFunctionalTest.class.getResource("file_1_sequence.txt").toURI());
			File f2 = new File(JemultiplexerFunctionalTest.class.getResource("file_2_sequence.txt").toURI());

			File outdir = f1.getParentFile();

			String[] argv = new String[] {
					"F1="+f1.getAbsolutePath(), 
					"F2="+f2.getAbsolutePath(), 
					"BF="+bcFile.getAbsolutePath(),
					"XT=0", //XT must be turned to 0 
					"Q=0", //make sure quality is not considered
					"O="+outdir.getAbsolutePath(), 
					"GZ=false",
					"BPOS=READ_1",
					"UN=false"
					,"SAME_HEADERS=false"
					,"ASYNC=false"
			};

			Jemultiplexer j = new Jemultiplexer();

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
				Assert.assertTrue("Wrong barcode : "+h, h.equals("CACTGT"));
			}
			//since a single BC was used, the same should be found in PE files
			m = fetchBarcodesInHeader(new File(outdir, "sample1_2_PE.txt")); //CACTGT is the BC
			for (String h : m.values()) {
				Assert.assertTrue("Wrong barcode : "+h, h.equals("CACTGT"));
			}

			
			cleanUpResultFiles(); 
			
			// SAME headers
			argv = new String[] {
					"F1="+f1.getAbsolutePath(), 
					"F2="+f2.getAbsolutePath(), 
					"BF="+bcFile.getAbsolutePath(),
					"XT=0", //XT must be turned to 0 
					"Q=0", //make sure quality is not considered
					"O="+outdir.getAbsolutePath(), 
					"GZ=false",
					"BPOS=READ_1",
					"UN=false"
					,"ASYNC=false"
					,"SAME_HEADERS=true"
			};

			j = new Jemultiplexer();

			//parse
			Assert.assertEquals(0, j.instanceMain(argv));


			//check header format
			Map<String, String[]> ma = fetchBarcodesInSameHeader(new File(outdir, "sample1_1_PE.txt")); //CACTGT is the BC
			for (String[] h : ma.values()) {
				Assert.assertTrue("Wrong barcode length : "+h.length, h.length == 1);
				Assert.assertTrue("Wrong barcode", h[0].equals("CACTGT"));
			}
			//since a single BC was used, the same should be found in PE files
			ma = fetchBarcodesInSameHeader(new File(outdir, "sample1_2_PE.txt")); //CACTGT is the BC
			for (String[] h : ma.values()) {
				Assert.assertTrue("Wrong barcode length: "+h.length, h.length == 1);
				Assert.assertTrue("Wrong barcode", h[0].equals("CACTGT"));
			}
			
			
		} catch (URISyntaxException e) {
			throw new RuntimeException(e);
		}

	}


	@Test
	public void testNoIndexFilesPERead2Only(){

		//bc in read 2 only

		try {
			File bcFile = new File(JemultiplexerFunctionalTest.class.getResource("barcodes_PE_read2only.txt").toURI()) ;
			File f1 = new File(JemultiplexerFunctionalTest.class.getResource("file_1_sequence.txt").toURI());
			File f2 = new File(JemultiplexerFunctionalTest.class.getResource("file_2_sequence.txt").toURI());

			File outdir = f1.getParentFile();

			String[] argv = new String[] {
					"F1="+f1.getAbsolutePath(), 
					"F2="+f2.getAbsolutePath(), 
					"BF="+bcFile.getAbsolutePath(),
					"XT=0", //XT must be turned to 0 
					"Q=0", //make sure quality is not considered
					"O="+outdir.getAbsolutePath(), 
					"GZ=false",
					"BPOS=READ_2",
					"UN=false"
					,"SAME_HEADERS=false"
					,"ASYNC=false"
			};

			Jemultiplexer j = new Jemultiplexer();

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
			Map<String, String> m = fetchBarcodesInHeader(new File(outdir, "sample1_1_PE.txt")); //GTATAG is the BC from READ_2
			for (String h : m.values()) {
				Assert.assertTrue("Wrong barcode", h.equals("GTATAG"));
			}
			//since a single BC was used, the same should be found in PE files
			m = fetchBarcodesInHeader(new File(outdir, "sample1_2_PE.txt")); //GTATAG is the BC
			for (String h : m.values()) {
				Assert.assertTrue("Wrong barcode", h.equals("GTATAG"));
			}

		} catch (URISyntaxException e) {
			throw new RuntimeException(e);
		}

	}


	@Test
	public void testNoIndexFilesPEBothRead(){

		//bc in read 2 only

		try {
			File bcFile = new File(JemultiplexerFunctionalTest.class.getResource("barcodes_PE.txt").toURI()) ;
			File f1 = new File(JemultiplexerFunctionalTest.class.getResource("file_1_sequence.txt").toURI());
			File f2 = new File(JemultiplexerFunctionalTest.class.getResource("file_2_sequence.txt").toURI());

			File outdir = f1.getParentFile();

			String[] argv = new String[] {
					"F1="+f1.getAbsolutePath(), 
					"F2="+f2.getAbsolutePath(), 
					"BF="+bcFile.getAbsolutePath(),
					"XT=0", //XT must be turned to 0 
					"Q=0", //make sure quality is not considered
					"O="+outdir.getAbsolutePath(), 
					"GZ=false",
					"BPOS=BOTH",
					"BM=BOTH",
					"BRED=false",
					"UN=false"
					,"SAME_HEADERS=false"
					,"ASYNC=false"

			};

			Jemultiplexer j = new Jemultiplexer();

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
			//CACTGT:GTATAG is the combined BC used but each reads receives its own BC
			Map<String, String> m = fetchBarcodesInHeader(new File(outdir, "sample1_1_PE.txt")); 
			for (String h : m.values()) {
				Assert.assertTrue("Wrong barcode", h.equals("CACTGT"));
			}
			//the same shoudl be found in PE files
			m = fetchBarcodesInHeader(new File(outdir, "sample1_2_PE.txt")); 
			for (String h : m.values()) {
				Assert.assertTrue("Wrong barcode", h.equals("GTATAG"));
			}

		} catch (URISyntaxException e) {
			throw new RuntimeException(e);
		}


	}
	
	@Test
	public void testNoIndexFilesPEMix(){

		//bc in read 2 only

		try {
			File bcFile = new File(JemultiplexerFunctionalTest.class.getResource("barcodes_PE_read1only.txt").toURI()) ;
			File f1 = new File(JemultiplexerFunctionalTest.class.getResource("file_1_sequence.txt").toURI());
			File f2 = new File(JemultiplexerFunctionalTest.class.getResource("file_2_sequence.txt").toURI());

			File outdir = f1.getParentFile();

			String[] argv = new String[] {
					"F1="+f1.getAbsolutePath(), 
					"F2="+f2.getAbsolutePath(), 
					"BF="+bcFile.getAbsolutePath(),
					"XT=0", //XT must be turned to 0 
					"Q=0", //make sure quality is not considered
					"O="+outdir.getAbsolutePath(), 
					"GZ=false",
					"BPOS=BOTH",
					"BM=READ_1",
					"BCLEN=6:10",
					"UN=false"
					,"SAME_HEADERS=false"
					,"ASYNC=false"

			};

			Jemultiplexer j = new Jemultiplexer();

			//parse
			Assert.assertEquals(0, j.instanceMain(argv));


			//check result file content, no need to check them all (if it works for one BC, it works for the other one as well)
			Assert.assertTrue("Dmltplexd file sample1_1_PE.txt does not contain expected reads", 
					checkReadsInFile(outdir, "sample1_1_PE.txt", "expout/CACTGT_1_noclip_top6reads.txt") );
			
			Assert.assertTrue("Dmltplexd file sample1_1_PE.txt does not contain reads of expected length", 
					checkReadLength(new File(outdir, "sample1_1_PE.txt"), 94) );
			
			
			Assert.assertTrue("Dmltplexd file sample1_2_PE.txt does not contain expected reads", 
					checkReadsInFile(outdir, "sample1_2_PE.txt", "expout/GTATAG_2_noclip_top6-CACTGT-reads.txt") );

			Assert.assertTrue("Dmltplexd file sample1_2_PE.txt does not contain reads of expected length", 
					checkReadLength(new File(outdir, "sample1_2_PE.txt"), 90) );
			
			
			//check header format
			//only read 1 is used for BM but another 10 base long bc is defined on read 2
			Map<String, String> m = fetchBarcodesInHeader(new File(outdir, "sample1_1_PE.txt")); 
			for (String h : m.values()) {
				//here the 6 base long sample encoding BC is found
				Assert.assertTrue("Wrong barcode", h.equals("CACTGT"));
			}
			//a 10 base long sequence should be found on this side, the first 6 bases should be GTATAG
			m = fetchBarcodesInHeader(new File(outdir, "sample1_2_PE.txt")); 
			for (String h : m.values()) {
				Assert.assertTrue("Wrong barcode", h.startsWith("GTATAG"));
				Assert.assertTrue("Wrong barcode length", h.length() == 10);
			}


		} catch (URISyntaxException e) {
			throw new RuntimeException(e);
		}


	}
	
	
	@Test
	public void testNoIndexFilesPEMix2(){

		//bc in read 2 only

		try {
			File bcFile = new File(JemultiplexerFunctionalTest.class.getResource("barcodes_PE_read1only.txt").toURI()) ;
			File f1 = new File(JemultiplexerFunctionalTest.class.getResource("file_1_sequence.txt").toURI());
			File f2 = new File(JemultiplexerFunctionalTest.class.getResource("file_2_sequence.txt").toURI());

			File outdir = f1.getParentFile();

			String[] argv = new String[] {
					"F1="+f1.getAbsolutePath(), 
					"F2="+f2.getAbsolutePath(), 
					"BF="+bcFile.getAbsolutePath(),
					"XT=1:0",  
					"Q=0", //make sure quality is not considered
					"O="+outdir.getAbsolutePath(), 
					"GZ=false",
					"BPOS=BOTH",
					"BM=READ_1",
					"BCLEN=6:10",
					"ZT=8:5",
					"UN=false"
					,"SAME_HEADERS=false"
					,"ASYNC=false"

			};

			Jemultiplexer j = new Jemultiplexer();

			//parse
			Assert.assertEquals(0, j.instanceMain(argv));


			//check result file content, no need to check them all (if it works for one BC, it works for the other one as well)
			Assert.assertTrue("Dmltplexd file sample1_1_PE.txt does not contain expected reads", 
					checkReadsInFile(outdir, "sample1_1_PE.txt", "expout/CACTGT_1_noclip_top6reads.txt") );
			
			Assert.assertTrue("Dmltplexd file sample1_1_PE.txt does not contain reads of expected length", 
					checkReadLength(new File(outdir, "sample1_1_PE.txt"), 85) );
			
			
			Assert.assertTrue("Dmltplexd file sample1_2_PE.txt does not contain expected reads", 
					checkReadsInFile(outdir, "sample1_2_PE.txt", "expout/GTATAG_2_noclip_top6-CACTGT-reads.txt") );

			Assert.assertTrue("Dmltplexd file sample1_2_PE.txt does not contain reads of expected length", 
					checkReadLength(new File(outdir, "sample1_2_PE.txt"), 85) );
			
			
			//check header format ; XT and ZT optiosn have no influence on barcodes in headers
			//only read 1 is used for BM but another 10 base long bc is defined on read 2
			Map<String, String> m = fetchBarcodesInHeader(new File(outdir, "sample1_1_PE.txt")); 
			for (String h : m.values()) {
				//here the 6 base long sample encoding BC is found
				Assert.assertTrue("Wrong barcode", h.equals("CACTGT"));
			}
			//a 10 base long sequence should be found on this side, the first 6 bases should be GTATAG
			m = fetchBarcodesInHeader(new File(outdir, "sample1_2_PE.txt")); 
			for (String h : m.values()) {
				Assert.assertTrue("Wrong barcode", h.startsWith("GTATAG"));
				Assert.assertTrue("Wrong barcode length", h.length() == 10);
			}



		} catch (URISyntaxException e) {
			throw new RuntimeException(e);
		}


	}





}














