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
package org.embl.gbcs.je.jedropseq;

import java.io.File;
import java.net.URISyntaxException;

import org.apache.log4j.Logger;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

public class JedropseqTest {

	

	private static Logger log = Logger.getLogger(JedropseqTest.class);


	static final String READ1_FILENAME = "dropseq_test_read1.txt";
	static final String READ2_FILENAME = "dropseq_test_read2.txt";
	static final String RES_FILENAME = "dropseq_test_read1_tagged.txt";
	
	//@After
	public void cleanUpResultFiles(){
		try {
			File f1 = new File(JedropseqTest.class.getResource(READ1_FILENAME).toURI());
			

			File outdir = f1.getParentFile();
			String [] fnames = new String []{
					RES_FILENAME
			};

			for (String fname : fnames) {
				
				File f = new File(outdir, fname) ;
				
				if(f.exists()){
					f.delete();
				}
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
	public void testClipper(){
		try {
			File f1 = new File(JedropseqTest.class.getResource(READ1_FILENAME).toURI());
			File f2 = new File(JedropseqTest.class.getResource(READ2_FILENAME).toURI());
			
			File res = new File(f1.getParent(), RES_FILENAME);
			//File res = new File( RES_FILENAME);

			String[] argv = new String[] {
					"F1="+f1.getAbsolutePath(), 
					"F2="+f2.getAbsolutePath(),
					"LEN=12",
					"ULEN=8",
					"O="+ res.getAbsolutePath(),
					"GZ=false"
					
			};

			Jedropseq j = new Jedropseq();
			Assert.assertFalse(j.TEST_MODE_STOP_AFTER_PARSING);
			
			//parse
			Assert.assertEquals(0, j.instanceMain(argv));

			
			
		} catch (URISyntaxException e) {
			throw new RuntimeException(e);
		}

	}
	
}
