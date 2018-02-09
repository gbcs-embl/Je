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
package org.embl.gbcs.je.jeclipper;

import java.io.File;
import java.net.URISyntaxException;

import org.apache.log4j.Logger;
import org.embl.gbcs.je.jeclipper.Jeclipper;
import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

public class JeclipperTest {

	

	private static Logger log = Logger.getLogger(JeclipperTest.class);


	@After
	public void cleanUpResultFiles(){
		try {
			File f1 = new File(JeclipperTest.class.getResource("file_1_sequence.txt").toURI());

			File outdir = f1.getParentFile();
			String [] fnames = new String []{
					"out_1.txt", "out_2.txt"
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
			File f1 = new File(JeclipperTest.class.getResource("file_1_sequence.txt").toURI());
			File f2 = new File(JeclipperTest.class.getResource("file_2_sequence.txt").toURI());


			File outdir = f1.getParentFile();

			String[] argv = new String[] {
					"F="+f1.getAbsolutePath(), 
					"F="+f2.getAbsolutePath(),
					"RL=<BARCODE1:12><SAMPLE1:x>",
					"RL=<BARCODE2:8><SAMPLE2:x>",
					"OL=1:B1B2:S1",
					"OL=2:B1B2:S2",
					"O="+outdir.getAbsolutePath(), 
					"GZ=false"
					
			};

			Jeclipper j = new Jeclipper();
			Assert.assertFalse(j.TEST_MODE_STOP_AFTER_PARSING);
			
			//parse
			Assert.assertEquals(0, j.instanceMain(argv));

			
			
		} catch (URISyntaxException e) {
			throw new RuntimeException(e);
		}

	}
	
}
