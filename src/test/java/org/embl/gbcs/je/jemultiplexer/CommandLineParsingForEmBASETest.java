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
import java.util.HashMap;
import java.util.Map;

import org.embl.cg.utilitytools.utils.ExceptionUtil;
import org.embl.cg.utilitytools.utils.FileUtil;
import org.embl.cg.utilitytools.utils.parser.csv.CSVLine;
import org.embl.gbcs.embase.api.EmBaseDatabaseFacade;
import org.embl.gbcs.embase.api.exception.EmBASEConnectionException;
import org.embl.gbcs.embase.api.model.NGSLibrary;
import org.junit.Assert;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class CommandLineParsingForEmBASETest {


	private static Logger log = LoggerFactory.getLogger(CommandLineParsingForEmBASETest.class);

	//connection params for testing
	public static final String DB_USER = "baseread";
	public static final String DB_PWD = "qb47ue92";
	public static final String DB_URL = "jdbc:mysql://gbcs/base"; 


	@Test
	public final void testSE() {

		try {
			String fileSE = "/g/furlong/incoming/2013-05-10-C20R8ACXX/C20R8ACXX_mesoRNAseq_test_13s003363-1-2_Schor_lane713s003363_1_sequence.txt.gz";
			Map<String, NGSLibrary> libs = new HashMap<String, NGSLibrary>();
			int id = 5569;
			libs.put("RNAseq_test_S3_68h_1", new NGSLibrary(id++, "RNAseq_test_S3_68h_1", "AAGTGC") );
			libs.put("RNAseq_test_S3_68h_2", new NGSLibrary(id++, "RNAseq_test_S3_68h_2", "TTGCGG") );
			libs.put("RNAseq_test_S3_68h_3", new NGSLibrary(id++, "RNAseq_test_S3_68h_3", "ATTATA") );
			libs.put("RNAseq_test_S3_68h_4", new NGSLibrary(id++, "RNAseq_test_S3_68h_4", "TACAAG") );
			libs.put("RNAseq_test_S6_68h_1", new NGSLibrary(id++, "RNAseq_test_S6_68h_1", "CCACTC") );
			libs.put("RNAseq_test_S6_68h_2", new NGSLibrary(id++, "RNAseq_test_S6_68h_2", "CCGTAT") );
			libs.put("RNAseq_test_S6_68h_3", new NGSLibrary(id++, "RNAseq_test_S6_68h_3", "GATGCT") );
			libs.put("RNAseq_test_U6_68h_1", new NGSLibrary(id++, "RNAseq_test_U6_68h_1", "GGAGAA") );
			libs.put("RNAseq_test_U6_68h_2", new NGSLibrary(id++, "RNAseq_test_U6_68h_2", "TCCGTC") );
			libs.put("RNAseq_test_U6_68h_3", new NGSLibrary(id++, "RNAseq_test_U6_68h_3", "CGTACG") );

			
			String[] argv = new String[] {
					"F1="+fileSE, 
					"USE_EMBASE=true"
			};

			
			EmBaseDatabaseFacade f = null;
			
			try{
				f = EmBaseDatabaseFacade.createEmBaseDatabaseFacade(null, DB_USER, DB_USER, DB_PWD);
			}catch (EmBASEConnectionException e) {
				log.info("No EmBASE connection available, skipping all EMBASE tests");
				return;
			}
			
			Jemultiplexer j = new Jemultiplexer();
			j.TEST_MODE_STOP_AFTER_PARSING = true;
			//parse
			Assert.assertEquals(0, j.instanceMain(argv));

			//check j values are correctly positioned
			Assert.assertTrue(!j.RUNNING_PAIRED_END);
			
			//this is always true in EMBASE mode
			Assert.assertTrue(j.CREATE_MD5_FILE);
			//this is always true in EMBASE mode			
			Assert.assertTrue(j.KEEP_UNASSIGNED_READ);
			//this is always true in EMBASE mode			
			Assert.assertTrue(j.GZIP_OUTPUTS);
			//this is always false in EMBASE mode			
			Assert.assertTrue(!j.STATS_ONLY);
			
			//check barcode file content
			File bc = j.BARCODE_FILE;
			log.info("BARCODE_FILE="+j.BARCODE_FILE.getAbsolutePath());
			for(CSVLine l : FileUtil.loadLinesFromCSVHeadedFile(bc, "\t", false)){
				Assert.assertEquals(3, l.getColNum());
				String n = l.getValue(0);
				String bar = l.getValue(1);
				String filep = l.getValue(2);
				
				Assert.assertTrue( libs.containsKey(n));
				NGSLibrary lib = libs.get(n);
				Assert.assertEquals( lib.getBarcode(), bar);
				lib.setSampleDir( EmBaseDatabaseFacade.getFakeSampleDir(lib.getName(), lib.getId()) );
				
				f.setUsername("girardot");
				
			}
			
			log.info("OUTPUT_DIR="+j.OUTPUT_DIR.getAbsolutePath());
			log.info("UNASSIGNED_FILE_NAME_1="+j.UNASSIGNED_FILE_NAME_1);
			if(j.RUNNING_PAIRED_END)
				log.info("UNASSIGNED_FILE_NAME_2="+j.UNASSIGNED_FILE_NAME_2);
			log.info("METRICS_FILE_NAME="+j.METRICS_FILE_NAME);
			
			
		} catch (Exception e) {
			Assert.fail("Should not have thrown exception "+ExceptionUtil.getStackTrace(e));
		}
	}
	
	
	
}