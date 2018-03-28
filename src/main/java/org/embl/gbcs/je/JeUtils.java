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

import java.util.Set;
import java.util.TreeSet;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.developpez.adiguba.shell.ProcessConsumer;

import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.util.FastqQualityFormat;
import htsjdk.samtools.util.QualityEncodingDetector;
import htsjdk.samtools.util.SolexaQualityConverter;

public class JeUtils {

	
	private static Logger log = LoggerFactory.getLogger(JeUtils.class);

	
	/*
	 * BC : for sample barcode, raw or corrected, with QT to store its quality string
	 */
	public static final String SAMTAG_BC = "BC"; 
	
	/*
	 * QT : Phred quality of the sample-barcode sequence in the BC (or RT) tag
	 */
	public static final String SAMTAG_QT = "QT"; 
	
	
	/*
	 * RX : Sequence bases of the (possibly corrected) unique molecular identifier
	 */
	public static final String SAMTAG_RX = "RX"; 
	
	/*
	 * QX : Quality score of the unique molecular identifier in the RX tag
	 */
	public static final String SAMTAG_QX = "QX"; 
	
	/*
	 * OX : Original unique molecular barcode bases
	 */
	public static final String SAMTAG_OX = "OX"; 
	
	/*
	 * BZ : Phred quality of the unique molecular barcode bases in the OX tag
	 */
	public static final String SAMTAG_BZ = "BZ"; 
	
	/*
	 * MI : Molecular identifier; a string that uniquely identifies the molecule from which the record was derived
	 */
	public static final String SAMTAG_MI = "MI"; 
	
	
	/**convert a string of quality numbers in Phred Scale (each quality has 2 char) to the Standard Phred + 33 encoding
	 * ie 
	 * @param s
	 * @return
	 */
	public static String toBytesThenPhred(String s) {
		byte [] arr = new byte [s.length()/2];
		int i =0;
		for(String t : s.split("(?<=\\G.{2})")) {
			arr[i] = Byte.parseByte(t);
			i++;
		}
		return SAMUtils.phredToFastq(arr);
	}
	
	
	/** 
	 * Based on the type of quality scores coming in, converts them to a numeric byte[] in phred scale. 
	 */
	public static void convertQualityToPhred(byte[] quals, final FastqQualityFormat version) {
		switch (version)  {
		case Standard:
			SAMUtils.fastqToPhred(quals);
			break ;
		case Solexa:
			SolexaQualityConverter.getSingleton().convertSolexaQualityCharsToPhredBinary(quals);
			break ;
		case Illumina:
			SolexaQualityConverter.getSingleton().convertSolexa_1_3_QualityCharsToPhredBinary(quals);
			break ;
		}
	}
	
	
	 /**
     * Looks at fastq input(s) and attempts to determine the proper quality format
     *
     * Closes the reader(s) by side effect
     *
     * @param readers readers on the input fastq files
     * @param expectedQuality If provided, will be used for sanity checking. If left null, autodetection will occur
     */
    public static FastqQualityFormat determineQualityFormat(final FastqReader [] readers, final FastqQualityFormat expectedQuality) {
        final QualityEncodingDetector detector = new QualityEncodingDetector();

        //add all fastq readers
        detector.add(QualityEncodingDetector.DEFAULT_MAX_RECORDS_TO_ITERATE, readers);
        //close all readers
        for (FastqReader reader : readers) {
        		reader.close();
		}
        //
        final FastqQualityFormat qualityFormat =  detector.generateBestGuess(QualityEncodingDetector.FileContext.FASTQ, expectedQuality);
        //in case there is no expected quality and different options were possible, warn user
        if (detector.isDeterminationAmbiguous()) {
            log.warn("Making ambiguous determination about fastq's quality encoding; more than one format possible based on observed qualities.");
        }
        
        	log.info(String.format("Auto-detected quality format as: %s.", qualityFormat));

        return qualityFormat;
    }
	
	
	/**
	 * @return the result of executing whoami on the underlying OS 
	 */
	public static String whoami() {
		String command = "whoami" ;
        Process myprocess = null;
        try {
            myprocess = Runtime.getRuntime().exec(command);
            ProcessConsumer pc = new ProcessConsumer(myprocess);
            StringBuffer err = new StringBuffer();
            StringBuffer out = new StringBuffer();
            
            int rc = pc.output( out ) 
                .error( err ) 
                .consume();
            
            if (rc != 0) {
                return null;
            }
            return	out.toString().replaceAll("\n", "");
            
            
        } catch (Exception e) {
			return null;
		}finally{
            if(myprocess!=null){
            	myprocess.destroy();
				myprocess = null;
            }
        }
         
		
	}

	public static int barcodeSlotCount(ReadLayout[] readLayouts) {
		
		return barcodeBlockUniqueIdSet(readLayouts).size();
	}
	
	public static Set<Integer> barcodeBlockUniqueIdSet(ReadLayout[] readLayouts) {
		Set<Integer> allIds = new TreeSet<Integer>();
		for (ReadLayout rl : readLayouts) {
			allIds.addAll(rl.getBarcodeBlockUniqueIds());
		}
		return allIds;
	}
	
	
	
	
	
	
	
}
