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

import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.FastqQualityFormat;
import htsjdk.samtools.util.QualityEncodingDetector;
import htsjdk.samtools.util.SolexaQualityConverter;

import java.io.File;
import java.util.Arrays;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class JeTry {
	private static Logger log = LoggerFactory.getLogger(JeTry.class);

	public JeTry() {
		// TODO Auto-generated constructor stub
	}

	public static void main(String[] args) {
		String f = "/g/furlong/incoming/2017-06-28-000000000-B954B/000000000-B954B_precap_allpromV2_17s003049-1-1_Ghavi-helm_lane117s003049_1_sequence.txt.gz";
		FastqReader reader = new FastqReader(new File(f));
		FastqQualityFormat QUALITY_FORMAT = QualityEncodingDetector.detect(100000, reader);
		log.info(String.format("Auto-detected quality encoding format as: %s.", QUALITY_FORMAT));
		
		FastqRecord r  =reader.iterator().next();
		System.out.println(r.toFastQString());
		System.out.println(r.getBaseQualityString());

		byte [] bites = Arrays.copyOf( r.getBaseQualityString().getBytes() , r.getBaseQualityString().getBytes().length);
		System.out.println(Arrays.toString( bites ));
		SAMUtils.fastqToPhred(bites);
		System.out.println(Arrays.toString( bites ));
		
		bites = Arrays.copyOf( r.getBaseQualityString().getBytes() , r.getBaseQualityString().getBytes().length);
		System.out.println(Arrays.toString( bites ));
		SolexaQualityConverter.getSingleton().convertSolexa_1_3_QualityCharsToPhredBinary(bites);
		System.out.println(Arrays.toString( bites ));
		
		
	}

}
