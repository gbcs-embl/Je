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

import htsjdk.samtools.Defaults;
import htsjdk.samtools.fastq.AsyncFastqWriter;
import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.util.Md5CalculatingOutputStream;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.zip.GZIPOutputStream;

public class JemultiplexerFastqWriterFactory {

	boolean useAsyncIo = Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS;
	

	/** Sets whether or not to use async io (i.e. a dedicated thread per writer. */
    public void setUseAsyncIo(final boolean useAsyncIo) { this.useAsyncIo = useAsyncIo; }

	
	public FastqWriter newWriter(final File out, boolean gzip, boolean createMd5) {
		
		if(!gzip)
			return newWriter(out, createMd5);
					
					
		GZIPOutputStream go;
		try {
			if(createMd5)
				go = new GZIPOutputStream(
						new Md5CalculatingOutputStream(
								new FileOutputStream(out),
								new File(out.getAbsolutePath() + ".md5"))
						);
			else
				go = new GZIPOutputStream(new FileOutputStream(out));
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
        final FastqWriter writer = new BasicFastqWriter( new PrintStream(go) );
        if (useAsyncIo) {
            return new AsyncFastqWriter(writer, AsyncFastqWriter.DEFAULT_QUEUE_SIZE);
        }
        else {
            return writer;
        }
    }
	
	
	public FastqWriter newWriter(final File out, boolean createMd5) {
		if(!createMd5)
			return newWriter(out);
		
		
		Md5CalculatingOutputStream go;
		try {
			go =new Md5CalculatingOutputStream(
					new FileOutputStream(out),
					new File(out.getAbsolutePath() + ".md5"));

		} catch (IOException e) {
			throw new RuntimeException(e);
		}
        final FastqWriter writer = new BasicFastqWriter( new PrintStream(go) );
	    if (useAsyncIo) {
	        return new AsyncFastqWriter(writer, AsyncFastqWriter.DEFAULT_QUEUE_SIZE);
	    }
	    else {
	        return writer;
	    }
	}
	
	public FastqWriter newWriter(final File out) {
	    final FastqWriter writer = new BasicFastqWriter(out);
	    if (useAsyncIo) {
	        return new AsyncFastqWriter(writer, AsyncFastqWriter.DEFAULT_QUEUE_SIZE);
	    }
	    else {
	        return writer;
	    }
	}
	
}
