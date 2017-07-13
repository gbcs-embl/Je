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

import com.developpez.adiguba.shell.ProcessConsumer;

public class JeUtils {

	
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
