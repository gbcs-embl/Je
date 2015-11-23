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
package org.embl.gbcs.je.jeduplicates;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.lang3.StringUtils;
import org.embl.cg.utilitytools.utils.stats.impl.SimpleNumberRandomizer;
import org.embl.gbcs.je.jeduplicates.MarkDuplicatesWithMolecularCode;
import org.junit.Assert;
import org.junit.Test;

public class SplitByBarcodesTest {

	
	@Test
	public void allCodesShouldBeSuccessfullyMatchedToTheirOrginalCodes(){
		
		MarkDuplicatesWithMolecularCode m = new MarkDuplicatesWithMolecularCode();
		
		/*
		 * get some barcoded objects with errors
		 * only ask up to 2 MMs so that all errors should fall back to its original barcode 
		 */
		int barcodelength = 8;
		String[] originalbarcodes = new String[]{
				StringUtils.repeat("A", barcodelength),
				StringUtils.repeat("T", barcodelength),
				StringUtils.repeat("G", barcodelength),
				StringUtils.repeat("C", barcodelength),
		};
		Map<String, List<FakeReadEndWithMolecularCode>> map = generateForBarcodes(originalbarcodes, 2 , 8);
		/*
		 * l contains 3 times (0,1,2 MMs) 10 objects ie 30 for each barcode
		 */
		int generatedRandomCodesPerOriCode = 30;
		
		
		Assert.assertTrue("barcode generation code is buggy", map.size() == originalbarcodes.length);
		List<FakeReadEndWithMolecularCode> l = new ArrayList<FakeReadEndWithMolecularCode>();
		System.out.println("Created lists: ");
		for (Entry<String, List<FakeReadEndWithMolecularCode>> e : map.entrySet()) {
			Assert.assertTrue("barcode generation code is buggy", e.getValue().size() == generatedRandomCodesPerOriCode);
			printListContent(e.getKey(), e.getValue());
			l.addAll(e.getValue());
		}
		
		
		/*
		 * TEST NOW WITH NO PRESET CODE LIST
		 * 
		 */
		int mm = 2;
		int maxN = 2;
		Map<String, List<FakeReadEndWithMolecularCode>> ll = m.splitDuplicatesByMolecularCodeGroupWithoutPredefinedCodeList(
				l, 
				mm, // 2 MM allowed
				maxN); // max 2 Ns

		//note that we expect 5 lists due to the UNDEF  
		int expectedUNDEFEntriesNumPerOriCode = (mm+1)  * (3 - maxN); // 3 is the max number of N that the method generateForBarcodes() includes in the similated data
		int expectedUNDEFEntriesNum = expectedUNDEFEntriesNumPerOriCode * originalbarcodes.length; // 3 is the max number of N that the method generateForBarcodes() includes in the similated data
		Assert.assertTrue("Did not get right number of lists back", ll.size() == (originalbarcodes.length+1));
		
		//check content
		System.out.println("After code splitting: ");
		for (Entry<String, List<FakeReadEndWithMolecularCode>> e : ll.entrySet()) {
			if(e.getKey().equals("UNDEF")){
				printListContent(e.getKey(), e.getValue());
				Assert.assertEquals("Did not get expected object number back in list for "+e.getKey(), e.getValue().size() , expectedUNDEFEntriesNum );
			}else{

				for (FakeReadEndWithMolecularCode read : e.getValue()) {
					Assert.assertEquals("Object with code "+read.getMolecularCode()+" wrongly assigned to "+e.getKey()+" instead of "+read.originalcode , read.originalcode, e.getKey());
				}
				Assert.assertEquals("Did not get expected object number back in list for "+e.getKey(), e.getValue().size() , (generatedRandomCodesPerOriCode-expectedUNDEFEntriesNumPerOriCode) );
			}
		}
		
		
		/*
		 * TEST WITH PRESET CODE LIST ; same expectations as above
		 * 
		 */
		ll = m.splitDuplicatesByMolecularCodeGroupWithPredefinedCodeList(
				l, 
				mm, // 2 MM allowed
				maxN,// max 2 Ns
				barcodelength, 
				new TreeSet<String>(Arrays.asList(originalbarcodes))
				); 

		Assert.assertTrue("Did not get right number of lists back", ll.size() == (originalbarcodes.length+1));
		
		//check content
		System.out.println("After code splitting: ");
		for (Entry<String, List<FakeReadEndWithMolecularCode>> e : ll.entrySet()) {
			if(e.getKey().equals("UNDEF")){
				Assert.assertEquals("Did not get expected object number back in list for "+e.getKey(), e.getValue().size() , expectedUNDEFEntriesNum );
			}else{

				for (FakeReadEndWithMolecularCode read : e.getValue()) {
					Assert.assertEquals("Object with code "+read.getMolecularCode()+" wrongly assigned to "+e.getKey()+" instead of "+read.originalcode , read.originalcode, e.getKey());
				}
				Assert.assertEquals("Did not get expected object number back in list for "+e.getKey(), e.getValue().size() , (generatedRandomCodesPerOriCode-expectedUNDEFEntriesNumPerOriCode) );
			}
		}
		
		
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	private void printListContent(String key,
			List<FakeReadEndWithMolecularCode> value) {
		System.out.println("list key is "+key+" content : ");
		for (FakeReadEndWithMolecularCode c : value) {
			System.out.print(c.getMolecularCode() + " ");
		}
		System.out.println();
		
	}


	public static Map<String, List<FakeReadEndWithMolecularCode>> generateForBarcodes(String [] bcs, int upToMMNumber, int codelen){
		
		Map<String, List<FakeReadEndWithMolecularCode>> m = new LinkedHashMap<String, List<FakeReadEndWithMolecularCode>>();
		
		int i = 0;
		for(String bc : bcs){
			List<FakeReadEndWithMolecularCode> l = new ArrayList<FakeReadEndWithMolecularCode>();
			for(int mm = 0 ; mm <= upToMMNumber; mm++){
				// add 10 objects with the exact code but mm mismathces; with some Ns that should not be considered during comparison
				l.add( new FakeReadEndWithMolecularCode(i++, bc,injectMMs(bc, mm)) );
				l.add( new FakeReadEndWithMolecularCode(i++, bc, injectMMs(bc, mm)) );
				l.add( new FakeReadEndWithMolecularCode(i++, bc,injectMMs(bc, mm)) );
				l.add( new FakeReadEndWithMolecularCode(i++, bc, injectMMs(bc, mm)) );
				l.add( new FakeReadEndWithMolecularCode(i++, bc,injectMMs(bc, mm)) );
				l.add( new FakeReadEndWithMolecularCode(i++, bc, injectMMs(bc, mm)) );
				l.add( new FakeReadEndWithMolecularCode(i++, bc,injectMMs(bc, mm)) );
				l.add( new FakeReadEndWithMolecularCode(i++, bc, injectNsAndMMs(bc, 1, mm)) );
				l.add( new FakeReadEndWithMolecularCode(i++, bc, injectNsAndMMs(bc, 2, mm)) );
				l.add( new FakeReadEndWithMolecularCode(i++, bc, injectNsAndMMs(bc, 3, mm)) );
			}
			m.put(bc, l);
		}
	
		return m;
		
	}

	
	
	private static String injectNsAndMMs(String bc, int numN, int numMMs) {
		SimpleNumberRandomizer rdm = new SimpleNumberRandomizer();
		StringBuilder newcode = new StringBuilder(bc);
		int injectedNs = 0;
		for(int pos : rdm.getRandomIntegers(false, bc.length(), numN+numMMs) ){
			if(injectedNs<numN){
				newcode.setCharAt(pos, 'N');
				injectedNs++;
			}else{
				//get random replacement
				List<Character> fourbases = new ArrayList<Character>( Arrays.asList( new Character[]{'A', 'C', 'G', 'T' }) );
				Character b = newcode.charAt(pos);
				fourbases.remove(b);
				Collections.shuffle(fourbases);
				newcode.setCharAt(pos, fourbases.get(0));
			}
		}
		return newcode.toString();
	}
	
	
	
	private static String injectMMs(String bc, int num) {
		if(num ==0 ) return bc;
		SimpleNumberRandomizer rdm = new SimpleNumberRandomizer();
		StringBuilder newcode = new StringBuilder(bc);
		for(int pos : rdm.getRandomIntegers(false, bc.length(), num) ){
			//get random replacemetn
			List<Character> fourbases = new ArrayList<Character>( Arrays.asList( new Character[]{'A', 'C', 'G', 'T' }) );
			Character b = newcode.charAt(pos);
			fourbases.remove(b);
			Collections.shuffle(fourbases);
			newcode.setCharAt(pos, fourbases.get(0));
		}
		return newcode.toString();
	}
	
	
	private static String injectNs(String bc, int num) {
		SimpleNumberRandomizer rdm = new SimpleNumberRandomizer();
		StringBuilder newcode = new StringBuilder(bc);
		for(int pos : rdm.getRandomIntegers(false, bc.length(), num) ){
			newcode.setCharAt(pos, 'N');
		}
		return newcode.toString();
	}
	
	
}















