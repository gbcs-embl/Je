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

import java.util.List;

public interface MolecularCodedReadEnds {

	
	/**
	 * Set the molecular code of the underlying object
	 * @param code the code as a String coded using ATGCN
	 */
	public void setMolecularCode(String code);
	
	/**
	 * Get the molecular code of the underlying object
	 * @return the code as a String coded using ATGCN
	 */
	public String getMolecularCode();
	
	
	/**
	 * Get the molecular code of the underlying object as a list of codes of identical length
	 * @param codelen length of an individual code. Each code is supposed to have the samelength
	 * @return the molecularCodes
	 */
	public List<String> getMolecularCodes(int codelen);
	
}
