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

import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

/**
 * 
 * Just a class to gather at a unique place all the configuration parameters.
 * Since these props should not be visible to users, we though it is better to embed 
 * them this way within the code instead of exposing them in a config file. The drawback is 
 * of course than one has to republish the jars upon config changes. On the other hand 
 * it should not be very frequent! 
 * 
 * @author girardot
 *
 */
public class ApplicationConfiguration {

	
	public static String DB_USER = null;
	public static String DB_PWD = null;
	public static String DB_URL = null; 
	public static String CONTACT_EMAIL = null;
	
	
	
	
	public static void init() throws IOException {
		InputStream is = ApplicationConfiguration.class.getResourceAsStream("/jemultiplexer.properties"); //at the root
		Properties p = new Properties();
		p.load(is);
		DB_URL = p.getProperty("embase.url");
		DB_USER = p.getProperty("embase.user");
		DB_PWD = p.getProperty("embase.pwd");
		CONTACT_EMAIL = p.getProperty("email");
	}
	
	
	
}
