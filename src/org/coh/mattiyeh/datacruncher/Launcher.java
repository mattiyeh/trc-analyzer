package org.coh.mattiyeh.datacruncher;

import java.io.IOException;
import java.net.InetAddress;

public class Launcher {
	
	private static String hostname = "";

	public static void main(String[] args) throws NumberFormatException, IOException {
	
		InetAddress addr;
	    addr = InetAddress.getLocalHost();
	    hostname = addr.getHostName();
		
		DataCruncher dc = new DataCruncher();
		dc.go();
	
	}
	
	public static String getHostname() {
		return hostname;
	}

}
