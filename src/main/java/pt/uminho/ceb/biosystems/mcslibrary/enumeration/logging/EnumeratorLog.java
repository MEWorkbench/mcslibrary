package pt.uminho.ceb.biosystems.mcslibrary.enumeration.logging;

import java.util.HashMap;

import pt.uminho.ceb.biosystems.mcslibrary.enumeration.implementation.CPLEXIntegratedEnumerator;

public class EnumeratorLog {
	int currentPrintEntry = 0;
	private long starttime;
	private HashMap<Integer,LogEntry> entries;
	
	public EnumeratorLog(CPLEXIntegratedEnumerator en) {
		starttime = System.currentTimeMillis();
		entries = new HashMap<Integer, LogEntry>();
	}
	
	public void addEntry(String entryType, String entryText) {
		long curTime = System.currentTimeMillis() - this.starttime;
		entries.put(entries.size(), new LogEntry(entryType, entryText, this, curTime));
		printNext();
	}

	private void printNext() {
		System.out.println(entries.get(currentPrintEntry));
		currentPrintEntry ++;
	}
	
}
