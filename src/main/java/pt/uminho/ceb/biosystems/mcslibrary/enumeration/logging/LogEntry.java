package pt.uminho.ceb.biosystems.mcslibrary.enumeration.logging;

public class LogEntry {
	private String entryType;
	private String entryText;
	private long time;
	
	public LogEntry(String entryType, String entryText, EnumeratorLog l, long time) {
		this.entryText = entryText;
		this.entryType = entryType;
		this.time = time;
	}
	
	public String getText() {
		return this.entryText;
	}
	
	public String getType() {
		return this.entryType;
	}

	public long getTime() {
		return time;
	}
	
	public String toString() {
		return "["+((double)time/1000)+"s]"+"["+entryType+"]: "+entryText;
	}

}
