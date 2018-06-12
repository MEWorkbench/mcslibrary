package pt.uminho.ceb.biosystems.mcslibrary.datastructures.dataframe.components;

import java.util.List;

public class Series {
	private List<String> series;
	private String label;

	public Series(List<String> series, String name){
		this.series = series;
		this.label = name;
	}
	
	public int size(){
		return this.series.size();
	}
	
	public String getElement(int i){
		return this.series.get(i);
	}
	
	public String getLabel(){
		return this.label;
	}
	
	public List<String> getSeries(){
		return series;
		
	}
}
