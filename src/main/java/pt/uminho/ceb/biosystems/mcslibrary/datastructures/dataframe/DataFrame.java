package pt.uminho.ceb.biosystems.mcslibrary.datastructures.dataframe;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;

import pt.uminho.ceb.biosystems.mcslibrary.datastructures.dataframe.components.Series;

public class DataFrame {
	private int columns;
	private int rows;
	
	LinkedHashMap<String, Series> dataMap;
	public DataFrame(List<List<String>> data, boolean header) {
		int start = 0;
		if (header) {
			start++;
		}
		this.columns = data.get(0).size();
		dataMap = new LinkedHashMap<>();
		for (int i = 0; i < columns; i++) {
			List<String> s = new ArrayList<>();
			for (int j = start; j < data.size(); j++) {
				s.add(data.get(j).get(i));
			}
//			System.out.println(s);
			if (header) {
				dataMap.put(data.get(0).get(i), new Series(s, data.get(0).get(i)));
			} else {
				dataMap.put(String.valueOf(i), new Series(s, String.valueOf(i)));
			}
		}
		this.rows = data.size() - start;
		
	}
	
	public Series getColumn(String name){
		return dataMap.get(name);
	}
		
	public List<String> getColumnNames(){
		return new ArrayList<String>(dataMap.keySet());
	}

	/**
	 * @return the rows
	 */
	public int numberOfRows() {
		return rows;
	}

	public String getItem(int rows, int columns){
		return dataMap.get(getColumnNames().get(columns)).getElement(rows);
	}
	
	/**
	 * @return the columns
	 */
	public int numberOfColumns() {
		return columns;
	}

}
