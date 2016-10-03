package pt.uminho.ceb.biosystems.mcslibrary.utilities;

import org.jfree.data.category.CategoryDataset;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RefineryUtilities;

public class ChartUtilities {

	public static CategoryDataset createCategoryDataset(String[] labels, String[] series, double[][] values){
		DefaultCategoryDataset ds = new DefaultCategoryDataset();
		for (int i = 0; i < values.length; i++) {
			for (int j = 0; j < values[0].length; j++) {
				ds.addValue(values[i][j], series[i], labels[j]);
			}
		}
		return ds;
	}

	public static XYDataset createXYDataset(double[][] xValues, double[][] yValues, String[] seriesName) {
		XYSeriesCollection ds = new XYSeriesCollection();
		for (int i = 0; i < seriesName.length; i++) {
			XYSeries series = new XYSeries(seriesName[i]);
			for (int j = 0; j < xValues[i].length; j++) {
				series.add(xValues[i][j], yValues[i][j]);
			}
			ds.addSeries(series);
		}
		return ds;

	}
	public static void displayPlot(BarChart chart){
	      chart.pack();
	      RefineryUtilities.centerFrameOnScreen(chart);
	      chart.setVisible(true);
	   }

	public static void displayPlot(ScatterPlot chart){
	      chart.pack();
	      RefineryUtilities.centerFrameOnScreen(chart);
	      chart.setVisible(true);
	   }


}
