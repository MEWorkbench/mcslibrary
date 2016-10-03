package pt.uminho.ceb.biosystems.mcslibrary.utilities;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.category.CategoryDataset;
import org.jfree.ui.ApplicationFrame;

public class BarChart extends ApplicationFrame{
		private static final long serialVersionUID = -1533068624614638098L;

		public BarChart(String chartTitle, String applicationTitle, CategoryDataset dataset) {
		      super( applicationTitle );
		      JFreeChart barChart = ChartFactory.createBarChart(
		         chartTitle,
		         "Category",
		         "Score",
		         dataset,
		         PlotOrientation.HORIZONTAL,
		         true, true, false);

		      ChartPanel chartPanel = new ChartPanel( barChart );
		      chartPanel.setPreferredSize(new java.awt.Dimension( 560 , 367 ) );
		      setContentPane( chartPanel );
		}
}
