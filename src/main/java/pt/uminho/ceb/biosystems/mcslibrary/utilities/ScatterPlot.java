package pt.uminho.ceb.biosystems.mcslibrary.utilities;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYDataset;
import org.jfree.ui.ApplicationFrame;

public class ScatterPlot extends ApplicationFrame{
	/**
	 *
	 */
	private static final long serialVersionUID = 2851566693067145867L;

	public ScatterPlot(String title, String applicationTitle, String xlab, String ylab, XYDataset dataset) {
	      super( applicationTitle );
	      JFreeChart barChart = ChartFactory.createScatterPlot(
	    		  title,
	    		  xlab,
	    		  ylab,
	    		  dataset, PlotOrientation.VERTICAL, rootPaneCheckingEnabled, rootPaneCheckingEnabled, rootPaneCheckingEnabled);

	      ChartPanel chartPanel = new ChartPanel( barChart );
	      chartPanel.setPreferredSize(new java.awt.Dimension( 560 , 367 ) );
	      setContentPane( chartPanel );
	      
	}
}
