/*
 * To the extent possible under law, the ImageJ developers have waived
 * all copyright and related or neighboring rights to this tutorial code.
 *
 * See the CC0 1.0 Universal license for details:
 *     http://creativecommons.org/publicdomain/zero/1.0/
 */

package bio.coil.CoilEdinburgh;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.ImageWindow;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.gui.ProfilePlot;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.measure.ResultsTable;
import ij.plugin.Duplicator;
import ij.plugin.ZProjector;
import ij.plugin.filter.Analyzer;
import ij.plugin.frame.RoiManager;
import ij.process.FloatPolygon;
import ij.process.ImageStatistics;
import io.scif.services.DatasetIOService;
import io.scif.services.FormatService;

import net.imagej.Dataset;
import net.imagej.ImageJ;
import net.imagej.ops.OpService;
import net.imagej.roi.ROIService;
import net.imglib2.Point;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

import org.scijava.command.Command;
import org.scijava.convert.ConvertService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;
import java.io.IOException;
import java.io.File;

import java.util.*;
import java.util.List;

import javax.swing.JOptionPane;

/**
 *
 * 
 */
@Plugin(type = Command.class, menuPath = "Plugins>Users Plugins>Fungal Tip Track")
public class Fungal_Track<T extends RealType<T>> implements Command {
	@Parameter
	File Imagefile;
	
	@Parameter
	File labkitFile;
	
    @Parameter
    private FormatService formatService;

    @Parameter
    private DatasetIOService datasetIOService;

    @Parameter
    private UIService uiService;
    
    @Parameter
    private ConvertService cs;

    @Parameter
    private OpService ops;
    
    @Parameter
    static ImageJ ij;

    @Parameter
    private ROIService roiService;
 
    
    @Override
    public void run() {
    	Dataset currentData = null;	
    	try {
			currentData = ij.scifio().datasetIO().open(Imagefile.getPath());
		} catch (IOException e) {
			e.printStackTrace();
		}

		@SuppressWarnings("unchecked")
		Img<T> image = (Img<T>) currentData.getImgPlus();  
		
		List<IntervalView<T>>convertedImage = SplitChannels(image);
		
		IdentifyTip(convertedImage); 
        IJ.run("Close All", "");
        new WaitForUserDialog("Finished", "Plugin Finished").show();
    }

    public List<IntervalView<T>> SplitChannels(Img<T> image){
    	List<IntervalView<T>> theChannels = new ArrayList<IntervalView<T>>();
    	
    	// defining which dimension of the image represents channels
   	   	int channelDim = 2;
    
    	// how many channels do we have?
   	   	long numChannels = image.dimension(channelDim);
   	   	for (int channelIndex = 0; channelIndex < numChannels; channelIndex++) {
   	   		IntervalView<T> inputChannel = ij.op().transform().hyperSliceView(image, channelDim, channelIndex); 
   	   		theChannels.add(inputChannel);
   	   	}
    	
    	return theChannels;
    }
    
    public void IdentifyTip(List<IntervalView<T>> convertedImage) {
    	
    	//Select the brightfield area
    	IntervalView<T> brightfieldImage = convertedImage.get(0);
    	ImagePlus impProjected = PrepImages(brightfieldImage);
    	impProjected.show();
    	int impProjectedID = impProjected.getID();
    	//Select the same fluorescence area
    	IntervalView<T> fluorescenceImage = convertedImage.get(1);
    	ImagePlus impProjectedFluorescence = PrepImages(fluorescenceImage);
    	impProjectedFluorescence.show();
    	int impProjectedFluorescenceID = impProjectedFluorescence.getID();
    	
    	//Select all the interesting fungal tips
    	List<ImagePlus> tipSelection = new ArrayList<>();
    	List<ImagePlus> tipFluorescentSelection = new ArrayList<>();
    	RoiManager rm = RoiManager.getRoiManager();
    	IJ.setTool("rectangle");
    	String theAnswer = "y";
    	do {
    		IJ.selectWindow(impProjectedID);
    		new WaitForUserDialog("Draw ROI", "Draw an ROI around a fungal tip").show();
    		rm.addRoi(null);
    		theAnswer = JOptionPane.showInputDialog("Do you want to ROI another tip, y/n ");
    	}while(theAnswer.equals("y"));
    	//Duplicate the fungal tips
    	for (int c=0;c<rm.getCount();c++) {
    		IJ.selectWindow(impProjectedID);
    		rm.select(c);
    		ImagePlus tempTip = new Duplicator().run(impProjected, 1, 3);
    		IJ.selectWindow(impProjectedFluorescenceID);
    		rm.select(c);
    		ImagePlus tempFluoroTip = new Duplicator().run(impProjectedFluorescence, 1, 3);
    		
    		tipSelection.add(tempTip.duplicate());
    		tipFluorescentSelection.add(tempFluoroTip.duplicate());
    		
    	}
    	
    	//Convert selected tips to masks using Labkit
    
    	ImagePlus [] impCollection = new ImagePlus [impProjected.getNChannels()];
    	//Get labkit file location to create string command
    	String name = labkitFile.getAbsolutePath();	
    	for (int a=0;a<tipSelection.size();a++) {
    		ImagePlus maskTemp = tipSelection.get(a);
    		maskTemp.show();
    		String impName = maskTemp.getTitle();	
    		String inputString = "input=" + impName +" segmenter_file=" +name+" use_gpu=true";	
    		IJ.run("Segment Image With Labkit", inputString);
    		ImagePlus maskedImage = WindowManager.getCurrentImage();
    		IJ.setAutoThreshold(maskedImage, "Default dark");
    		IJ.run(maskedImage, "Analyze Particles...", "size=2000-Infinity show=Masks include stack");
    		ImagePlus temp = WindowManager.getCurrentImage();
    		impCollection[a] = temp.duplicate();
    		maskTemp.changes=false;
    		maskTemp.close();
    		maskedImage.changes=false;
    		maskedImage.close();	
    		temp.changes=false;
    		temp.close();
    	}      
    	
    	clearRoiManager(impProjectedID);

    	for (int i=0;i<tipSelection.size();i++) {
    		ImagePlus hyphalTipImage = impCollection[i];
    		hyphalTipImage.show();
    		for (int s=1;s<hyphalTipImage.getNSlices()+1;s++) {
    			hyphalTipImage.setSlice(s);    			
    			IJ.setAutoThreshold(hyphalTipImage, "Default");
    			IJ.run(hyphalTipImage, "Analyze Particles...", "size=3000-Infinity show=Masks include add slice");
    			ImagePlus temp = WindowManager.getCurrentImage();
    			int tempID = temp.getID();
    			IJ.run(temp, "Gaussian Blur...", "sigma=3 slice");
    			Prefs.blackBackground = false;
    			IJ.run(temp, "Convert to Mask", "");
     	        IJ.run(temp, "Skeletonize", "");
    	        clearRoiManager(tempID);
    	        IJ.setAutoThreshold(temp, "Default");
    	        IJ.run(temp, "Analyze Particles...", "  show=Nothing add");
    	        Roi lineRoi = rm.getRoi(0);
    	        //Cycle through roi manager to pick up broken lines
    	        int longLine = 0;
    	        for (int z=0;z<rm.getCount();z++) {
    	        	rm.select(z);
    	        	FloatPolygon numPoints = lineRoi.getFloatPolygon();
    	        	longLine = longLine + numPoints.npoints;
    	        }
    	      
    	        int lineLength = longLine;
    	        
    	        //Find tip direction for fluorescence measurements
    	        int[] dims = temp.getDimensions();
    	        Roi lineProfile = rm.getRoi(0);
    	        float [] xPos = lineProfile.getFloatPolygon().xpoints;
    	        float [] yPos = lineProfile.getFloatPolygon().ypoints;
    	        int [] tipCoords = FindTip(xPos, yPos,dims);
    	        
    	        //Extract same frame and area in fluorescence channel and measure it
    	        ImagePlus fluorescentTemp = tipFluorescentSelection.get(i);
    	        
    	        MeasureFluorescence(temp, tempID, fluorescentTemp,rm);
    	        clearRoiManager(tempID);
    	        
    	        temp.changes=false;
    	        temp.close();
    		}
    	}
    	
    }
    
    private ImagePlus PrepImages(IntervalView<T> imageToProject) {
    	    	
    	imageToProject = (IntervalView<T>) Views.moveAxis(imageToProject, 2,3);
    	ImagePlus impTips = ImageJFunctions.wrap(imageToProject,"Working_Image");
    	ImagePlus impProjected = ZProjector.run(impTips,"Max_Working_Image");
    	
    	return impProjected;
    }
  
    private void MeasureFluorescence(ImagePlus temp, int tempID, ImagePlus fluorescentTemp, RoiManager rm) {
    	
    	IJ.selectWindow(tempID);
    	IJ.run(temp, "Select None", "");
    	Prefs.blackBackground = false;
    	for(int d=0;d<4;d++) {
    		IJ.run(temp, "Dilate", "");
    	}
        IJ.setAutoThreshold(temp, "Default");
		IJ.run(temp, "Analyze Particles...", "  show=Nothing add"); 
		fluorescentTemp.show();
		int fluorescentTempID = fluorescentTemp.getID();
		IJ.selectWindow(fluorescentTempID);
		rm.select(1);
        IJ.run(fluorescentTemp, "Area to Line", "");
    	IJ.run(fluorescentTemp, "Plot Profile", "");
    	
    	Roi plotValues = fluorescentTemp.getRoi();
    	java.awt.Point[] linePoints = plotValues.getContainedPoints();
    	
    	 // Get the active plot
    	PlotWindow pw;

		ImagePlus imp = WindowManager.getCurrentImage();   //Get Plot Window
		if (imp==null){ 
			IJ.error("There are no plots open."); return;   //Check for open line profile plots
		}    
		ImageWindow win = imp.getWindow();
		
			pw = (PlotWindow)win;
			float[] yvalues = pw.getYValues();	 //Get Y values from plot window
			float[] xvalues = pw.getXValues();   //Get X values from plot window
		
       
    	int stop =0;
    }
 
   
    public int [] FindTip(float[] xPos,float[] yPos, int[] dims) {

    	int [] tipCoords = new int [2];
    	
    	//is startpoint right
    	float xPoint = 0;
    	float yPoint = 0;
    	RoiManager rm = RoiManager.getRoiManager();
    	if(xPos[0]>(dims[0]/2)){
    		float tempIntX = dims[0];
    		float tempIntY = dims[1];
    		float [] coordTip = RightTip(xPos, yPos,tempIntX,tempIntY);
    		
    
    		Roi[] test = rm.getRoisAsArray();
        	FloatPolygon coords = test[0].getInterpolatedPolygon();
       
        	float [] lineXPos = coords.xpoints;
	    	float [] lineYPos = coords.ypoints;
	    	float [] coordLineTip = RightTip(lineXPos, lineYPos,tempIntX,tempIntY);
	    	tipCoords[0]=(int) coordLineTip[0];
	    	tipCoords[1]=(int) coordLineTip[1];
     
    	}
    	
    	//is startingpoint left
    	xPoint = 0;
    	yPoint = 0;
    	if(xPos[0]<(dims[0]/2)) {
    		float tempIntX = 0;
    		float tempIntY = 0;
    		float [] coordTip = LeftTip(xPos, yPos,tempIntX,tempIntY);
    
    		Roi[] test = rm.getRoisAsArray();
        	FloatPolygon coords = test[0].getInterpolatedPolygon();

        	float [] lineXPos = coords.xpoints;
	    	float [] lineYPos = coords.ypoints;
	    	float [] coordLineTip = LeftTip(lineXPos, lineYPos,tempIntX,tempIntY);
	    	tipCoords[0]=(int) coordLineTip[0];
	    	tipCoords[1]=(int) coordLineTip[1];
    
    	}
    
		return tipCoords;
    	
    }
    
    public float[] RightTip(float[] xPos,float[] yPos, float tempIntX, float tempIntY) {
    	
    	float [] coordTip = new float[2];
    	for (int a=0;a<xPos.length;a++) {
			if (xPos[a]<tempIntX) {
				tempIntX=xPos[a];
				tempIntY=yPos[a];
			}
		}
    	coordTip[0] = tempIntX;
    	coordTip[1] = tempIntY;
    	return coordTip;
		
    }
    
    public float[] LeftTip(float[] xPos,float[] yPos, float tempIntX, float tempIntY) {
    	float [] coordTip = new float[2];
    	for (int a=0;a<xPos.length;a++) {
			if (xPos[a]>tempIntX) {
				tempIntX=xPos[a];
				tempIntY=yPos[a];
			}
		}
    	coordTip[0] = tempIntX;
    	coordTip[1] = tempIntY;
    	return coordTip;
    }
   
    
    public void ClearResults(){
		
		ResultsTable emptyrt = new ResultsTable();	
		emptyrt = Analyzer.getResultsTable();
		int valnums = emptyrt.getCounter();
		for(int a=0;a<valnums;a++){
			IJ.deleteRows(0, a);
		}
	}
    
    /*
    * Function to remove all regions from 
    * the region of interest manager
    */
    public void clearRoiManager(int impID) {
    	IJ.selectWindow(impID);
        ImagePlus imp = IJ.getImage();
        RoiManager rm = RoiManager.getRoiManager();
      
        if(rm.getCount()>0) {
        	rm.runCommand(imp,"Deselect");
        	rm.runCommand(imp,"Delete");
        }
    }


    public static void main(final String... args) throws Exception {
        // create the ImageJ application context with all available services
        final ImageJ ij = new ImageJ();
        ij.ui().showUI();
        ij.command().run(Fungal_Track.class, true);
    }

}

