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
import ij.WindowManager;
import ij.gui.Roi;
import ij.plugin.ChannelSplitter;
import ij.plugin.ZProjector;
import ij.plugin.frame.RoiManager;
import ij.process.ImageProcessor;
import io.scif.bf.BioFormatsFormat;
import io.scif.services.DatasetIOService;
import io.scif.services.FormatService;
import net.imagej.ImageJ;
import net.imagej.ops.OpService;
import net.imagej.roi.ROIService;
import net.imglib2.type.numeric.RealType;
import org.apache.commons.io.FilenameUtils;
import org.scijava.command.Command;
import org.scijava.convert.ConvertService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.Date;

/**
 *
 * 
 */
@Plugin(type = Command.class, menuPath = "Plugins>Users Plugins>Sara Plugin")
public class Sara_Plugin<T extends RealType<T>> implements Command {

	@Parameter
	private FormatService formatService;

	@Parameter
	private DatasetIOService datasetIOService;

	@Parameter
	private UIService uiService;

	@Parameter
	private OpService ops;

	@Parameter
	private ROIService roiService;

	@Parameter(label = "Open Folder: ", style="directory")
	public File filePath;

	@Parameter(label = "Environment Path: ", style = "directory")
	public File envpath;

	@Parameter(label = "Cellpose Model Path: ")
	public File modelpath;

	RoiManager roiManager;
	String filename;
	double pixelsize;


	@Override
    public void run() {

		if(RoiManager.getInstance()!=null){
			roiManager= RoiManager.getInstance();}else{
			roiManager = new RoiManager();
		}
		roiManager.reset();

		//Open folder
		File[] files = filePath.listFiles();

		//For each .czi file
		for (File file : files) {
			if (file.toString().contains(".czi") && !file.toString().contains("Overlay")) {

				//Open file and get filename and filepath
				IJ.run("Bio-Formats Importer", "open=["+file.getAbsolutePath()+"] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
				ImagePlus imp = WindowManager.getCurrentImage();
				pixelsize = imp.getCalibration().pixelHeight;
				filename = FilenameUtils.removeExtension(file.getName());

				//Z-Project
				ImagePlus impZ = ZProjector.run(imp, "Average");
				ImagePlus[] split = ChannelSplitter.split(impZ);
				split[0].show();
				split[1].show();
				split[2].show();

				//Run Cellpose on channel 3
				int size = (int) ( 14/pixelsize);
				Cellpose_Wrapper cpw = new Cellpose_Wrapper(modelpath.getPath(), envpath.getPath(), size, split[2]);
				cpw.run(true);
				Roi[] cellRois = roiManager.getRoisAsArray();
				roiManager.reset();

				//Threshold spots on channel 3 (Yen)
				Cell[] cells = getSpots(split[2], cellRois);

				// Measure spots and background area and intensity per cell in channels 1 and 3
				double[][] output = new double[cells.length][10];
				String[] titles = new String[10];
				for (int i = 0; i< cells.length; i++){
					titles[0] = "Cell Area (um^2)";
					output[i][0] = cells[i].totalArea()*pixelsize*pixelsize;
					titles[1] = "Cell Intensity Ch0";
					output[i][1] = cells[i].totalIntensity(split[0]);
					titles[2] = "Cell Intensity Ch1";
					output[i][2] = cells[i].totalIntensity(split[1]);
					titles[3] = "Cell Intensity Ch2";
					output[i][3] = cells[i].totalIntensity(split[2]);

					titles[4] = "Background Area (um^2)";
					output[i][4] = cells[i].cellArea()*pixelsize*pixelsize;
					titles[5] = "Background Intensity Ch0";
					output[i][5] = cells[i].cellIntensity(split[0]);
					titles[6] = "Background Intensity Ch2";
					output[i][6] = cells[i].cellIntensity(split[2]);

					titles[7] = "Spots Area (um^2)";
					output[i][7] = cells[i].spotsArea()*pixelsize*pixelsize;
					titles[8] = "Spots Intensity Ch0";
					output[i][8] = cells[i].spotsIntensity(split[0]);
					titles[9] = "Spots Intensity Ch2";
					output[i][9] = cells[i].spotsIntensity(split[2]);
				}

				//Draw output image
				ImagePlus outputImage = outputImage(impZ, cells);
				outputImage.show();
				IJ.save(WindowManager.getCurrentImage(), Paths.get(filePath.toString(), filename) + "_Overlay.tif");

				//Make results file
				try {
					MakeResults(titles, output);
				} catch (IOException e) {
					e.printStackTrace();
				}

				IJ.run("Close All", "");

			}
		}
	}

	private Cell[] getSpots(ImagePlus imp, Roi[] rois){
		Cell[] cells = new Cell[rois.length];
		imp.show();
		for (int i= 0; i< rois.length; i++) {
			imp.setRoi(rois[i]);
			IJ.run(imp, "Subtract Background...", "rolling=25");
			IJ.setAutoThreshold(imp, "Yen dark");
			IJ.run(imp, "Analyze Particles...", "size=10-Infinity pixel exclude add");
			cells[i] = new Cell(rois[i]);
			Roi[] spotROIs = roiManager.getRoisAsArray();
			roiManager.reset();
			for (int j=0 ; j< spotROIs.length; j++){
				if(cells[i].containsNucleii(spotROIs[j])){
					cells[i].addSpot(spotROIs[j]);
				}
			}
		}
		return cells;
	}

	private ImagePlus outputImage(ImagePlus imp, Cell[] cells){
		ImagePlus[] split = ChannelSplitter.split(imp.duplicate());
		split[0].show();
		split[1].show();
		split[2].show();
		ImagePlus overlay = IJ.createImage("Overlay", "16-bit", split[1].getWidth(), split[1].getHeight(),
				split[1].getNChannels(), split[1].getNSlices(), split[1].getNFrames());
		overlay.show();
		drawNumbers(overlay, cells);
		IJ.run("Merge Channels...", "c1=["+split[0].getTitle()+"] c2=["+split[1].getTitle()+
				"] c3=["+split[2].getTitle()+"] c4=[Overlay] create");
		return WindowManager.getCurrentImage();
	}

	private void drawNumbers(ImagePlus ProjectedWindow, Cell[] cells) {
		ImageProcessor ip = ProjectedWindow.getProcessor();
		Font font = new Font("SansSerif", Font.BOLD, 30);
		ip.setFont(font);
		ip.setColor(Color.white);
		for(int i = 0; i< cells.length;i++) {
			String cellnumber = String.valueOf(i + 1);
			ip.draw(cells[i].roi);
			ip.drawString(cellnumber, (int) cells[i].roi.getContourCentroid()[0], (int) cells[i].roi.getContourCentroid()[1]);
			ProjectedWindow.updateAndDraw();
			for(int j = 0; j< cells[i].spots.size();j++) {
				ip.draw(cells[i].spots.get(j));
				ProjectedWindow.updateAndDraw();
			}
		}
	}

	private String[] splitFilename(String filename){
		int indexEnd = filename.length()-1;
		if(filename.contains("-Airy")){
			indexEnd = filename.indexOf("-Airy");
		}
		return new String[]{filename.substring(0, indexEnd-2), filename.substring(indexEnd-1, indexEnd)};
	}

	public void MakeResults(String[] titles, double[][] results) throws IOException {
		Date date = new Date(); // This object contains the current date value
		SimpleDateFormat formatter = new SimpleDateFormat("dd-MM-yyyy, hh:mm:ss");
		String[] splitFilename = splitFilename(filename);
		String CreateName = Paths.get(String.valueOf(filePath), splitFilename[0] + "_Results.csv").toString();
		boolean exists = Files.exists(Paths.get(CreateName));
		try {
			FileWriter fileWriter = new FileWriter(CreateName, true);
			BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);
			//bufferedWriter.newLine();
			if(!exists) {
				bufferedWriter.write(formatter.format(date));
				bufferedWriter.newLine();
				//bufferedWriter.newLine();

				String titleString = "File, Cell, ";
				for (String title : titles) {
					titleString = titleString + title + ",";
				}
				bufferedWriter.write(titleString);//write header 1
				bufferedWriter.newLine();
			}
			for (int i =0; i < results.length; i++) {//for each slice create and write the output string
				bufferedWriter.newLine();
				String resultsString = splitFilename[1]+","+ (i+1)+",";
				for (double result:results[i]){
					resultsString = resultsString + result + ",";
				}
				bufferedWriter.write(resultsString);
			}
			bufferedWriter.close();
		} catch (IOException ex) {
			System.out.println(
					"Error writing to file '" + CreateName + "'");
		}
	}

    public static void main(final String... args) throws Exception {
        // create the ImageJ application context with all available services
        final ImageJ ij = new ImageJ();
        ij.ui().showUI();
        ij.command().run(Sara_Plugin.class, true);
    }

}

