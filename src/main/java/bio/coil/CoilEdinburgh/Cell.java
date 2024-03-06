package bio.coil.CoilEdinburgh;

import ij.ImagePlus;
import ij.gui.Roi;
import ij.process.ImageProcessor;

import java.util.ArrayList;

public class Cell {
    Roi roi;
    ArrayList<Roi> spots;

    public Cell(Roi roi) {
        this.roi = roi;
        this.spots = new ArrayList<>();
    }

    public void addSpot(Roi newRoi){
        spots.add(newRoi);
    }

    public double totalArea() { return  roi.getStatistics().area;}

    public double spotsArea(){
        double spotArea = 0;
        for(int i = 0; i< spots.size();i++){
            spotArea = spotArea + spots.get(i).getStatistics().area;
        }
        return spotArea;
    }

    public double cellArea(){
        return totalArea() - spotsArea();
    }

    public double spotsIntensity(ImagePlus imp){
        ImageProcessor ip = imp.getProcessor();
        double spotIntensity = 0;
        for(int i = 0; i< spots.size();i++){
            ip.setRoi(spots.get(i));
            spotIntensity = spotIntensity + ip.getStatistics().mean*ip.getStatistics().area;
        }
        return spotIntensity/spotsArea();
    }

    public double totalIntensity(ImagePlus imp){
        ImageProcessor ip = imp.getProcessor();
        ip.setRoi(roi);
        return ip.getStatistics().mean;
    }

    public double cellIntensity(ImagePlus imp){
        double totalInt = totalIntensity(imp)*totalArea();
        double spotsInt = spotsIntensity(imp)*spotsArea();
        return (totalInt - spotsInt)/cellArea();
    }

    public int nSpots(){
        return spots.size();
    }

    public boolean containsSpot(Roi newRoi){
        int x = (int) newRoi.getContourCentroid()[0];
        int y = (int) newRoi.getContourCentroid()[1];
        if(roi.contains(x,y)){
            return true;
        }
        return false;
    }

}
