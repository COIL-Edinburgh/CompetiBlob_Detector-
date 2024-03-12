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
import ij.io.FileSaver;
import ij.process.ImageStatistics;
import org.apache.commons.io.FileUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class Cellpose_Wrapper {
    ImagePlus imp;
    final String modelpath;
    final String envPath;
    final int diameter;

    public Cellpose_Wrapper(String modelpath, String envPath, int diameter, ImagePlus imp) {
        this.modelpath = modelpath;
        this.envPath = envPath;
        this.diameter = diameter;
        this.imp = imp;
    }

    TempDirectory temp;

    public ImagePlus run(boolean ROIs){

    temp = new TempDirectory("temp");
    purgeDirectory(temp.getPath().toFile());

    saveImageWithoutUI(imp, temp.getPath().toString());

    //Construct Commands for cellpose
    List<String> cmd = getCmd(diameter,modelpath);

    //Send to cellpose
    try {
        CmdTest(cmd);//.toString().replace(",", ""));
    } catch (Exception e) {
        e.printStackTrace();
    }

    //return image
    ImagePlus masks = IJ.openImage(Paths.get(temp.getPath().toString(),"Temp_cp_masks.tif").toString());

    if(ROIs) {
        masks.show();
        getROIsfromMask(masks);
        //getRois(masks);
    }

    temp.deleteOnExit();
        return masks;

}

private void getROIsfromMask(ImagePlus mask) {
    mask.show();
    ImageStatistics stats = mask.getStatistics();
    //For each ROI (intensity per cell mask is +1 to intensity
    for (int i = 1; i < stats.max + 1; i++) {
        //Set the threshold for the cell and use analyse particles to add to ROI manager
        IJ.setThreshold(mask, i, i);
        IJ.run(mask, "Analyze Particles...", "exclude add");
    }
}

private  List<String> getCmd(int diameter, String model){
    List<String> cmd = new ArrayList<>();
    List<String> start_cmd=  Arrays.asList("cmd.exe", "/C");
    cmd.addAll(start_cmd);
    List<String> condaCmd = Arrays.asList("CALL", "conda.bat", "activate", envPath);
    cmd.addAll(condaCmd);
    cmd.add("&");
    List<String> cellpose_args_cmd = Arrays.asList("python", "-Xutf8", "-m", "cellpose");
    cmd.addAll(cellpose_args_cmd);
    List<String> options = Arrays.asList("--diameter",diameter+"","--verbose", "--pretrained_model", model,"--save_tif","--use_gpu", "--dir",temp.getPath().toString());
    cmd.addAll(options);
    if (imp.getNChannels()==2){
        List<String> options2 = Arrays.asList("--chan", 2+"","--chan2", 1+"");
        cmd.addAll(options2);
    }
    return cmd;
}

void purgeDirectory(File dir) {
    for (File file: dir.listFiles()) {
        if (file.isDirectory())
            purgeDirectory(file);
        file.delete();
    }
}

public class TempDirectory {
    final Path path;

    public TempDirectory(String prefix) {
        try {
            //Path currentRelativePath = Paths.get("");
            path = Files.createTempDirectory(prefix);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public Path getPath() {
        return path;
    }

    public void deleteOnExit() {
         try {
                FileUtils.deleteDirectory(path.toFile());
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }

    }



public void CmdTest(List<String> cmd) throws Exception {
    ProcessBuilder builder = new ProcessBuilder(cmd);
    builder.redirectErrorStream(true);
    Process p = builder.start();
    BufferedReader r = new BufferedReader(new InputStreamReader(p.getInputStream()));
    String line;
        while (true) {
            line = r.readLine();
            if (line == null) { break; }
            System.out.println(line);
        }
    }

private static void saveImageWithoutUI(ImagePlus imagePlus, String filePath) {
    FileSaver fileSaver = new FileSaver(imagePlus);
    fileSaver.saveAsTiff(Paths.get(filePath,"Temp.tif").toString() ); // You can use other formats like PNG, TIFF, etc.
    // Optionally, you can also close the image after saving
    //imagePlus.close();
}
}
