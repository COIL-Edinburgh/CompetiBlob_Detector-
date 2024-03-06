##Sara_Plugin

A plugin to investigate the relationship between Colocalisation of the protein of interest with DNA and 
transfection levels of a second protein. 

###Input:
The plugin takes as input a folder of three channel, Z-stack .czi files;

   - Ch0 = Protein of interest
   - Ch1 = Transfected protein
   - Ch0 = Nuclear Stain

The plugin requires the path to the cellpose environment and the cellpose models to be used to be provided. We tested
using the nucleitorch_3 model.

###Method:

- The plugin performs an "Average" Z-projection on all channels.
- The plugin uses the cellpose model with a diameter of 14 microns  to segment the nuclei in Ch2. 
- A 25 pixel rolling-ball background subtraction is applied to Ch2 then a Yen threshold is applied to each nucleus to 
  segment out bright spots. 
- The area and intensity of the nucleus, the bright spots within the nucleus and the background within the nucleus are 
reported for Ch0 and Ch2.
- The intensity of each nucleus in Ch1 (level of transfection) is reported.

###Output:

For each file processed a filename_Overlay.tif will be saved with the detected spots and numbered nuclei drawn in a 4th 
channel.

If the input files are named according to 'Genotype_N-Airyscan Processing.czi' the plugin will create results files 
labelled Genotype_Results.csv. With each repetition (N) numbered in the file. If the files do not contain the String 
"-Airy" a separate results file will be created for each file. The _results .csv file contains for each nucleus;

File | Cell | Cell Area | Cell Intensity Ch0 | Cell Intensity Ch1 | Cell Intensity Ch2 | Spots Area | 
Spots Intensity Ch0 | Spots Intensity Ch2 | Background Area | Background Intensity Ch0 | Background Intensity Ch2 |




