MitoTouch Matlab code
------------------------------------------------------------

Our images taken with the Nikon Spinning Disk confocal microscope are classified according to the following rationale:
Cell type > stress > Date of experiment > microscopy field

Tree definition:

-Cell type: corresponds to the cell line used

-Stress: corresponds to stress applied to cells:
-Nostress (basal condition)
-SD: serum deprivation (elongation of mitochondria)
-UVB500: UVB irradiation at 500mJ/cm^-2 (fission of mitochondria)

-Date: Date the image was taken

-Microscopy field: Images taken at a given position

In each field folder, there are 4 images in .tiff format corresponding to images taken at the same position on the channels:
- Green (CellMask Green): cellular contour
- Blue (Hoechst): nucleus
- Red (MitoTracker Red): mitochondria
- Merge: superposition of the three previous images.

Images can be viewed using Fiji ImageJ software --> https://fiji.sc/
------------------------------------------------------------

MitoTouch features:

All files ending in "MERGE.tif" that are included in a folder named "Champ - Copie (X) " will be processed in the selectionned tree directory.

It is higly recommended to test MITOTOUCH parameters on 1 image of the the selectionned tree directory before processing analysis.

Load test image:
Allows you to load and display a test image, in order to refine the selection parameters

Test: 
Allows you to reload the test image after changing the parameters

Nuclei threshold :
Ability of the software to detect the nuclei of the cells, which will be surrounded in yellow. Increasing it means that the selection will be stricter.
A value that is too low can lead to selecting elements that are not nuclei, a value that is too high can lead to not selecting nuclei.

Cell threshold : 
Ability of the software to detect the contour of the cells, which will be surrounded in pink. Increasing it means that the selection will be stricter.
A value that is too low can lead to selecting elements that are not part of the cells, a value that is too high can lead to not selecting certain cells. 

Cell erode : 
Determines the contour accuracy of the cell contour selection. A higher value will have finer contours 

Filter border cell :
Allows to select not complete cells that must be present in the image for the software to select it. 
A value of 0% will select all partial cells, while a value of 100% will only select cells completely in the image.

Export mito based data : Check this option if you want to have parameters for each mitochondria instead of a mean of mitochondrial parameters per cell.

Display : Allows to check the selection on the test image in order to check the accuracy of the parameters

Original: Displays the original image	
Nuclei: Displays the selected nuclei 	
Cell: Displays the contour of the selected cells
Mito : 	Displays mitochondria
Skeleton: Displays the skeleton of the mitochondria 	 
Final: Displays all elements

Load test image:
Allows you to load and display a test image, in order to refine the selection parameters

Test: 
Allows you to reload the test image after changing the parameters

Run : 
Allows you to perform processing under the selected folder. All files ending in "MERGE.tif" that are included in a folder named "Champ - Copie (X) " will be selected for processing
CAUTION, this also includes the entire descending file tree of the selected folder.

Contact:
-----------------------
abdel.aouacheria@umontpellier.fr & victor.racine@quantacell.com 

Parallel :
Distributes the operation over the different cores of the computer to increase processing speed
CAUTION, due to the power required, a minimum of 16 GB of RAM is recommended. 

At the end of the analysis a "result" folder in format will be created in the selectionned tree directory. It will be composed of all the images processed by MITOTOUCH and a .csv file containing all the numerical results.
