# MRI_Muscle_Cross-Sectional_Area_PTOA
Matlab program for calculating extensor and flexor muscle cross-sectional areas from MRI 3T T1-FFE images of thigh muscles in subjects with first-time ACL rupture with and without concomitant meniscal injury. 

The M-file reads a T1 FFE image file, crops the image, and thresholds the image. The user then selects a pixel within the subcutaneous fat and femur in the left and right thighs. The program separately finds all connected muscle and subcutaneous fat tissues. The femur is filled and used to exclude the femur and marrow from the noncontractile elements within the muscles.

The program prompts the user to create a polygon region of interest around the flexor muscles, and then the extensor muscles. This is used to divide the muscles into two main muscle group areas. A plot of just the muscles is used to verify the division of the muscle into extensors and flexors before continuing with the program.

The cross-sectional areas for the muscles, subcutaneous fat and noncontractile elements are displayed in the command window. The program also outputs the results to a spreadsheet, mthresh_ptoa_*.xlsx, in the directory PTOA_Muscle_CSA_* where * is the visit name (Baseline, Y1 [Year 1], or Y2 [Year 2]).

Plots of the raw image, threshold histogram, muscles (extensors, flexors and total), subcutaneous fat and noncontractile elements are written to a Postscript (older versions of Matlab) or PDF (newer versions of Matlab) file mthresh_ptoa_*_**.ps/.pdf, where * is the visit name (Baseline, Y1, or Y2) and ** is the subject number, into the results directory. Note that both the results spreadsheet and the plots are all in the same directory.

See Polygon_ROI_Guide_PTOA.pdf for tips on creating the polygon ROI. See musc_thresh_ptoa_Guide.pdf for a guide to using the program.

See comments in musc_thresh_ptoa.m for more information.

Notes on the use of the program.

    This program is for the PTOA study. The directory structure must include subdirectories for each visit (PTOA * MRI Files, where * is the visit name [Baseline, Y1, or Y2]). The directory structure is used to identify the visit of the MRI images.

    M-file function roi_mov.m must be in the current path or directory. I recommend putting both musc_thresh-ptoa.m and roi.mov.m in the same directory.

    The output MS-Excel spreadsheet, mthresh_ptoa_*.xlsx, can NOT be open in another program (e.g. MS-Excel) while using this program.

    Running the program for the same MRI image will result in duplicate data for that image in the spreadsheet. The spreadsheet should be checked for any duplicate data before any statistical analyzes.
