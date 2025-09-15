# Microchip2DPostParticleTracking

## MATLAB Project for 2-D post-experimental, blob and object-aligned detection and tracking of microchips using a monocular camera setup

This code was run using MATLAB R2024b (using the Java-based desktop and graphics system) on Windows 10/11. File paths need to be modified for Unix-based platforms (macOS, Linux). MATLAB versions R2025a and above, which uses a new environment based on HTML and JavaScript, may cause unintended graphical issues.

## How to setup the MATLAB Project
Make sure Git and Git-LFS are installed. The Git repository uses Git LFS to store/manage the large sample video file. Ensure Git LFS is appropriately installed to checkout the video file; otherwise, the video file will show up as a small pointer file.
1. Open MATLAB (current version tested for MATLAB R2024b).
2. Navigate to the 'HOME' tab, and then click 'New' > 'Project', 'From Git'.
    1. Ensure 'Source control tool' uses Git
    2. For the 'Repository path', enter: "https://github.com/raykyeung/Microchip2DPostParticleTracking"
    3. For the 'Sandbox", select a suitable location for the local repository (ex. "C:\Users\USERNAME\LocalGitRepo\https://github.com/raykyeung/Microchip2DPostParticleTracking")
    4. Click 'Retrieve' and allow MATLAB to pull all files. Project dependencies and tracking are handled by the .prj file and 'resources' folder.

## Structure of the MATLAB Project
### Main scripts (under BatchTestScripts)

**Batch1pTracking_BlobGaussianMixture.m** - Main script to perform particle tracking using blob detection and a Gaussian mixture model. Briefly, the script performs the following:
1. Calls B2KVideoSummaryBGM.m to summarize videos and identify the start and end frames to analyze.
2. Calls B2KCoordinateSystemBGM.m to perform microchannel wall detection and image calibration.
3. Calls B2KParticleTrackingBGM.m to perform object-oriented particle tracking.
4. Export all data and figures (default under 'Downloads')
5. Save log file containing dependency information and required program files

**Batch1pTracking_ObjectOriented.m** - Main script to perform particle tracking using blob detection and a Gaussian mixture model. Briefly, the script performs the following:
1. Calls B2KVideoSummaryOO.m to summarize videos and identify the start and end frames to analyze.
2. Calls B2KCoordinateSystemOO.m to perform microchannel wall detection and image calibration.
3. Calls B2KParticleTrackingOO.m to perform object-oriented particle tracking.
4. Export all data and figures (default under 'Downloads')
5. Save log file containing dependency information and required program files

### Functions
**B2KDistbw2Lines.m** - Function to evaluate the distance between two lines in order to determine the distance between detected channel walls.
**B2KMinBoundParallelogram.m** - Function to evaluate the minimum enclosing parallelogram.
**subpixelEdges.m** - Function to perform subpixel edge detection based on the partial area effect. Used to accurately locate the edges of the particle.
**B2Khoughlines.m** - Function to perform Hough transform to evaluate the Hough lines corresponding to the channel wall.

**B2KExpParam.m** - Function to evaluate experimental parameters from video name.
**B2KFlags.m** and **B2KFlagsBGM** - Functions to provide booleans for different run settings.
**B2KMyFigures.m** - Function to assist with figure handle generation.
**B2KVideoSummaryOO.m** - Function to perform video summarization based on background subtraction using defined image statistics.
**B2KVideoSummaryBGM.m** - Function to perform video summarization based on a robust Gaussian mixture model.
**B2KRotateandCrop.m** - Function to handle rotation and cropping for videos with slightly angled channels.
**B2KCoordinateSystemOO.m** and **B2KCoordinateSystemBGM.m** - Functions to perform microchannel wall detection and image calibration. Maps data to defined world coordinates.
**B2KParticleTrackingOO.m** - Function to perform 2-D particle tracking based on subpixel edge detection and a mapped minimum enclosing parallelogram. Uses a constant-velocity Kalman filter for particle tracking prediction.
**B2KParticleTrackingBGM.m** - Function to perform 2-D particle tracking based on a robust Gaussian mixture model for foreground detection and blob detection of the particle. Uses a constant-velocity Kalman filter for particle tracking prediction.

### Classes
#### +B2KConstants
Classes for storing particle and fluid properties:
* @Acetonitrile.m
* @Ethanol.m
* @Glycerol.m
* @Isopropanol.m
* @Methanol.m
* @pChip.m
* @Water.m

### Data Videos
**1p_Water_R1_221026_H1.000_W1.500_25mL-min_1500FPS_2.** - Sample video of a flowing microchip in a straight, rectangular microfluidic channel (top-down, fronto-parallel view; channel height of 1.00 mm; channel width of 1.50 mm; water solvent; 25 mL/min flow rate; 1500 FPS capture). 

## Copyright
Copyright (c) 2025 Raymond K. Yeung

This product includes software developed at the University of California, Riverside (https://www.ucr.edu/) for NSF CPS Grant Award #1740052 (https://www.nsf.gov/awardsearch/showAward?AWD_ID=1740052).

This product is licensed under the terms of the Apache 2.0 license.

The project contains the following code and derivatives of 3rd-party resources:
* Function, Distbw2Lines.m, Version 1.0.0 by Gravish based on Dan Sunday's Practical Geometry Algorithms written in C++, under MIT license, for evaluating distance between two lines.
* Function suite, MinBoundSuite, Version 1.2.0.0 by John D'Errico from MATLAB File Exchange (https://www.mathworks.com/matlabcentral/fileexchange/34767-a-suite-of-minimal-bounding-objects?s_tid=srchtitle) under MIT License, for finding the minimum bounding parallelogram around detected edges.
* Function suite, Accurate subpixel edge location, Version 2.13.3 by Agustin Trujillo-Pino from MATLAB File Exchange (https://www.mathworks.com/matlabcentral/fileexchange/48908-accurate-subpixel-edge-location?s_tid=srchtitle) under MIT License, for finding the subpixel edges of the particle.
* Function, houghline, MATLAB R2024b, by MathWorks based on Rafael C. Gonzalez, Richard E. Woods, Steven L. Eddins, "Digital Image Processing Using MATLAB", Prentice Hall, 2003.