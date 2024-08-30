# Self-Supervised Classification of Surfaces using Reflectance Transformation Imaging
This repository contains a MATLAB script that performs Reflectance Transformation Imaging (RTI) analysis on a set of images, extracts Region of Interest (ROI) signals, and applies clustering techniques like Self-Organizing Maps (SOM) and k-means to classify and visualize these signals.

## Table of Contents

- [Introduction](#introduction)
- [Dependencies](#dependencies)
- [Usage](#usage)
- [File Structure](#file-structure)
- [Code Overview](#code-overview)
- [Clustering and Visualization](#clustering-and-visualization)
- [License](#license)

## Introduction

This script is designed to work with RTI datasets, typically containing multiple images that capture surface details under varying lighting conditions. The script reads these images, allows users to define a Region of Interest (ROI), extracts pixel intensity signals from the ROI across the image set, and applies clustering techniques to identify patterns and group similar signals.

## Dependencies

The script is written in MATLAB and requires the following:

- MATLAB R2021a or later
- Image Processing Toolbox
- Neural Network Toolbox (for SOM)
- Statistics and Machine Learning Toolbox (for k-means)

## Usage

1. **Set Paths**: Update the paths in the script to point to your dataset:
   - `testimage_path`: Path to the initial image to select ROI.
   - `RTI_acq`: Path to the folder containing RTI images.
   - `lpfile_path`: Path to the light position file (not used directly in this script but kept for reference).

2. **Run the Script**: Execute the MATLAB script. It will prompt you to select a ROI interactively. The script then reads all images in the specified folder, extracts signals from the ROI, and applies clustering algorithms.

3. **Visualize Results**: The script provides various visualizations including signal plots for each cluster and overlay masks on the original image based on SOM and k-means clustering results.

## File Structure

- `main_script.m`: The main MATLAB script for processing RTI images and clustering ROI signals.

## Code Overview

### Paths for RTI Data

The following MATLAB code sets the paths for various data sources:

```matlab
% PATHS for real data
testimage_path = 'Path to the test image used for processing';
RTI_acq = 'Directory path containing RTI acquisition data.';
lpfile_path = 'Path to the .lp file with light directions';
```

 The following code snippet is used for reading and displaying an image, allowing interactive selection of a region of interest (ROI) via a rectangle, and then creating and displaying a binary mask based on that ROI.

```matlab
I = imread(testimage_path);  % Read the image file
imshow(I);  % Display the image
h = imrect;  % Interactively draw a rectangle
% Get the position of the rectangle
roiPosition = wait(h);
% Create a binary mask for the ROI
mask = false(size(I, 1), size(I, 2));
mask(round(roiPosition(2)):round(roiPosition(2)+roiPosition(4)), round(roiPosition(1)):round(roiPosition(1)+roiPosition(3))) = true;
% Display the mask (optional)
figure;
imshow(mask);
```


