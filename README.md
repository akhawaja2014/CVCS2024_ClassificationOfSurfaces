# Self-Supervised Classification of Surfaces using Reflectance Transformation Imaging
CVCS 2024 Code





## Dependencies

The script is written in MATLAB and requires the following:

- MATLAB R2021a or later
- Image Processing Toolbox
- Neural Network Toolbox (for SOM)
- Statistics and Machine Learning Toolbox (for k-means)




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


