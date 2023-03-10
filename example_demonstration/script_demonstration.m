% BASIC IMAGE ANALYSIS SCRIPT
% This script is intended as an example to help you analyze the pictures you have acquired
% during the "Diffraction, interference and Fourier filtering" lab. Notice that you cannot
% simply run the algorithm "as it is", you will need to change a few things to make it
% work with your images (therefore you will have to understand it, contact your lab
% supervisor if you have questions). You can also write your own algorithm making
% different steps if you want (this here is just for inspiration).

close all; clear all;

%% Import the image to Matlab
myImage1 = imread('doubleSlit.jpg');

%% Display the image as a figure
figure(1)
imagesc(myImage1);
title('Step 1: original image');

%% Select the RGB color channel to use (in this case the blue, as it is less saturated)
myImage2 = myImage1(:,:,3);
figure(2)
imagesc(myImage2);
title('Step2: using the blue channel only');

%% Apply a median filter to reduce hotspots
myImage3 = medfilt2(myImage2);
figure(3)
imagesc(myImage3);
title('Step 3: after median filtering');

%% Rotate the image so that the diffraction pattern is horizontal
myImage4 = imrotate(myImage3, 2);
figure(4)
imagesc(myImage4);
title('Step 4: after rotation');

%% Cut a region of interest within the image, where the diffraction pattern is contained
myImage5 = myImage4(388:491, 50:1270);
figure(5)
imagesc(myImage5);

%% Integrate over the first dimension (i.e. vertically) to get the total counts per pixel.
data1 = sum(myImage5);
figure(6);
plot(data1);
title('Step 6: vertically integrate the signal');

%% Substract the background baseline
baseline = 2800 + (3600/length(data1)) * (1:length(data1));
data2 = data1 - baseline;
figure(7);
plot(data2);

%% Calibrate the pixel with space. Using the mm paper in the image I see that 420 px
% correspond to 50 mm. Therefore, 8.4 px/mm (I leave the error calculations to you :p ).
% Notice that I have centered the 0 at the main peak.
xAxisMm = ((1 : length(data2)) - 702.575) / 8.4;
figure;
plot(xAxisMm, data2);
xlabel('Distance (mm)');
ylabel('Counts (a.u.)');
title('Step 7: Calibrate the data from signal-vs-pixels into signal-vs-distance');

%% Calibrate the space with angle. For this example, I used a distance between the
% diffraction slit and the screen of 1 meter. Normalize and plot the results nicely.
xAxisDeg = atand(xAxisMm / 1000);
f = figure;
data3 = data2 / max(data2);
plot(xAxisDeg, data3);
xlim([-4 , 4]);
ylim([0 , 1]);
xlabel('Angle (°)');
ylabel('Normalized signal');
title('Final step: normalized signal vs angle');
grid on;
