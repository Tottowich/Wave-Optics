% BASIC IMAGE ANALYSIS SCRIPT
% This script is intended as an example to help you analyze the pictures you have acquired
% during the "Diffraction, interference and Fourier filtering" lab. Notice that you cannot
% simply run the algorithm "as it is", you will need to change a few things to make it
% work with your images (therefore you will have to understand it, contact your lab
% supervisor if you have questions). You can also write your own algorithm making
% different steps if you want (this here is just for inspiration).

close all; clear all;
%%
lamb = 633*10^-9;
b = 4*10^-5;
a = 1.25*10^-4;
d = 1/(144*10^3);
L = 0.22;
new_d_max = @(y,m) lamb*m*L./y;
new_b_min = @(y,m) lamb*m*L./y;
new_a_min = @(y,m) lamb*(m-1/2)*L./y;
rel_err = @(y,m) (new_b(y,m)-b)/b;
show_steps = true;
%% Import the image to Matlab
myImage1 = imread('Images/Grating/grating.jpg');

%% Display the image as a figure
if show_steps
    figure(1)
    imagesc(myImage1);
    title('Step 1: original image');
end

%% Select the RGB color channel to use (in this case the blue, as it is less saturated)
myImage2 = myImage1(:,:,1);
if show_steps
    figure(2)
    imagesc(myImage2);
    title('Step2: using the blue channel only');
end

%% Apply a median filter to reduce hotspots
myImage3 = medfilt2(myImage2);
if show_steps
    figure(3)
    imagesc(myImage3);
    title('Step 3: after median filtering');
end

%% Rotate the image so that the diffraction pattern is horizontal
myImage4 = imrotate(myImage3, 4);
if show_steps
    figure(4)
    imagesc(myImage4);
    title('Step 4: after rotation');
end

%% Cut a region of interest within the image, where the diffraction pattern is contained
myImage5 = myImage4(1000:1600, :);%Double Slit
%myImage5 = myImage4(1550:1700, :);
if show_steps
    figure(5)
    imagesc(myImage5);
end

%% Integrate over the first dimension (i.e. vertically) to get the total counts per pixel.
data1 = sum(myImage5);
if show_steps
    figure(6);
    plot(data1);
    title('Step 6: vertically integrate the signal');
end

%% Substract the background baseline
baseline = 2800 + (3600/length(data1)) * (1:length(data1));
data2 = data1 - baseline;
if show_steps
    figure(7);
    plot(data2);
end

%% Calibrate the pixel with space. Using the mm paper in the image I see that 420 px
% correspond to 50 mm. Therefore, 8.4 px/mm (I leave the error calculations to you :p ).
% Notice that I have centered the 0 at the main peak.
S = (833-616)/5; %px/mm
xAxisMm = ((1 : length(data2)) - (2590+1427)/2) / S; % Single Slit
if show_steps
    figure;
    plot(xAxisMm, data2);
    xlabel('Distance (mm)');
    ylabel('Counts (a.u.)');
    title('Step 7: Calibrate the data from signal-vs-pixels into signal-vs-distance');
end

%% Calibrate the space with angle. For this example, I used a distance between the
% diffraction slit and the screen of 1 meter. Normalize and plot the results nicely.
xAxisDeg = atand(xAxisMm / (L*1000));
f = figure;
data3 = data2 / max(data2);
plot(xAxisDeg, data3);
xlim([-11 , 11]);
ylim([0 , 1]);
xlabel('Angle (°)');
ylabel('Normalized signal');
title('Final step: normalized signal vs angle');
grid on;
% b
%%
visualize = true;
% Interference : a
[PKS,LOCS] = find_max(data3,xAxisMm,1,10,visualize);
y_interference = LOCS*10^-3;
if visualize
    grid on;
    xlabel('mm - displacement');
    ylabel('Normalized signal');
    title('Normalized signal vs angle with Min peaks');
    xlim([-60 , 60]);
    ylim([0 , 1]);
end
left_hand = sum(y_interference<0);
right_hand = sum(y_interference>=0);
m = [-left_hand:-1 1:right_hand];
grating_constants = new_d_max(y_interference,m);
d_mean = mean(grating_constants);
d_std = std(grating_constants);
%% Interference Error:
N_inter = length(slit_spacings);
dL = 5*10^-3;
dP = 25;

dy = dP/(S*10^3);

mabs = abs(m);
da = sqrt((dL*lamb*(mabs-1/2)./y_interference).^2+((mabs-1/2).*L*dy*lamb./(y_interference.^2)).^2);
da_mean = mean(da);
da_std = std(da);
fprintf("Interference\n");
fprintf("a = %.3e +- %.3e\n",a_mean,da_mean);
fprintf("a_std = %.3e \n",a_std);
fprintf("Error std = %.3e \n\n",da_std);
%%
% Diffraction : b
%{
[PKS,LOCS] = find_min(data3,xAxisMm,2.5,inf,visualize);
y_diffraction = LOCS*10^-3;
if visualize
    grid on;
    xlabel('mm - displacement');
    ylabel('Normalized signal');
    title('Normalized signal vs angle with Min peaks');
    xlim([-60 , 60]);
    ylim([0 , 1]);
end
left_hand = sum(y_diffraction<0);
right_hand = sum(y_diffraction>=0);
m_diff = [-left_hand:-1 1:right_hand];
slit_widths = new_d_max(y_diffraction,m_diff);
d_mean = mean(slit_widths);
d_std = std(slit_widths);
%% Diffraction Error:
N_diff = length(slit_widths);
mabs = abs(m_diff);
db =sqrt((dL*lamb*mabs./y_diffraction).^2+(mabs.*L*dy*lamb./(y_diffraction.^2)).^2);
db_mean = mean(db);
db_std = std(db);
fprintf("Diffraction\n");
fprintf("b = %.3e +- %.3e\n",b_mean,db_mean);
fprintf("b_std = %.3e \n",b_std);
fprintf("Error std = %.3e \n",db_std);
%}



