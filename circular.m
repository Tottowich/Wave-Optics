
% BASIC IMAGE ANALYSIS SCRIPT
% This script is intended as an example to help you analyze the pictures you have acquired
% during the "Diffraction, interference and Fourier filtering" lab. Notice that you cannot
% simply run the algorithm "as it is", you will need to change a few things to make it
% work with your images (therefore you will have to understand it, contact your lab
% supervisor if you have questions). You can also write your own algorithm making
% different steps if you want (this here is just for inspiration).

close all; clear all;
%%
lamb = 632.8*10^-9;
b = 1.5*10^-4;
L = 0.972;
minima = [1.220,2.233,3.238];
new_b = @(y,m) abs(lamb.*minima(m)*L./y);%sind(th);
new_rel_err = @(y,m) (new_b(y,m)-b)./b;
rel_err = @(app,real) abs((app-real)./real);
%% Import the image to Matlab
path = 'Images/CircularApature/circularApature.jpg';
channel = 1;
rotation = 2.7;
sec_x = 1500:3000;
sec_y = 1300:1310;
baseline_sign = 0; % No baseline subtraction or addition
shift = (604+929)/2;
S = (135/5);
L = 0.972;
show_steps = true;
mm_limit = [-100,100];
deg_limit = [-2,2];
[xAxisDeg, xAxisMm, data3] = image_process(path,channel,rotation,sec_x, ...
               sec_y,baseline_sign,shift,S,L,show_steps,mm_limit,deg_limit);

%%
% Diffraction : b
visualize = true;
[PKS,LOCS] = find_min(data3,xAxisMm,1,inf,visualize);
y_diffraction = LOCS*10^-3;
if visualize
    grid on;
    xlabel('mm - displacement');
    ylabel('Normalized signal');
    title('Normalized signal vs angle with Min peaks');
    xlim([-30 , 30]);
    ylim([0 , 1]);
end
if show_steps
    plot(xAxisDeg,data3);
    hold on
    degs = atand(y_diffraction/L);
    scatter(degs,PKS,'r^','filled');
    hold off
    grid on;
    xlabel('Angle (°)');
    ylabel('Normalized signal');
    title('Normalized signal vs angle with Min peaks');
    xlim([-2 , 2]);
    ylim([0 , 1]);
end
left_hand = sum(y_diffraction<0);
right_hand = sum(y_diffraction>=0);
m = [-left_hand:-1 1:right_hand];
mabs = abs(m);
valid_minima = mabs<=3; % Could only find the m-values for minima 1,2 & 3.
y_diffraction = y_diffraction(valid_minima);
mabs = mabs(valid_minima);
slit_widths = new_b(y_diffraction,mabs);
d_mean = mean(slit_widths);
d_std = std(slit_widths);
%% Diffraction Error:
N_diff = length(slit_widths);
dL = 5*10^-3;
dP = 25;
dy = dP/(S*10^3);
mabs_circ = minima(mabs);
dd =sqrt((dL*lamb*mabs_circ./y_diffraction).^2+(mabs_circ.*L*dy*lamb./(y_diffraction.^2)).^2);
dd_mean = mean(dd);
dd_std = std(dd);
fprintf("Diffraction Circular\n");
fprintf("d = %.3e +- %.3e\n",d_mean,dd_mean);
fprintf("d_std = %.3e \n",d_std);
fprintf("Error std = %.3e \n",dd_std);
fprintf("Relative Error = %.3e \n\n",rel_err(d_mean,b));



%% Appendix
%{
%%
myImage1 = imread('Images/CircularApature/circularApature.jpg');

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
myImage4 = imrotate(myImage3, 2.7);
if show_steps
    figure(4)
    imagesc(myImage4);
    title('Step 4: after rotation');
end

%% Cut a region of interest within the image, where the diffraction pattern is contained
%myImage5 = myImage4(1450:1550, 1000:3000); Double Slit
myImage5 = myImage4(1300:1310, 1500:3000);
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
data2 = data1;% + baseline;
if show_steps
    figure(7);
    plot(data2);
end

%% Calibrate the pixel with space. Using the mm paper in the image I see that 420 px
% correspond to 50 mm. Therefore, 8.4 px/mm (I leave the error calculations to you :p ).
% Notice that I have centered the 0 at the main peak.
xAxisMm = ((1 : length(data2)) - (604+929)/2) / (135/5); % Single Slit
S = (135/5);
if show_steps
    figure;
    plot(xAxisMm, data2);
    xlabel('Distance (mm)');
    ylabel('Counts (a.u.)');
    title('Step 7: Calibrate the data from signal-vs-pixels into signal-vs-distance');
end

%% Calibrate the space with angle. For this example, I used a distance between the
% diffraction slit and the screen of 1 meter. Normalize and plot the results nicely.
xAxisDeg = atand(xAxisMm / 1000*L);
f = figure;
data3 = data2 / max(data2);
plot(xAxisDeg, data3);
xlim([-2 , 2]);
ylim([0 , 1]);
xlabel('Angle (°)');
ylabel('Normalized signal');
title('Final step: normalized signal vs angle');
grid on;
% b
%}

