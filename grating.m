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
%%
path = 'Images/Grating/grating.jpg';
channel = 1;
rotation = 4;
sec_x = "all";
sec_y = 1000:1600;
baseline_sign = 0; % No baseline subtraction or addition
shift = (2590+1427)/2;
S = (833-616)/5; %px/mm
L = 0.22;
show_steps = true;
mm_limit = [-60,60];
deg_limit = [-12,12];
[xAxisDeg, xAxisMm, data3] = image_process(path,channel,rotation,sec_x, ...
               sec_y,baseline_sign,shift,S,L,show_steps,mm_limit,deg_limit);
%% Import the image to Matlab

myImage1 = imread('Images/Grating/grating.jpg');

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

if show_steps
    plot(xAxisDeg,data3);
    hold on
    degs = atand(y_interference/L);
    scatter(degs,PKS,'r^','filled');
    hold off
    grid on;
    xlabel('Angle (°)');
    ylabel('Normalized signal');
    title('Normalized signal vs angle with Min peaks');
    xlim([-12 , 12]);
    ylim([0 , 1]);
end

left_hand = sum(y_interference<0);
right_hand = sum(y_interference>=0);
m = [-left_hand:-1 1:right_hand];
grating_constants = new_d_max(y_interference,m);
d_mean = mean(grating_constants);
d_std = std(grating_constants);
%% Interference Error:
N_inter = length(grating_constants);
dL = 5*10^-3;
dP = 25;

dy = dP/(S*10^3);

mabs = abs(m);
dd = sqrt((dL*lamb*(mabs-1/2)./y_interference).^2+((mabs-1/2).*L*dy*lamb./(y_interference.^2)).^2);
dd_mean = sqrt(mean(dd.^2));
dd_std = std(dd);
fprintf("Interference\n");
fprintf("a = %.3e +- %.3e\n",d_mean,dd_mean);
fprintf("a_std = %.3e \n",d_std);
fprintf("Error std = %.3e \n\n",dd_std);
