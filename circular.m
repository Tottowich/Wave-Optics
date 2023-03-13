
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

thetas_analytical = deg2rad(-2.5:0.001:2.5);
k = pi./lamb;

%inter = cos(alph).^2;
X = 2*k*(b/2)*sin(thetas_analytical);
R = (2*besselj(1,X)./X).^2;
R = R./max(R);
plot(rad2deg(thetas_analytical),R);
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
offset = 400;
[xAxisDeg, xAxisMm, data3] = image_process(path,channel,rotation,sec_x, ...
               sec_y,baseline_sign,shift,S,L,show_steps,mm_limit,deg_limit,offset);

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
    plot(rad2deg(thetas_analytical),R,"m")
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


