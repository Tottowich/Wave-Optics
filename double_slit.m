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
b = 4*10^-5;%;4*10^-5;
a = 1.25*10^-4;%1.25*10^-4;
L = 0.972;
N = 2;
new_b_min = @(y,m) lamb*m*L./y;
new_a_min = @(y,m) lamb*(m-1/2)*L./y;
thetas_analytical = deg2rad(-2.5:0.001:2.5);
Beta = @(theta) pi*b*sin(theta)./lamb;
Alpha = @(theta) pi*a*sin(theta)./lamb;
beta = Beta(thetas_analytical);
alph = Alpha(thetas_analytical);
diff = (sin(beta)./beta).^2;
%inter = cos(alph).^2;
inter = (sin(N.*alph)./sin(alph)).^2;%.*(mod(alph,pi)~=0)+N.^2*(mod(alph,pi)==0);
R = diff.*inter;
R = R./max(R);
plot(rad2deg(thetas_analytical),R);
rel_err = @(y,m) (new_b(y,m)-b)/b;
%% Import the image to Matlab
path = 'Images/DoubleSlitNew/doubleSlit.jpg';
channel = 1;
rotation = 3.6;
sec_x = 1000:3000;
sec_y = 1500:1550;
baseline_sign = 0; % No baseline subtraction or addition
shift = (782+1412)/2;
S = (90/5);
L = 0.972;
show_steps = true;
mm_limit = [-100,100];
deg_limit = [-2,2];
offset = 1700;
[xAxisDeg, xAxisMm, data3] = image_process(path,channel,rotation,sec_x, ...
               sec_y,baseline_sign,shift,S,L,show_steps,mm_limit,deg_limit,offset);
%% Interference : a
visualize = true;
[PKS,LOCS] = find_min(data3,xAxisMm,0.5,2.5,visualize);
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
    plot(rad2deg(thetas_analytical),R,"m");
    hold off
    grid on;
    xlabel('Angle (°)');
    ylabel('Normalized signal');
    title('Normalized signal vs angle with Min peaks');
    xlim([-3.2, 3.2]);
    ylim([0 , 1]);
    legend(["Exp - values","Minima","Analytical"])  
end
left_hand = sum(y_interference<0);
right_hand = sum(y_interference>=0);
m_inter = [-left_hand:-1 1:right_hand];
slit_spacings = new_a_min(y_interference,m_inter);
a_mean = mean(slit_spacings);
a_std = std(slit_spacings);


%% Interference Error:
N_inter = length(slit_spacings);
dL = 5*10^-3;
dP = 25;

dy = dP/(S*10^3);

mabs_inter = abs(m_inter);
da = sqrt((dL*lamb*(mabs_inter-1/2)./y_interference).^2+((mabs_inter-1/2).*L*dy*lamb./(y_interference.^2)).^2);
da_mean = mean(da.^2);
da_std = std(da);

fprintf("Interference\n");
fprintf("a = %.3e +- %.3e\n",a_mean,da_mean);
fprintf("a_std = %.3e \n",a_std);
fprintf("Error std = %.3e \n\n",da_std);
%%
% Diffraction : b

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
slit_widths = new_b_min(y_diffraction,m_diff);
b_mean = mean(slit_widths);
b_std = std(slit_widths);
%% Diffraction Error:
N_diff = length(slit_widths);
mabs_diff = abs(m_diff);
db =sqrt((dL*lamb*mabs_diff./y_diffraction).^2+(mabs_diff.*L*dy*lamb./(y_diffraction.^2)).^2);
db_mean = sqrt(mean(db.^2));
db_std = std(db);
fprintf("Diffraction\n");
fprintf("b = %.3e +- %.3e\n",b_mean,db_mean);
fprintf("b_std = %.3e \n",b_std);
fprintf("Error std = %.3e \n",db_std);
