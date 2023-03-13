function [xAxisDeg, xAxisMm, data3] = image_process(path,channel,rotation,sec_x,sec_y,baseline_sign,shift,S,L,show_steps,mm_limit,deg_limit)
    % BASIC IMAGE ANALYSIS SCRIPT
    % This script is intended as an example to help you analyze the pictures you have acquired
    % during the "Diffraction, interference and Fourier filtering" lab. Notice that you cannot
    % simply run the algorithm "as it is", you will need to change a few things to make it
    % work with your images (therefore you will have to understand it, contact your lab
    % supervisor if you have questions). You can also write your own algorithm making
    % different steps if you want (this here is just for inspiration).
    %% Import the image to Matlab
    arguments
        path (1,1) string
        channel (1,1) double = 3;
        rotation (1,1) double = 0;
        sec_x = [1000,3000];
        sec_y (:,1) double = [1450,1550];
        baseline_sign (1,1) double = -1;
        shift (1,1) double = 0;
        S (1,1) double = 8.4;
        L (1,1) double = 1;
        show_steps (1,1) logical = false;
        mm_limit (:,1) double = [-100,100];
        deg_limit (:,1) double = [-2,2];
    end
    % Print input arguments
    fprintf("path: %s",path)
    fprintf("channel: %d",channel)
    fprintf("rotation: %d",rotation)
    fprintf("sec_x: %s",sec_x)
    fprintf("sec_y: %d",sec_y)
    fprintf("baseline_sign: %d",baseline_sign)
    fprintf("shift: %d",shift)
    fprintf("S: %d",S)
    fprintf("L: %d",L)
    fprintf("show_steps: %d",show_steps)
    fprintf("mm_limit: %d",mm_limit)
    fprintf("deg_limit: %d",deg_limit)
    
    
    myImage1 = imread(path);

    %% Display the image as a figure
    if show_steps
        figure(1)
        imagesc(myImage1);
        title('Step 1: original image');
    end

    %% Select the RGB color channel to use (in this case the blue, as it is less saturated)
    myImage2 = myImage1(:,:,channel);
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
    myImage4 = imrotate(myImage3, rotation);
    if show_steps
        figure(4)
        imagesc(myImage4);
        title('Step 4: after rotation');
    end

    %% Cut a region of interest within the image, where the diffraction pattern is contained
    %myImage5 = myImage4(1450:1550, 1000:3000); Double Slit
    disp(strcmp(sec_x,"all"))
    disp(sec_x)
    if strcmp(sec_x,"all")==1
        myImage5 = myImage4(sec_y, :);
    else
        myImage5 = myImage4(sec_y, sec_x);
    end
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
    data2 = data1 + baseline_sign*baseline;
    if show_steps
        figure(7);
        plot(data2);
    end

    %% Calibrate the pixel with space. Using the mm paper in the image I see that 420 px
    % correspond to 50 mm. Therefore, 8.4 px/mm (I leave the error calculations to you :p ).
    % Notice that I have centered the 0 at the main peak.
    xAxisMm = ((1 : length(data2)) - shift) / S; % Single Slit
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
    if show_steps
        plot(xAxisDeg, data3);
        % xlim([-2 , 2]);
        ylim([0 , 1]);
        xlabel('Angle (°)');
        ylabel('Normalized signal');
        title('Final step: normalized signal vs angle');
        grid on;
    end
end