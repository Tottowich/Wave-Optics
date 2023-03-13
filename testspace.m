%%
% Example image:
RGB = imread("./Images/UnFiltered/sharp.jpg");
% Physically filtered image:
RGB_filtered = imread("./Images/Filtered/cropped.png");
cols = 3;
rows = 2;
subplot(rows,cols,1);
radius = 25;
imshow(RGB)
% Gray scale:
title("Input Image")
fourier_plane = fftshift(fft2(RGB));
fabs = abs(fourier_plane);
fourier_img = fabs*255/max(fabs,[],'all');


[N,M] = size(RGB(:,:,1));
[X,Y] = meshgrid(1:N,1:M);

subplot(rows,cols,4);
imshow(centercrop(fourier_img,N,M))
title("Fourier Plane")

centerX = ceil(N/2);
centerY = ceil(M/2);
distFromCenter =sqrt((X-centerX).^2+(Y-centerY).^2);
lowpass = (distFromCenter <= radius)';
fourier_filtered = fourier_plane.*lowpass;
fourier_removed = fourier_plane.*(lowpass==0);

fourier_img_filtered = fourier_img.*lowpass;
fourier_img_removed = fourier_img.*(lowpass==0);


subplot(rows,cols,5);
imshow(centercrop(fourier_img_removed,N,M))
title("Fourier Plane Noise")
subplot(rows,cols,6);
imshow(centercrop(fourier_img_filtered,N,M))
title("Fourier Plane Filtered")

filtered_img = real(ifft2(ifftshift(fourier_filtered)))/255;
removed_img = real(ifft2(ifftshift(fourier_removed)))/255;
subplot(rows,cols,2)
imshow(removed_img)
title("Removed Noise")

ax = subplot(rows,cols,3);
imshow(filtered_img)
sharp = imsharpen(filtered_img,'Amount',1);
subplot(rows,cols,3)
imshow(sharp)
title("Filtered Image")
sgtitle("Fourier Transform")
ax.TitleHorizontalAlignment = 'left'; 
%%
% Subplot of the cropped pre filtered image and the cropped post filtered image:
figure(1)
subplot(1,2,1)
imshow(RGB_filtered)
title("Physically Filtered")
subplot(1,2,2)
imshow(sharp)
title("Post Filtered")
sgtitle("Fourier Transform")
% Subplot of the cropped pre filtered image and input image:
figure(2)
subplot(1,2,1)
imshow(RGB)
title("Input Image")
subplot(1,2,2)
imshow(RGB_filtered)
title("Physically Filtered")
sgtitle("Fourier Transform")


function I = centercrop(I,N,M)

    x_cent = ceil(N/2);
    y_cent = ceil(M/2);
    size_of_cropped_img = 200;
    %I2 = imcrop(I,rect) crops the image I. rect is a four-element position vector of the
    %form [xmin ymin width height] that specifies the size and position of the crop rectangle. 
    %imcrop returns the cropped image, I2.
    xmin = x_cent-size_of_cropped_img/2;
    ymin = y_cent-size_of_cropped_img/2;
    %I = imcrop(I,[xmin+1 ymin+1 size_of_cropped_img-1 size_of_cropped_img-1]);
    I = I(xmin+1:xmin+size_of_cropped_img,ymin+1:ymin+size_of_cropped_img,:);
end