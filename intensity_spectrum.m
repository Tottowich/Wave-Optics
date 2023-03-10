img = imread("Images/DoubleSlit/IMG-0729.jpg");
imshow(img)
I = improfile();
%%
red = I(:,:,3);
plot(red)
%red = I(:,:,1);