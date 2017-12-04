function [mag,theta] = gradientMagnitude(im,sigma)

% image = im2single(imread(im));
% BW = edge(image,'canny');
% filter = fspecial('laplacian');
% GFilterImg = imgfilter(image,filter);
% GFilterImg = imgaussfilt(image,sigma);
GFilterImg = imgaussfilt(im,sigma);


% red = image(:,:,1); % Red channel\
% green = image(:,:,2); % Green channel
% blue = image(:,:,3); % Blue channel
red = GFilterImg(:,:,1); % Red channel
green = GFilterImg(:,:,2); % Green channel
blue = GFilterImg(:,:,3); % Blue channel

% red = BW(:,:,1); % Red channel
% green = BW(:,:,2); % Green channel
% blue = BW(:,:,3); % Blue channel

% figure, imshow(just_green), title('Red channel')

convMatx = [1,0,-1;
           1,0,-1;
           1,0,-1];
convMaty = [1,1,1;
           0,0,0;
           -1,-1,-1];
%% Red Channel       
Gx_red = imfilter(red, convMatx);
Gy_red = imfilter(red, convMaty);
G_red = sqrt(Gx_red.^2+Gy_red.^2);

% Ix = imfilter(red, double(Gx_red));
% Iy = imfilter(red, double(Gy_red));
% G_red = sqrt(Ix.^2+Iy.^2);

%% Blue Channel
Gx_blue = imfilter(blue, convMatx);
Gy_blue = imfilter(blue, convMaty);
G_blue = sqrt(Gx_blue.^2+Gy_blue.^2);

% Ix = imfilter(blue, double(Gx_blue));
% Iy = imfilter(blue,double(Gy_blue));
% G_blue = sqrt(Ix.^2+Iy.^2);
%% Green Channel
Gx_green = imfilter(green, convMatx);
Gy_green = imfilter(green, convMaty);
G_green = sqrt(Gx_green.^2+Gy_green.^2);

% Ix = imfilter(green, double(Gx_green));
% Iy = imfilter(green,double(Gy_green));
% G_green = sqrt(Ix.^2+Iy.^2);
%% Find the Max pixel
G_size = zeros(size(G_red));
theta = zeros(size(G_red));
mag = zeros(size(G_red));
for i = 1:numel(G_size)
     pixelGrad = [G_red(i),G_blue(i),G_green(i)];
     mag(i) = max(pixelGrad);
     if mag(i) == G_red(i)
         theta(i) = atan2(Gy_red(i),Gx_red(i));
     end
     
     if mag(i) == G_blue(i)
         theta(i) = atan2(Gy_blue(i),Gx_blue(i));
     end
     
     if mag(i) == G_green(i)
         theta(i) = atan2(Gy_green(i),Gx_green(i));
     end
end

end