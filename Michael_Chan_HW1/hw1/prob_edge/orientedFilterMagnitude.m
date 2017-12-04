function [mag,theta] = orientedFilterMagnitude(im)
% 
% image = im2single(imread(im));
% BW = edge(image,'canny');
% GFilterImg = imgaussfilt(image);
GFilterImg = imgaussfilt(im);
orientation = linspace(0,2*pi,20);

red = GFilterImg(:,:,1); % Red channel
green = GFilterImg(:,:,2); % Green channel
blue = GFilterImg(:,:,3); % Blue channel


convMatx = [1,0,-1;
           1,0,-1;
           1,0,-1];
convMaty = [1,1,1;
           0,0,0;
           -1,-1,-1];
       
% rotMat = [1 0 0; 0 cos(phi) -sin(phi);0 sin(phi) cos(phi)];
for element = 1:numel(orientation)
%% Red Channel    
rotMat = [1 0 0; 0 cos(orientation(element)) -sin(orientation(element));
    0 sin(orientation(element)) cos(orientation(element))];

Gx_red = imfilter(red, convMatx);
Gy_red = imfilter(red, convMaty);
Gx_red_rot = imfilter(Gx_red,rotMat);
Gy_red_rot = imfilter(Gy_red,rotMat);
G_red_rot = sqrt(Gx_red_rot.^2+Gy_red_rot.^2);

%% Blue Channel
Gx_blue = imfilter(blue, convMatx);
Gy_blue = imfilter(blue, convMaty);
Gx_blue_rot = imfilter(Gx_blue,rotMat);
Gy_blue_rot = imfilter(Gy_blue,rotMat);
G_blue_rot = sqrt(Gx_blue_rot.^2+Gy_blue_rot.^2);

%% Green Channel
Gx_green = imfilter(green, convMatx);
Gy_green = imfilter(green, convMaty);
Gx_green_rot = imfilter(Gx_green,rotMat);
Gy_green_rot = imfilter(Gy_green,rotMat);
G_green_rot = sqrt(Gx_green_rot.^2+Gy_green_rot.^2);


%% Find the Max pixel
G_size = zeros(size(Gx_green));
theta = zeros(size(Gx_green));
mag = zeros(size(Gx_green));
for i = 1:numel(G_size)
pixelGrad = [G_red_rot(i),G_blue_rot(i),G_green_rot(i)];
mag(i) = max(pixelGrad);
end
for i = 1:numel(G_size)
     pixelGrad = [G_red_rot(i),G_blue_rot(i),G_green_rot(i)];
     mag(i) = max(pixelGrad);
     if mag(i) == G_red_rot(i)
         theta(i) = cos(orientation(element)).*Gy_red(i)+sin(orientation(element)).*Gx_red(i);
     end
     
     if mag(i) == G_blue_rot(i)
         theta(i) = cos(orientation(element)).*Gy_blue(i)+sin(orientation(element)).*Gx_blue(i);
     end
     
     if mag(i) == G_green_rot(i)
         theta(i) = cos(orientation(element)).*Gy_green(i)+sin(orientation(element)).*Gx_green(i);
     end
end



% maxRB = max(G_red_rot,G_blue_rot);
% mag = max(maxRB,G_green_rot);
% 
% end

% figure(1)
% imshow(mag);
% 
% figure(2)
% imagesc(theta,[0,2*pi]);
% nonmax(image,theta)

end