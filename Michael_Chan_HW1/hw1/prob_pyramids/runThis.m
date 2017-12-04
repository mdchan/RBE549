
function [G, L] = pyramidsGL(im, N)
% [G, L] = pyramidsGL(im, N)
% Creates Gaussian (G) and Laplacian (L) pyramids of level N from image im.
% G and L are cell where G{i}, L{i} stores the i-th level of Gaussian and Laplacian pyramid, respectively. 
% a = 0.375;
% w = [(1/4)-(a/2),1/4,a,1/4,(1/4)-(a/2)];
close all;

image = im2single(imread(im));
cutoff_frequency = 7; 
filter = fspecial('Gaussian', cutoff_frequency*4+1, cutoff_frequency);

% % first level instatiation
G{1} = image;
lowFilterImg = imfilter(image,filter);
Laplacian = image - lowFilterImg;
L{1} = Laplacian;

% % runs up to n iterations
for level = 1:N-1
    
% Gaussian
Gauss = impyramid(G{level}, 'reduce');
G{(level+1)} = Gauss;

%Laplacian;
blurImg = imfilter(G{level},filter);
Laplacian = G{level} - blurImg;
L{(level+1)} = Laplacian;
end
end

function displayPyramids(G, L)

Pyramids = tight_subplot(2,5,[.01 .03],[.1 .01],[.01 .01])
for ii = 1:5; 
    axes(Pyramids(ii)); 
Pyramids(ii) = imshow(G{ii});
end
for ii = 6:10; 
    axes(Pyramids(ii)); 
Pyramids(ii) == imshow(L{(ii-5)});
end

end

function displayFFT(im, minv, maxv)
% Displays FFT images

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE BELOW. Use imfilter() to create 'low_frequencies' and
% 'high_frequencies' and then combine them to create 'hybrid_image'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove the high frequencies from image1 by blurring it. The amount of
% blur that works best will vary with different image pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%low_frequencies = 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove the low frequencies from image2. The easiest way to do this is to
% subtract a blurred version of image2 from the original version of image2.
% This will give you an image centered at zero with negative values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%high_frequencies = 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine the high frequencies and low frequencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%hybrid_image = 

%% Visualize and save outputs
% figure(1); imshow(low_frequencies)
% figure(2); imshow(high_frequencies + 0.5);
% vis = vis_hybrid_image(hybrid_image);
% figure(3); imshow(vis);
% imwrite(low_frequencies, 'low_frequencies.jpg', 'quality', 95);
% imwrite(high_frequencies + 0.5, 'high_frequencies.jpg', 'quality', 95);
% imwrite(hybrid_image, 'hybrid_image.jpg', 'quality', 95);
% imwrite(vis, 'hybrid_image_scales.jpg', 'quality', 95);