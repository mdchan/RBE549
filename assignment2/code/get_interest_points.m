% Local Feature Stencil Code
% Written by James Hays for CS 4476/6476 @ Georgia Tech

% Returns a set of interest points for the input image

% 'image' can be grayscale or color, your choice.
% 'feature_width', in pixels, is the local feature width. It might be
%   useful in this function in order to (a) suppress boundary interest
%   points (where a feature wouldn't fit entirely in the image, anyway)
%   or (b) scale the image filters being used. Or you can ignore it.

% 'x' and 'y' are nx1 vectors of x and y coordinates of interest points.
% 'confidence' is an nx1 vector indicating the strength of the interest
%   point. You might use this later or not.
% 'scale' and 'orientation' are nx1 vectors indicating the scale and
%   orientation of each interest point. These are OPTIONAL. By default you
%   do not need to make scale and orientation invariant local features.
function [x,y] = get_interest_points(image, feature_width)

% Implement the Harris corner detector (See Szeliski 4.1.1) to start with.
% You can create additional interest point detector functions (e.g. MSER)
% for extra credit.

% If you're finding spurious interest point detections near the boundaries,
% it is safe to simply suppress the gradients / corners near the edges of
% the image.

% The lecture slides and textbook are a bit vague on how to do the
% non-maximum suppression once you've thresholded the cornerness score.
% You are free to experiment. Here are some helpful functions:
%  BWLABEL and the newer BWCONNCOMP will find connected components in 
% thresholded binary image. You could, for instance, take the maximum value
% within each component.
%  COLFILT can be used to run a max() operator on each sliding window. You
% could use this to ensure that every interest point is at a local maximum
% of cornerness.

% find the size of the image
[imrow,imcol] = size(image);
threshold = .25;
siz = [imrow,imcol];
sig = 1;

% derivative matrices
convMatx = [1,0,-1;
           1,0,-1;
           1,0,-1];
convMaty = [1,1,1;
           0,0,0;
           -1,-1,-1];
   
% find the x and y derivatives
Lx = imfilter(image, convMatx,'same');
Ly = imfilter(image, convMaty, 'same');

% put together the different parts of the equation
Lx2 = Lx .* Lx;
Ly2 = Ly .* Ly;
Lxy = Lx .* Ly;

% harris filter
H = [-2 -1 0 1 2];

% gaussian filter
% H = fspecial('gaussian',max(1,fix(6*sig)), sig);

% apply harris filter
Gx = imfilter(Lx2,H,'same');
Gy = imfilter(Ly2,H,'same');
Gxy = imfilter(Lxy,H,'same');

har = Gx .* Gy - Gxy.^2 - 0.005*(Gx + Gy).^2;
% har = (Gx.*Gy - Gxy.^2) ./ (Gx + Gy+eps);


radius = 1;
order = 2*radius +1;
mx = ordfilt2(har,order^2,ones(order));
harris_points = (har == mx)&(har >= threshold);
[y,x] = find(harris_points);

% this take the harris points and trims them so that we can take the 
% 16x16 window around them
badx = find(x<=20|x>=imrow-20);
x(badx) = [];
y(badx) = [];
bady = find(y<=20|y>=imcol-20);
y(bady) = [];
x(bady) = [];

%% this is the code for plotting harris points
% figure
% imshow(image),hold on
% plot(x,y,'ys'),

end

