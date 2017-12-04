% % Local Feature Stencil Code
% % Written by James Hays for CS 4476/6476 @ Georgia Tech
%
% % Returns a set of feature descriptors for a given set of interest points.
%
% % 'image' can be grayscale or color, your choice.
% % 'x' and 'y' are nx1 vectors of x and y coordinates of interest points.
% %   The local features should be centered at x and y.
% % 'feature_width', in pixels, is the local feature width. You can assume
% %   that feature_width will be a multiple of 4 (i.e. every cell of your
% %   local SIFT-like feature will have an integer width and height).
% % If you want to detect and describe features at multiple scales or
% % particular orientations you can add input arguments.
%
% % 'features' is the array of computed features. It should have the
% %   following size: [length(x) x feature dimensionality] (e.g. 128 for
% %   standard SIFT)
%
function [features] = get_features(image, x, y, feature_width)
%
% % To start with, you might want to simply use normalized patches as your
% % local feature. This is very simple to code and works OK. However, to get
% % full credit you will need to implement the more effective SIFT descriptor
% % (See Szeliski 4.1.2 or the original publications at
% % http://www.cs.ubc.ca/~lowe/keypoints/)
%
% % Your implementation does not need to exactly match the SIFT reference.
% % Here are the key properties your (baseline) descriptor should have:
% %  (1) a 4x4 grid of cells, each feature_width/4. 'cell' in this context
% %    nothing to do with the Matlab data structue of cell(). It is simply
% %    the terminology used in the feature literature to describe the spatial
% %    bins where gradient distributions will be described.
% %  (2) each cell should have a histogram of the local distribution of
% %    gradients in 8 orientations. Appending these histograms together will
% %    give you 4x4 x 8 = 128 dimensions.
% %  (3) Each feature should be normalized to unit length
% %
% % You do not need to perform the interpolation in which each gradient
% % measurement contributes to multiple orientation bins in multiple cells
% % As described in Szeliski, a single gradient measurement creates a
% % weighted contribution to the 4 nearest cells and the 2 nearest
% % orientation bins within each cell, for 8 total contributions. This type
% % of interpolation probably will help, though.
%
% % You do not have to explicitly compute the gradient orientation at each
% % pixel (although you are free to do so). You can instead filter with
% % oriented filters (e.g. a filter that responds to edges with a specific
% % orientation). All of your SIFT-like feature can be constructed entirely
% % from filtering fairly quickly in this way.
%
% % You do not need to do the normalize -> threshold -> normalize again
% % operation as detailed in Szeliski and the SIFT paper. It can help, though.
%
% % Another simple trick which can help is to raise each element of the final
% % feature vector to some power that is less than one.
%
%
%
% [height,width] = size(image);
% filter = fspecial('gaussian',[3 3],2);
% imagef = imfilter(image,filter);
% % filter the image
% [Gmag,Gdir] = imgradient(imagef);
% % find mag and grad of gradient
% binang = [45,90,100,135,180,-135,-90,-45,0];
% % create bins
% for i = 1:size(x)
%     current_x = x(i,:);
%     current_y = y(i,:);
%     magG = Gmag(current_x-8:current_x+7,current_y-8:current_y+7);
%     dirG = Gdir(current_x-8:current_x+7,current_y-8:current_y+7);
% %     make 16x16 windows around feature point (mag and theta)
% %     feature = [45 90 100 135 180 -135 -90 -45 0];
% feature = [0 45 90 135 180 -135 -90 -45];
%     for j = 1 : 4 : size(magG,1)
%         for k = 1 : 4 : size(dirG,2)
%             mag = magG(j : j+3,k : k+3);
%             theta = dirG(j : j+3, k : k+3);
% %             makes 4x4 windows
%             Orientation = zeros(1,8);
% %          Bins the magnitudes based on orientation
%             for m =1:size(mag,1)
%                 for n=1:size(mag,2)
%                     if (theta(m,n)>=binang(1)) && (theta(m,n)<binang(2))
%                         Orientation(1) = Orientation(1) + mag(m,n);
%                     else
%                         if (theta(m,n)>=binang(3)) && (theta(m,n)<binang(4))
%                             Orientation(2) = Orientation(2)+ mag(m,n);
%                         else
%                             if (theta(m,n)>=binang(4)) && (theta(m,n)<binang(5))
%                                 Orientation(3) = Orientation(3)+ mag(m,n);
%                             else
%                                 if (theta(m,n)>=binang(5)) && (theta(m,n)<binang(6))
%                                     Orientation(4) = Orientation(4)+ mag(m,n);
%                                 else
%                                     if (theta(m,n)>=binang(6)) && (theta(m,n)<binang(7))
%                                         Orientation(5) = Orientation(5)+ mag(m,n);
%                                     else
%                                         if (theta(m,n)>=binang(7)) && (theta(m,n)<binang(8))
%                                             Orientation(6) = Orientation(6)+ mag(m,n);
%                                         else
%                                             if (theta(m,n)>=binang(8)) && (theta(m,n)<binang(9))
%                                                 Orientation(7) = Orientation(7)+ mag(m,n);
%                                             else
%                                                 Orientation(8) = Orientation(8)+ mag(m,n);
%                                             end
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%             feature = [feature Orientation];
% %             append feature to the feature matrix
%         end
%     end
%     features(i,:) = feature./norm(feature);
% %     normalize the features
% end
% features = features .^ .75;
%
% end
%
%
%
xpos = round(x);
ypos = round(y);
convMatx = [1,0,-1;
    1,0,-1;
    1,0,-1];
convMaty = [1,1,1;
    0,0,0;
    -1,-1,-1];
features = [];
[height,width] = size(image);
filter = fspecial('gaussian',[3 3],2);
image = imfilter(image,filter);
[Gmag,Gdir] = imgradient(image);
%16x16 patch around features
for k = 1:numel(xpos)
    
    %     calculate the 16x16 pixel image patch
    patchMag = Gmag(xpos(k)-8:xpos(k)+7,ypos(k)-8:ypos(k)+7);
    patchDir = Gdir(xpos(k)-8:xpos(k)+7,ypos(k)-8:ypos(k)+7);
    
    %     allocate bins
    binang = linspace(0,2*pi,9);
    %     for every pixel in the patch, calculate gradient orientation
    %     sum it up for the gradient orientation of image
    for i=1:numel(patchMag)
        %   mag
        mag = sqrt(gradx(i).^2+grady(i).^2);
        %   theta
        theta = abs((atan2(gradx(i),grady(i))));
        
        % bin the mags based on theta
        if (theta>binang(1)) && (theta<=binang(2))
            Orientation(1) = Orientation(1) + mag;
        end
        if (theta>binang(2)) && (theta<=binang(3))
            Orientation(2) = Orientation(2)+ mag;
        end
        if (theta>binang(3)) && (theta<=binang(4))
            Orientation(3) = Orientation(3)+ mag;
        end
        if (theta>binang(4)) && (theta<=binang(5))
            Orientation(4) = Orientation(4)+ mag;
        end
        if (theta>binang(5)) && (theta<=binang(6))
            Orientation(5) = Orientation(5)+ mag;
        end
        if (theta>binang(6)) && (theta<=binang(7))
            Orientation(6) = Orientation(6)+ mag;
        end
        if (theta>binang(7)) && (theta<=binang(8))
            Orientation(7) = Orientation(7)+ mag;
        end
        if (theta>binang(8)) && (theta<=binang(9))
            Orientation(8) = Orientation(8)+ mag;
        end
    end
end

%     now we have gradient orientation of the whole 16x16 patch (8x1
%     vector)

%     find the index of the max element in orientation
[~,I] = max(patchOrientation);
%     take a larger image around the feature point
bigPatch = image(xpos(k)-10:xpos(k)+10,ypos(k)-10:ypos(k)+10);
%     rotate the image patch by the largest mag, setting it to ground truth
%     cut down to be 16x16 again
normPatch = imrotate(bigPatch,binang(I));
normPatch(17:end,:) = [];
normPatch(:,17:end) = [];
% now we have normalized 16x16 patch


%     blur the patch through a gaussian
normPatchG = imgaussfilt(normPatch,.5);

%     find x and y derivatives
normgradx = imfilter(normPatchG, convMatx);
normgrady = imfilter(normPatchG, convMaty);

%     take 16x16 normPatch, and turn it into 16 4x4 matrices
%     break patch into 16 4x4 cells
N = 4*ones(1,4);
pixSq = mat2cell(normPatchG, N,N);
normgradxSplit = mat2cell(normgradx, N,N);
normgradySplit = mat2cell(normgradx, N,N);
pixFeature = [];
%    compute gradient for each pixel in normpatchG
for i=1:numel(pixSq)
    Orientation = zeros(1,8);
    for j = 1:numel(pixSq{i})
        % magnitude
        mag = sqrt(normgradxSplit{i}(j).^2+normgradySplit{i}(j).^2);
        % orientation
        theta = abs((atan2(normgradxSplit{i}(j),normgradySplit{i}(j))));
        % bin the mags based on theta
        if (theta>binang(1)) && (theta<=binang(2))
            Orientation(1) = Orientation(1) + mag;
        end
        if (theta>binang(2)) && (theta<=binang(3))
            Orientation(2) = Orientation(2)+ mag;
        end
        if (theta>binang(3)) && (theta<=binang(4))
            Orientation(3) = Orientation(3)+ mag;
        end
        if (theta>binang(4)) && (theta<=binang(5))
            Orientation(4) = Orientation(4)+ mag;
        end
        if (theta>binang(5)) && (theta<=binang(6))
            Orientation(5) = Orientation(5)+ mag;
        end
        if (theta>binang(6)) && (theta<=binang(7))
            Orientation(6) = Orientation(6)+ mag;
        end
        if (theta>binang(7)) && (theta<=binang(8))
            Orientation(7) = Orientation(7)+ mag;
        end
        if (theta>binang(8)) && (theta<=binang(9))
            Orientation(8) = Orientation(8)+ mag;
        end
    end
    %         this is a 1x128 matrix of the feature descriptor
    pixFeature = horzcat(pixFeature,Orientation);
end
% add the feature descriptor to the features matrix
% should be k x 128
features = vertcat(features,pixFeature);


end

% this normalizes the final feature vector?
features = features.^.75;

end