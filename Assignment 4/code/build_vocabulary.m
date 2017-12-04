% Starter code prepared by James Hays for Computer Vision

% This function will extract a set of feature descriptors from the training images,
% cluster them into a visual vocabulary with k-means,
% and then return the cluster centers.

% Notes:
% - To save computation time, we might consider sampling from the set of training images.
% - Per image, we could randomly sample descriptors, or densely sample descriptors,
% or even try extracting descriptors at interest points.
% - For dense sampling, we can set a stride or step side, e.g., extract a feature every 20 pixels.
% - Recommended first feature descriptor to try: HOG.

% Function inputs: 
% - 'image_paths': a N x 1 cell array of image paths.
% - 'vocab_size' the size of the vocabulary.

% Function outputs:
% - 'vocab' should be vocab_size x descriptor length. Each row is a cluster centroid / visual word.

function vocab = build_vocabulary( image_paths, vocab_size )

[imNum,~] = size(image_paths);
feats = [];

for i=1:length(image_paths)
    im = imread(image_paths{i});
    
    surfFeat = detectSURFFeatures(im);
%     [features,validPoints] = extractFeatures(im,surfFeat.selectStrongest(300));
    points = surfFeat.selectStrongest(300);
    [features,validPoints] = extractHOGFeatures(im,points);
    noSamp = size(features,1);
    randomPerm = randperm(size(features,1));
    sampPerc = round(noSamp/5);
    sample = randomPerm(1:sampPerc);
    extractFeat = features(sample,:);
    feats = vertcat(feats,extractFeat);
    
end
[IDX,vocab] = kmeans((feats),vocab_size);
end



