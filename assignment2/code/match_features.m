% Local Feature Stencil Code
% Written by James Hays for CS 4476/6476 @ Georgia Tech

% 'features1' and 'features2' are the n x feature dimensionality features
%   from the two images.
% If you want to include geometric verification in this stage, you can add
% the x and y locations of the features as additional inputs.
%
% 'matches' is a k x 2 matrix, where k is the number of matches. The first
%   column is an index in features1, the second column is an index
%   in features2.
% 'Confidences' is a k x 1 matrix with a real valued confidence for every
%   match.
% 'matches' and 'confidences' can empty, e.g. 0x2 and 0x1.
function [matches, confidences] = match_features(features1, features2)

% This function does not need to be symmetric (e.g. it can produce
% different numbers of matches depending on the order of the arguments).

% To start with, simply implement the "ratio test", equation 4.18 in
% section 4.1.3 of Szeliski. For extra credit you can implement various
% forms of spatial verification of matches.
% features1 = mat2cell(features1,1,128);
% features2 = mat2cell(features2,1,128);

a = .1;
matches = [];
confidences = [];

for i=1:size(features1)%k
    for j = 1:size(features2)%k
        %         to find the smallest distance find difference between position i
        %         against every element in feature 2
        %         SSD method
%                 d(j) = sum((features1(i,:)-features2(j,:)).^2);
        
        %         norm method
        d(j) = sum(abs(features1(i,:)-features2(j,:)));
        %         now we have a 1x128 vector with all the differences of i against
        %         features2
    end
%     sort it, 
    [B,I] = sort(d);
    d1 = B(1);
    d2 = B(2);
    ratio = d1/d2;
    if ratio < .7
        conf = exp(-a*d1);
        confidences = vertcat(confidences,conf);
        idx = [i,I(1)];
        matches = vertcat(matches,idx);
        
    end
    
end



% % Placeholder that you can delete. Random matches and confidences
% num_features = min(size(features1, 1), size(features2,1));
% matches = zeros(num_features, 2);
% matches(:,1) = randperm(num_features);
% matches(:,2) = randperm(num_features);
% confidences = rand(num_features,1);

% Sort the matches so that the most confident onces are at the top of the
% list. You should probably not delete this, so that the evaluation
% functions can be run on the top matches easily.
[confidences, ind] = sort(confidences, 'descend');
matches = matches(ind,:);