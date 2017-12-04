function pyramidsGL(im, N)
% [G, L] = pyramidsGL(im, N)
% Creates Gaussian (G) and Laplacian (L) pyramids of level N from image im.
% G and L are cell where G{i}, L{i} stores the i-th level of Gaussian and Laplacian pyramid, respectively. 

close all;
image = im2single(imread(im));
blurImg = imgaussfilt(image,2);
G{1} = image;
Laplacian = mat2gray(blurImg - G{1});
L{1} = Laplacian;

% runs up to n iterations
for level = 1:N-1
% Gaussian

G{level+1} = impyramid(blurImg,'reduce');

%Laplacian;
blurImg = imgaussfilt(G{level},1);
Laplacian = mat2gray(blurImg - G{level});
L{(level+1)} = Laplacian;
end
   
figure(1)
% This is the right code
Pyramids = tight_subplot(2,5,[.01 .03],[.1 .01],[.01 .01]);
for ii = 1:5; 
    axes(Pyramids(ii)); 
Pyramids(ii) = imshow((G{ii}));
colormap jet
end
for ii = 6:10; 
    axes(Pyramids(ii)); 
Pyramids(ii) = imshow((L{(ii-5)}));
%Pyramids(ii) = imagesc(log(abs(fftshift(fft2(L{ii-5})))),[0,1]);
colormap jet
end



figure(2)
% This is the right code
Pyramids = tight_subplot(2,5,[.01 .03],[.1 .01],[.01 .01]);
for ii = 1:5; 
    axes(Pyramids(ii)); 
% Pyramids(ii) = imshow((G{ii}));
Pyramids(ii) = imagesc(log(abs(fftshift(fft2(rgb2gray((G{ii})))))));
colormap default
end
for ii = 6:10; 
    axes(Pyramids(ii)); 
% Pyramids(ii) = imshow((L{(ii-5)}));
Pyramids(ii) = imagesc(log(abs(fftshift(fft2(rgb2gray(L{ii-5}))))));
colormap default
end


end


