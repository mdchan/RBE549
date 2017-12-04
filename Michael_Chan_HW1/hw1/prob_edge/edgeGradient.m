function bmap = edgeGradient(im)
c = edge(rgb2gray(im),'canny');
[m,t] = gradientMagnitude(im,2);
size(c)
size(m)
bmap = (c .* m);
end