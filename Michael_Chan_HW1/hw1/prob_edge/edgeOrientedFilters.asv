function bmap = edgeOrientedFilters(im)
c = edge(rgb2gray(im),'canny');
[m,t] = orientedFilterMagnitude(im);
size(c)
size(m)
bmap = (c .* m);
end