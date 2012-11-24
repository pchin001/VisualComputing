function [ pt2 ] = epipole( img1, img2 )
%EPIPOLE finds corresponding point in second image given first image
fMatrix = fundamentalMatrix(img1, img2); 
fMatrix = reshape(fMatrix, 3, 3); 
pt2 = linsolve(fMatrix * pt1, zeros(3, 1)); 
end

