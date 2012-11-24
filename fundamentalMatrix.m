function [ fmatrix ] = fundamentalMatrix( img1, img2 )
%FUNDAMENTALMATRIX Calculates the fundamental matrix of the given image set
%   Detailed explanation goes here

%% compute fundamental matrix 
% x*x2*f1 + y*x2*f2 + x2*f3 + x*y2*f4 + y*y2*f5 + y2*f6 + x*f7 + y*f8 + 1*f9 = 0 
% [ x*x2 y*x2 x2 x*y2 y*y2 y2 x y ] [f1] = -1
pts1 = ginput(8); 
pts2 = ginput(8); 
A = zeros(size(pts1, 1), 8);
    for i = 1:size(pts1, 1)
        A(i, :) = [pts1(i,1)*pts2(i,1), pts1(i, 2)*pts2(i, 1), ...
                   pts2(i,1), pts1(i, 1)*pts2(i, 2), pts1(i, 2)*pts2(i, 2), ...
                   pts2(i, 2), pts1(i, 1), pts1(i, 2)]; 
    end
    f = linsolve(A, zeros(8, 1)); 
end

