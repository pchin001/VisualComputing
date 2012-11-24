function [ results ] = DepthEstimation( C, data )
%DEPTHESTIMATION Summary of this function goes here
%   C1      Camera Calibration Matrices (3x4xN)
%   data    Corresponding points between images (2xMxN)


%% Check for errors

[~,~,numCameras] = size(C);
[~,numPts,numData] = size(data);

if numCameras ~= numData
    error('C and data need to have the same 3rd dimension.');
end

if(numCameras < 2)
    error('Need at least 2 cameras to calculate.');
end


%% Calculate the predictions for all 3D points

results = zeros(4,numPts);

for p = 1 : numPts
    %% Set up a linear system of equations for solving
    linsys = zeros(numCameras,4);
    for cam = 1 : numCameras
        row = 2 * cam - 1;
        % u-coordinate row
        linsys(row,1) = data(1,p,cam)*C(3,1,cam) - C(1,1,cam);
        linsys(row,2) = data(1,p,cam)*C(3,2,cam) - C(1,2,cam);
        linsys(row,3) = data(1,p,cam)*C(3,3,cam) - C(1,3,cam);
        linsys(row,4) = data(1,p,cam)*C(3,4,cam) - C(1,4,cam);
        % v-coordinate row
        linsys(row+1,1) = data(2,p,cam)*C(3,1,cam) - C(2,1,cam);
        linsys(row+1,2) = data(2,p,cam)*C(3,2,cam) - C(2,2,cam);
        linsys(row+1,3) = data(2,p,cam)*C(3,3,cam) - C(2,3,cam);
        linsys(row+1,4) = data(2,p,cam)*C(3,4,cam) - C(2,4,cam);
    end

    %% Solve the linear system with SVD
    [U,S,V] = svd(linsys);
    P = V(:,end);
    results(:,p) = P ./ P(end);

end



end

