function [ output_args ] = DepthEstimation( C, data )
%DEPTHESTIMATION Summary of this function goes here
%   C1      Camera Calibration Matrices (3x4xN)
%   data    Corresponding points between images (2x1xN)

%% Check for errors
[~,~,numCameras] = size(C);
[~,~,numData] = size(data);

if numCameras ~= numData
    error('C and data need to have the same 3rd dimension.');
end

if(numCameras < 2)
    error('Need at least 2 cameras to calculate.');
end


%% Set up a linear system of equations for solving





end

