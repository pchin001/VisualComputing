%% Depth Estimation
clc
clear all


%% Set necessary information

imgs = [
    '4177';
    '4178';
    '4180';
    '4199';
    '4201';
];

[numCams,~] = size(imgs);
numPts = 2;


%% Load data

data_filename = strcat('../imgset1/CART.mat');
load(data_filename);
Cdata = C(:,:,1:numCams);


%% Create test points

data = zeros(2,numPts,numCams);

% 1st col is [0 0 0], 2nd col is [200, 200, 0]
o4177 = [ 223 1564; 2614 930]';
o4179 = [202 1264; 2703 1026]';
o4813 = [2775 1215; 162 1192]';

% put it all into one variable
data(1:2,1:numPts,1) = o4177;
data(1:2,1:numPts,2) = o4179;
data(1:2,1:numPts,3) = o4813;


%% Estimate the 3D points

results = DepthEstimation(Cdata, data);
results


