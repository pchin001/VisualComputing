%% Part A: Camera Calibration
clc
clear all

%% Load the previously saved data
load('dscf4177_pixel.mat');
D2 = dscf4177.D2(1:21,:)';
D3 = dscf4177.D3(1:21,:)';

%% Solve using built-in inverses
% C = linsolve(D3',D2')'
C = (D2/D3)

%% Test the C matrix

% choose where to get the test point
% test3D = D3; % original input
% test2D = D2;
test3D = [0 0 0 1]'; % origin
test2D = [222 1563 1]';
% test3D = [20 180 0 1]'; % top left corner
% test2D = [1243 669 1]';
% test3D = [180 200 0 1]'; % top right corner
% test2D = [2455 885 1]';

% calculate and compare
newD2 = C * test3D;
test2D;

% normalize so that w = 1
newD2(1,:) = newD2(1,:) ./ newD2(3,:);
newD2(2,:) = newD2(2,:) ./ newD2(3,:);
newD2(3,:) = newD2(3,:) ./ newD2(3,:);

% find error, by euclidian distance
[~,count] = size(newD2);
diff = abs(test2D - newD2);
error = sqrt(diff(1,:).^2 + diff(2,:).^2);
avg_err = sum(error) ./ count




