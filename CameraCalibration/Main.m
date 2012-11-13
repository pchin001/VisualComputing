%% Part A: Camera Calibration
clc
clear all

width = 2848;
height = 2136;

%% Load the previously saved data
load('dscf4177_pixel.mat');
D2 = dscf4177.D2(1:6,:)';
D3 = dscf4177.D3(1:6,:)';

%% Load corners of checkboard data and 2 Yellow points
% D3 = [
%       0   0  0 1;
%     180  20  0 1;
%     200 200  0 1;
%       0 200  0 1;
%       0  48 67 1;
%      32  80 67 1;
% ]';
% D2 = [
%      225 1566 1;
%     2028 2097 1;
%     2610  933 1;
%     1185  566 1;
%      373  763 1;
%      820  673 1;
% ]';

%% Convert image coordinates to XY coordinates
D2(1,:) = D2(1,:) - width/2;
D2(2,:) = -D2(2,:) + height/2;

%% Solve using built-in inverses
% C = linsolve(D3',D2')'
% C = (D2/D3)
% C = C ./ C(end,end)

%% Create linear system (matrix)
[~,len] = size(D2);
linsys = zeros(len*2,12);
for i = 1:len
    row = 2*i - 1;
    % u coordinate row
    linsys(row,1:4) = -D3(:,i);
    linsys(row,9:12) = D2(1,i) .* D3(:,i);
    % v coordinate row
    linsys(row+1,5:8) = -D3(:,i);
    linsys(row+1,9:12) = D2(2,i) .* D3(:,i);
end

linsys

%% Solve using rref
ref = rref(linsys);

%% Solve using SVD
[U,S,V] = svd(linsys)
C = zeros(3,4);
C(1,:) = V(1:4,end);
C(2,:) = V(5:8,end);
C(3,:) = V(9:12,end);
% C = C ./ C(end,end);
C

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
% newD2 = C * test3D;
% test2D;

% normalize so that w = 1
% newD2(1,:) = newD2(1,:) ./ newD2(3,:);
% newD2(2,:) = newD2(2,:) ./ newD2(3,:);
% newD2(3,:) = newD2(3,:) ./ newD2(3,:);

% find error, by euclidian distance
% [~,count] = size(newD2);
% diff = abs(test2D - newD2);
% error = sqrt(diff(1,:).^2 + diff(2,:).^2);
% avg_err = sum(error) ./ count


%% Find Position and Orientation of each C-Matrix
A = C(1:3,1:3);

% find R and K using RQ-Decomposition
% this implements RQ Decomposition using QR Decomposition
ReverseRows = [0 0 1; 0 1 0 ; 1 0 0];
[R K] = qr((ReverseRows * A)');
K = ReverseRows * K' * ReverseRows
R = ReverseRows * R'

% find R and K using QR-Decomposition
% [R,K] = qr(A) % R is Rotation Matrix, K is Intrinsic Properies

% normalize K so that last cell is 1
K = K ./ K(end:end)

% find T
T = inv(R) * inv(K) * C(:,end) % KRT = C(last column)
% T = T ./ T(end)

