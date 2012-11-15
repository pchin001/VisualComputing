%% Part A: Camera Calibration
clc
clear all

convertXY = 0;
width = 2848;
height = 2136;


%% Load the sample test data 
% load('dscf4177_pixel.mat');
% D2 = dscf4177.D2(1:end,:)';
% D3 = dscf4177.D3(1:end,:)';


%% Load the set 1 data

imgs = [
    'DSCF4177';
    'DSCF4179';
    'DSCF4183';
    'DSCF4186';
    'DSCF4188';
    'DSCF4192';
];

[len,~] = size(imgs);

for i = 1:len
    load_file = strcat('ccdata/', imgs(i,:), '_2D_3D.mat');
    load(load_file);
    D2 = data.D2';
    D3 = data.D3';
    
    %% Convert image coordinates to XY coordinates
    if(convertXY)
        D2(1,:) = D2(1,:) - width/2;
        D2(2,:) = -D2(2,:) + height/2;
    end


    %% Calculate the intrinsic and extrinsic properties
    [C,A,R,T] = CameraCalibration(D2,D3);


    %% Test the C matrix

    % use original input as test points
    test3D = D3; % original input
    test2D = D2;

    % use corners as test points, test2D is from 'DSCF4177.jpg'
%     test3D(:,1) = [0 0 0 1]'; % origin
%     test2D(:,1) = [222 1563 1]';
%     test3D(:,2) = [20 180 0 1]'; % top left corner
%     test2D(:,2) = [1243 669 1]';
%     test3D(:,3) = [180 200 0 1]'; % top right corner
%     test2D(:,3) = [2455 885 1]';

    % calculate
    newD2 = C * test3D;

    % normalize so that w = 1
    newD2(1,:) = newD2(1,:) ./ newD2(3,:);
    newD2(2,:) = newD2(2,:) ./ newD2(3,:);
    newD2(3,:) = newD2(3,:) ./ newD2(3,:);

    % convert (x,y) back to (row,col)
    if(convertXY)
        newD2(1,:) = newD2(1,:) + width/2;
        newD2(2,:) = -(newD2(2,:) - height/2);
    end

    % compare
    % test2D
    % newD2

    % find error, by euclidian distance
    [~,count] = size(newD2);
    avgerr = abs(test2D - newD2); % diff
    avgerr = sqrt(avgerr(1,:).^2 + avgerr(2,:).^2); % euclidian
    avgerr = sum(avgerr) ./ count % average error

    %% print out A to verify
    A

end

