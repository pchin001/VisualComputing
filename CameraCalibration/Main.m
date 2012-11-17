%% Part A: Camera Calibration

% http://www.vision.caltech.edu/bouguetj/calib_doc/papers/heikkila97.pdf
% http://excelsior.cs.ucsb.edu/courses/cs290i_mvg/pdf/calibration.pdf

%% setup
clc
clear all

% set to 1 to see test points on image plot
errorCheck = 0;

% set to 1 to use XY coordinate system
% set to 0 to use image coordinate system
convertXY = 1;
imgw = 2848;
imgh = 2136;


%% Load the image set to use

% set 1
imgs = [
    'DSCF4177';
    'DSCF4179';
    'DSCF4183';
    'DSCF4186';
    'DSCF4188';
    'DSCF4192';
];

% set 2
% imgs = [
% ];

%% Initialize variables to save calculations from each image
[len,~] = size(imgs);
projection = reshape(repmat(eye(4,4),1,len), 4,4,len);
intrinsic = reshape(repmat(eye(4,4),1,len), 4,4,len);
rotation = reshape(repmat(eye(4,4),1,len), 4,4,len);
translation = reshape(repmat(eye(4,4),1,len), 4,4,len);
avgerr = zeros(1,len);


%% Begin calculation for each image
for i = 1:len
    load_file = strcat('ccdata/', imgs(i,:), '_2D_3D.mat');
    load(load_file);
    
    % use odd set of data to calculate C
    data2D = data.D2(1:end,:)';
    data3D = data.D3(1:end,:)';
    
    % choose some 3D test points
    test3D = [
          0,   0,  50, 1;
          0,   0,   0, 1;
        100, 100,   0, 1;
        100, 100,  50, 1;
    ]';
    
    % Convert image coordinates to XY coordinates
    if(convertXY)
        data2D(1,:) = data2D(1,:) - imgw/2;
        data2D(2,:) = -data2D(2,:) + imgh/2;
    end

    % Calculate the intrinsic and extrinsic properties
%     [C, A, R, T, ~, ~, ~, ~] = funcCalibrate(data3D(1:3,:)',data2D(1:2,:)');
    [C,A,R,T] = CameraCalibration(data2D,data3D);
    R = sign(C(3,4)) .* R;
%     T = sign(C(3,4)) .* T;
    
    % store this image's A to compare with other images' A
    projection(1:3,1:4,i) = C;
    intrinsic(1:3,1:3,i) = A;
    rotation(1:3,1:3,i) = R;
    translation(1:3,4,i) = T;
    
    % calculate prediction, and compare with image
    if(errorCheck)
        pred2D = C * test3D;
        pred2D = pred2D ./ repmat(pred2D(3,:),3,1); % normalize w = 1    
        imshow(strcat('../imgset1/', imgs(i,:), '.jpg'));
        hold on;
        plot(pred2D(1,:),pred2D(2,:),'Marker','p','Color',[.88 .48 0],'MarkerSize',20)
        pause;
    end
    
end



%% print out the calculations from each image

% mean and variance of A
intrinsic;
intrinsicMean = mean(intrinsic,3);
intrinsicStd = std(intrinsic,0,3);
intrinsicScale = intrinsicStd ./ intrinsicMean;

    
%% save the information
C = projection(1:3,1:4,:);
A = intrinsic(1:3,1:3,:)
R = rotation(1:3,1:3,:)
T = translation(1:3,4,:)
save_file = strcat('../imgset1/CART.mat');
save(save_file, 'C', 'A', 'R', 'T');


%% plot orientation and position of cameras

PlotCamerasInScene(rotation, translation);

