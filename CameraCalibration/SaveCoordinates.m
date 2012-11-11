%% SaveCoordinates.m
clc
clear all

%% Set 1
imgs = [
    'DSCF4177';
    'DSCF4179';
    'DSCF4183';
    'DSCF4186';
    'DSCF4188';
    'DSCF4192';
];

%% Set 2


%% Collect coordinates from images

n = 24;
[len,~] = size(imgs);

% for j = 1:len
%     load_file = strcat('../imgset1/', imgs(j,:), '.jpg');
%     imshow(load_file);
%     data.D2 = ginput(n);
%     data.D3(1:n, 3) = 1;
%     save_file = strcat('ccdata/', imgs(j,:), '.mat');
%     save(save_file, 'data')
% end

%% Load and Plot 2D to find 3D Correspondence

for j = 1:len
    load_img = strcat('../imgset1/', imgs(j,:), '.jpg');
    temp = imread(load_img);
    image(temp);
    load_mat = strcat('ccdata/', imgs(j,:), '.mat');
    load(load_mat);
    for i = 1:n 
        text(data.D2(n,1), data.D2(n,2), num2str(i));
    end
    
%     save_file = strcat('ccdata/', imgs(j,:), '_3D.mat');
%     save(save_file, 'data')
end

%% Hard-coded version of 3D Coordinates


% dscf4177_3D = [
%     % top of yellow block
%      0 48 67; % left
%     32 48 67; % center
%     32 80 67; % right
%     % top of center green block
%     16 16 29; % left
%     48 16 29; % center
%     48 48 29; % right
%     % top of red block
%      0  0 29; % left
%     64  0 29; % center
%     64 64 29; % right
%     % diagonal of checkerboard
%     4*20 2*20 0;
%     5*20 1*20 0;
%     6*20    0 0
%     % far corner of checkerboard
%      9*20  9*20 0;
%     10*20  9*20 0;
%     10*20 10*20 0
% ];

