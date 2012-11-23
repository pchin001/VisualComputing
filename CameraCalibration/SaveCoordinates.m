%% SaveCoordinates.m
clc
clear all

set = '2';

%% Set 1 Images
if set == '1'
imgs = [
    'DSCF4177';
    'DSCF4179';
    'DSCF4183';
    'DSCF4186';
    'DSCF4188';
    'DSCF4192';
];

n = 24; % number of points to collect per image
[len,~] = size(imgs);
end

%% Set 2 Images

if set == '2'
imgs = [
    'DSCF4184';
    'DSCF4187';
    'DSCF4189';
    'DSCF4195';
    'DSCF4199';
    'DSCF4201';
];

n = 24; % number of points to collect per image
[len,~] = size(imgs);
end

%% Collect coordinates from images
% 
% for j = 1:len
%     load_file = strcat('../imgset', set, '/', imgs(j,:), '.jpg');
%     imshow(load_file);
%     data.D2 = ginput(n);
%     data.D2(:,3) = 1;
%     data.D3(1:n,4) = 1;
%     save_file = strcat('ccdata', set, '/', imgs(j,:), '.mat');
%     save(save_file, 'data')
% end
% 

%% Load and Plot 2D to find 3D Correspondence

for j = 1:len
    load_img = strcat('../imgset', set, '/', imgs(j,:), '.jpg');
    load_mat = strcat('ccdata', set, '/', imgs(j,:), '.mat');
    load(load_mat);
    
    for i = 1:n
    imshow(load_img);
    text(data.D2(i,1), data.D2(i,2), num2str(i), 'color', [1 0 1]);
    pts = input(strcat(num2str(i), ': '));
    data.D3(i,1:4) = [ pts 1 ];
    end
    
    save_file = strcat('ccdata', set, '/', imgs(j,:), '_2D_3D.mat');
    save(save_file, 'data');
end

