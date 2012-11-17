function [ ] = PlotCamerasInScene( rotation, translation )
% 'rotation' and 'translation' need to be (4x4xn) matrices
% n is the number of cameras you want to plot
[~,~,n] = size(rotation);

%% assign plot properties
figure;
hold on;
xlabel('x');
ylabel('y');
zlabel('z');
axis([-200 400 -200 400 0 600]);


%% plot the scene references

% plot the planar checkerboard
plot3([0 200 200 0 0],[0 0 200 200 0],[0 0 0 0 0]);

% plot the blocks at the origin
trans = .85; % transparency
plotcube([64 64 19],[ 0  0  0],trans,[0 0 1]); % blue
plotcube([64 64 10],[ 0  0 19],trans,[1 0 0]); % red
plotcube([32 32 19],[16 16 29],trans,[0 1 0]); % center green
plotcube([32 32 19],[ 0 48 29],trans,[0 1 0]); % corner green
plotcube([32 32 19],[ 0 48 48],trans,[1 1 0]); % yellow


%% calculate and plot each camera

vlen = 100;
camera = zeros(4,2,3);

for i = 1:n
    % define xyz unit vectors
    camera(:,:,1) = [0 0 0 1;  vlen  0    0    1]';
    camera(:,:,2) = [0 0 0 1;   0   vlen  0    1]';
    camera(:,:,3) = [0 0 0 1;   0    0   vlen  1]';
    
    % for each unit vector
    for j = 3:3
        % rotate and translate
        camera(:,:,j) = translation(:,:,i)*rotation(:,:,i)*camera(:,:,j);
        % normalize
        camera(:,:,j) = camera(:,:,j) ./ repmat(camera(4,:,j),4,1);
        % plot arrow
        arrow(camera(1:3,1,j)', camera(1:3,2,j)');
    end
end

axis equal;


end

