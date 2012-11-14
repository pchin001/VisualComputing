%% Part A: Camera Calibration
clc
clear all

convertXY = 0;
width = 2848;
height = 2136;


%% Load the previously saved data
load('dscf4177_pixel.mat');
D2 = dscf4177.D2(1:end,:)';
D3 = dscf4177.D3(1:end,:)';


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
if(convertXY)
    D2(1,:) = D2(1,:) - width/2;
    D2(2,:) = -D2(2,:) + height/2;
end


%% Create linear system (matrix)
[~,len] = size(D2);
linsys = zeros(len*2,12);
for i = 1:len
    row = 2*i - 1;
    % u coordinate row
    linsys(row,1:4) = D3(:,i);
    linsys(row,9:12) = -D2(1,i) .* D3(1:4,i);
    % v coordinate row
    linsys(row+1,5:8) = D3(:,i);
    linsys(row+1,9:12) = -D2(2,i) .* D3(1:4,i);
end


%% Solve using SVD
[U,S,V] = svd(linsys,0);
if(length(V(:,end)) == 11)
    c_elems = [V(:,end); 1];
else
    c_elems = V(:,end);
end
C = reshape(c_elems,4,3)';

% Scale the Matrix with the scale factor
% C = C ./ sqrt(C(3,1)^2 + C(3,2)^2 + C(3,3)^2);
% C = C ./ C(end,end);
C


%% Test the C matrix

% use original input as test points
% test3D = D3; % original input
% test2D = D2;

% use corners as test points
test3D(:,1) = [0 0 0 1]'; % origin
test2D(:,1) = [222 1563 1]';
test3D(:,2) = [20 180 0 1]'; % top left corner
test2D(:,2) = [1243 669 1]';
test3D(:,3) = [180 200 0 1]'; % top right corner
test2D(:,3) = [2455 885 1]';

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
avgerr = sum(avgerr) ./ count; % average error


%% Find Position and Orientation of each C-Matrix
A = C(1:3,1:3);

% find R and K using RQ-Decomposition
% this implements RQ Decomposition using QR Decomposition
ReverseRows = [0 0 1; 0 1 0 ; 1 0 0];
[R K] = qr((ReverseRows * A)');
K = ReverseRows * K' * ReverseRows;
R = ReverseRows * R';

% find R and K using QR-Decomposition
% [R,K] = qr(A); % R is Rotation Matrix, K is Intrinsic Properies

% normalize K so that last cell is 1
K = K ./ K(end:end);

% print K and R for verification
K;
R;
asin(-R(3,1)) * 180 / pi; % rotation around y-axis

% find T
% T = inv(R) * inv(K) * C(:,end) % KRT = C(last column)
T = inv(A) * C(:,end) % A*T = C(last column)


%% Using online code to calculate T
% http://siddhantahuja.wordpress.com/2010/02/20/570/

% setup
M = C;

% begin code!

inFront=1;
if inFront
    s = sign(M(3,4));
else
    s = -sign(M(3,4));
end
% Thus, we can now calculate T_z or T(3)
T(3) = s*M(3,4);
% Create a 3x3 Rotation matrix and fill it with zeros
R = zeros(3,3);
% From equations, last row of the rotation matrix is the same as the first
% three elements of the last row of the calibration matrix
R(3,:)=s*M(3,1:3);
% Matrix M can be written as:
% M=( m1'   )
%   ( m2' m4)
%   ( m3'   )
% We can now calculate mi, where mi is a 3 element vector
m1 = M(1,1:3)';
m2 = M(2,1:3)';
m3 = M(3,1:3)';
m4 = M(1:3,4);
% Now, we can calculate the centres of projection u_0 and v_0
u_0 = m1'*m3;
v_0 = m2'*m3;
% Calculating the alpha values in u and v directions,
alpha_u=sqrt( m1'*m1 - u_0^2 );
alpha_v=sqrt( m2'*m2 - v_0^2 );
% We can now calculate the first and second rows of the rotation matrix
R(1,:) = s*(u_0*M(3,1:3) - M(1,1:3) ) / alpha_u;
R(2,:) = s*(v_0*M(3,1:3) - M(2,1:3) ) / alpha_v;
% We can also calculate the first and second elements of the Translation
% vector
T(1) = s*(u_0*M(3,4) - M(1,4) ) / alpha_u;
T(2) = s*(v_0*M(3,4) - M(2,4) ) / alpha_v;
T = T';
translationVector=T
% The rotation matrix R obtained with this estimation procedure is not guaranteed to be orthogonal.
% Therefore we calculate the rotation matrix that is closest to the estimated matrix (in the
% Frobenius norm sense). Let R = UDV'T then R = UV'T.
[U,D,V] = svd(R);
R = U*V';
rotationMatrix=R;
