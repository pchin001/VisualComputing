function [ C, A, R, T ] = CameraCalibration( D2, D3 )
%CAMERACALIBRATION Summary of this function goes here
%   D2,D3: each column is a point


%% Create linear system (matrix)
[~,len] = size(D2);
linsys = zeros(len*2,12);
for i = 1:len
    row = 2*i - 1;
    % u coordinate row
    linsys(row,1:4) = D3(:,i);
    linsys(row,9:12) = -D2(1,i) .* D3(:,i);
    % v coordinate row
    linsys(row+1,5:8) = D3(:,i);
    linsys(row+1,9:12) = -D2(2,i) .* D3(:,i);
end


%% Solve calibration matrix using SVD

% 12 unknowns
[U,S,V] = svd(linsys);
C = reshape(V(:,end),4,3)';


%% Scale the camera matrix with a scale factor

% this scaling is singularity free
% C = C ./ sqrt(C(3,1)^2 + C(3,2)^2 + C(3,3)^2);

% this scaling my produce a singular matrix if last cell is small
% C = C ./ C(end,end);


%% Find A,R,T using Majumder and Wiki
B = C(1:3,1:3);

% find A and R using RQ-Decomposition from online
[wikiA,wikiR] = rq(B);

% normalize A so that last cell is 1
wikiA = wikiA ./ wikiA(end:end);

% find T, where A*R*T = C(last column)
wikiT = inv(B) * C(:,end); % since B = A*R

% wikiT is the position of 3D origin (0,0,0) with respect to the camera's
% coordinate system. We reverse it here to find the camera's position with
% respect to the global 3D coordinate system.
wikiT = -wikiT;


%% Switch some angles to make it correct

% http://nghiaho.com/?page_id=846
% http://www.robertblum.com/articles/2005/02/14/decomposing-matrices

% The intrinsic properties matrix, A, is not the same across images. I
% suspect that RQ Decomposition is giving us a rotation matrix with one or
% more of the rotations negated. Fortunately, for the data set I am using,
% RQ Decomposition seems to consistently negate the X Rotation.

% decompose the rotation matrix into X,Y,Z rotations
thetaX = atan2(wikiR(3,2),wikiR(3,3));
thetaY = atan2(-wikiR(3,1),sqrt(wikiR(3,2)^2+wikiR(3,3)^2));
thetaZ = atan2(wikiR(2,1),wikiR(1,1));

% choose which rotation to negate
thetaX = -thetaX;
% thetaY = -thetaY;
% thetaZ = -thetaZ;

% recombine the rotation matrices
wikiR = Rz(thetaZ)' * Ry(thetaY)' * Rx(thetaX)';


%% The graveyard of test code

% thetaX = atan2(Rnew(3,2),Rnew(3,3)) / pi * 180;
% thetaY = atan2(-Rnew(3,1),sqrt(Rnew(3,2)^2+Rnew(3,3)^2)) / pi * 180;
% thetaZ = atan2(Rnew(2,1),Rnew(1,1)) / pi * 180;

% 
% % C = A * R;
% C(1:3,1:3);
% Cnew = wikiA*Rnew;
% % X = C(1:3,1:3) * inv(Cnew);
% X = R * inv(Rnew);
% % R = X * Rnew
% Ccheck = wikiA * X * Rnew;
% % Cnew = A * X * Rnew
% % C = X * Cnew
% % C = A * X * Rnew          replacement with Cnew
% % C = X * newA * newR       RQ Decompose Cnew
% [newR, newA] = qr(Cnew);
% % wikiA = wikiA ./ wikiA(end,end);
% newA = X * newA;
% % newA = newA ./ newA(end,end);
% newR;
% X * newA * newR


% thetaX = atan2(newR(3,2),newR(3,3))
% thetaY = atan2(-newR(3,1),sqrt(newR(3,2)^2+newR(3,3)^2))
% thetaZ = atan2(newR(2,1),newR(1,1))


    
%% Return

A = wikiA;
R = wikiR;
T = wikiT;


end

