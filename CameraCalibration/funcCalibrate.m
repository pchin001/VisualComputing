% *************************************************************************
% Title: Function-Produce Calibration Matrix and Average Pixel Error using
% Singular value decomposition (SVD)
% Author: Siddhant Ahuja
% Created: February 2010
% Copyright Siddhant Ahuja, 2010
% ***Inputs***
% 1. File3D: 3-D File containing world coordinates
% 2. File2D: 2-D File containing pixel coordinates
% ***Outputs***
% calibMatrix: 3x4 Calibration Matrix
% rotationMatrix: 3x3 Rotation Matrix
% translationVector: 3x1 Translation vector
% alpha_u: Parameter for scaling in X-direction
% alpha_v: Parameter for scaling in Y-direction
% u_0: Optical Centre of the image (x value)
% v_0: Optical Centre of the image (y value)
% reprojMatrix: Reprojection of 3D points to 2D using Calibration Matrix
% avgError_u: Average Pixel Error in x direction
% avgError_v: Average Pixel Error in y direction
% timeTaken: Time taken by the code
% Example Usage of Function: 
% [calibMatrix, rotationMatrix, translationVector, alpha_u, alpha_v, u_0,
% v_0, reprojMatrix, avgError_u, avgError_v, timeTaken]= funcCalibrate('3D.txt', '2D.txt');
% *************************************************************************
function [calibMatrix, A, rotationMatrix, translationVector, reprojMatrix, avgError_u, avgError_v, timeTaken] = funcCalibrate(File3D, File2D)
 
isPlot = 0;

% Assuming no Radial Distortion is present.
% Initialize the timer to calculate the time consumed.
tic;
matrix3D = File3D;
% Read 3-D file
% [matrix3D]=textread(File3D);
% Delete First row of the 3D Matrix as it only contains number of data points
% matrix3D(1,:)=[]; 
% Find the size of the 3D Matrix
[nR3D, nC3D]=size(matrix3D);
% Check to make sure matrix has 3 columns
if(nC3D~=3)
    error('The matrix for 3D points does not have 3 columns.');
end
% Separate out the X values from the matrix
X = matrix3D(:,1);
% Separate out the X values from the matrix
Y = matrix3D(:,2);
% Separate out the X values from the matrix
Z = matrix3D(:,3);
% Plot the points in 3-D
if(isPlot)
    figure;
    plot3(matrix3D(:,1),matrix3D(:,2),matrix3D(:,3),'.');
end
% Read 2-D file
matrix2D = File2D;
% [matrix2D]=textread(File2D);
% Delete First row of the 2D Matrix as it only contains number of data points
% matrix2D(1,:)=[]; 
% Find the size of the 2D Matrix
[nR2D, nC2D]=size(matrix2D);
% Check to make sure matrix has 2 columns
if(nC2D~=2)
    error('The matrix for 2D points does not have 2 columns.');
end
% Separate out the u values from the matrix
u_values=matrix2D(:,1);
% Separate out the v values from the matrix
v_values=matrix2D(:,2);
% Plot the points in 2-D
if(isPlot)
figure;
plot(matrix2D(:,1),matrix2D(:,2),'.');
end
% Check to make sure number of rows of the 3D Matrix is the same as the
% number of rows of the 2D Matrix
if(nR3D~=nR2D)
    error('Please make sure number of 3D and 2D points is the same.');
end
% ***Linear Solution of the Calibration Matrix***
% Let Calibtration Matrix be denoted by (M)
% Writing linear equations of the form AV=0, where A is a 2nx12 measurement
% matrix and V is a 12-element unknown vector.
% Create a Matrix of ones (o) with the same length as that of the u_values
% vector
o = ones(size(u_values));
% Create a Matrix of zeros (z) with the same length as that of the u_values
% vector
z = zeros(size(u_values));
% Populate the odd rows for A
AoddRows  = [ X Y Z o z z z z -u_values.*X -u_values.*Y -u_values.*Z -u_values ];
% Populate the even rows for A
AevenRows = [ z z z z X Y Z o -v_values.*X -v_values.*Y -v_values.*Z -v_values ];
% Concatenate odd and even rows of A
A=[AoddRows; AevenRows];
% Compute SVD on the matrix
[U, S, V] = svd(A,0);
% Assuming no noise, since the elements of the diagonal matrix S are in descending order, to
% get the eigenvectors corresponding to the smallest eignevalue, we can
% just grab the last column of matrix 
m = V(:,end);
% Construct the camera calibration matrix M
M = reshape(m,4,3)';
calibMatrix=M;
% Since the norm of Projection matrix M is equal to 1, we can calculate the
% absolute scale factor lambda
abs_lambda=sqrt(M(3,1)^2 + M(3,2)^2 + M(3,3)^2);
% Scale the Matrix with the scale factor
M = M / abs_lambda;
% In case the origin of the world frame is in front of the camera, we have
% s=lambda/abs_lambda, or From T_z=s*m_34, we can re-write s in the form of
% sign value of m_34. Here we assume it is in the front.
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
translationVector=T;
% The rotation matrix R obtained with this estimation procedure is not guaranteed to be orthogonal. 
% Therefore we calculate the rotation matrix that is closest to the estimated matrix (in the 
% Frobenius norm sense). Let R = UDV'T then R = UV'T.
[U,D,V] = svd(R);
R = U*V';
rotationMatrix=R;
% Reproject 3Dpoints to 2D points using the calibration matrix to calculate average pixel errors 
newMatrix2D=zeros(nR3D,2);
for i=1:1:nR3D
    num_u=M(1,1)*matrix3D(i,1) + M(1,2)*matrix3D(i,2) + M(1,3)*matrix3D(i,3) + M(1,4);
    num_v=M(2,1)*matrix3D(i,1) + M(2,2)*matrix3D(i,2) + M(2,3)*matrix3D(i,3) + M(2,4);
    den=M(3,1)*matrix3D(i,1) + M(3,2)*matrix3D(i,2) + M(3,3)*matrix3D(i,3) + M(3,4);
    newMatrix2D(i,1)=num_u/den;
    newMatrix2D(i,2)=num_v/den;
end
% Reprojected matrix
reprojMatrix=newMatrix2D;
% Calculate difference between the reprojectedMatrix and the original 2D
% Matrix
errorDiff=reprojMatrix-matrix2D;
% Calculate average pixel error in x direction
avgError_u=mean(errorDiff(:,1));
% Calculate average pixel error in y direction
avgError_v=mean(errorDiff(:,2));
% Stop the timer to calculate the time consumed.
timeTaken=toc;


A = [
  alpha_u   0           u_0;
  0         alpha_v     v_0;
  0         0           1;
];
