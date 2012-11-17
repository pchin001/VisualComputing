function [ output_args ] = funcCalibrate2( w, i )
%FUNCCALIBRATE2 Summary of this function goes here
% http://www.sci.utah.edu/~gerig/CS6320-S2012/Materials/hw1_abhishek_projectreport_b.pdf


%% world coordinates in cms

%z coordinates
wz = w(3,:);
% wz(1:6) = 7*2.8; wz(31:36) = 7*2.8;
% wz(7:12) = 6*2.8; wz(37:42) = 6*2.8;
% wz(13:18) = 5*2.8; wz(43:48) = 5*2.8;
% wz(19:24) = 4*2.8; wz(49:54) = 4*2.8;
% wz(25:30) = 3*2.8; wz(55:60) = 3*2.8;

%x coordinates
wx = w(1,:);
% wx(1:60) = 0;
% wx(1:6:30) = 0.9 + 6*2.8;
% wx(2:6:30) = 0.9 + 5*2.8;
% wx(3:6:30) = 0.9 + 4*2.8;
% wx(4:6:30) = 0.9 + 3*2.8;
% wx(5:6:30) = 0.9 + 2*2.8;
% wx(6:6:30) = 0.9 + 1*2.8;

%y coordinates
wy = w(2,:);
% wy(1:60) = 0;
% wy(31:6:60) = 0.9 + 1*2.8;
% wy(32:6:60) = 0.9 + 2*2.8;
% wy(33:6:60) = 0.9 + 3*2.8;
% wy(34:6:60) = 0.9 + 4*2.8;
% wy(35:6:60) = 0.9 + 5*2.8;
% wy(36:6:60) = 0.9 + 6*2.8;


%% image coordinates

ix = i(1,:);

% ix = [0.4020 0.5177 0.6333 0.7365 0.8294 0.9202 0.4062 0.5218 0.6333 0.7386 ...
% 0.8294 0.9202 0.4124 0.5239 0.6333 0.7365 0.8273 0.9202 0.4144 0.5300 ...
% 0.6312 0.7324 0.8315 0.9202 0.4144 0.5280 0.6353 0.7324 0.8294 0.9161 ...
% 1.1700 1.2712 1.3724 1.4942 1.6119 1.7419 1.1659 1.2650 1.3724 1.4859 ...
% 1.6057 1.7337 1.1618 1.2588 1.3662 1.4818 1.5995 1.7275 1.1577 1.2547 ...
% 1.3620 1.4735 1.5912 1.7171 1.1535 1.2506 1.3538 1.4694 1.5850 1.7130] * 1000;

iy = i(2,:);

% iy = [0.4222 0.4168 0.4168 0.4141 0.4168 0.4195 0.5731 0.5677 0.5570 0.5543 ...
% 0.5489 0.5462 0.7213 0.7079 0.6998 0.6863 0.6782 0.6701 0.8696 0.8534 ...
% 0.8345 0.8210 0.8076 0.7941 1.0124 0.9908 0.9720 0.9504 0.9342 0.9181 ...
% 0.4195 0.4222 0.4276 0.4276 0.4330 0.4384 0.5435 0.5516 0.5597 0.5704 ...
% 0.5785 0.5866 0.6728 0.6836 0.6917 0.7052 0.7186 0.7348 0.7968 0.8103 ...
% 0.8237 0.8399 0.8588 0.8722 0.9181 0.9342 0.9531 0.9720 0.9908 1.0151] * 1000;

n = length(i(2,:)); %number of points
P(1:2*n,1:12) = 0;
j=1;


%% construct matrix P

for i=1:2:120
P(i,1) = wx(j); P(i+1,5) = wx(j);
P(i,2) = wy(j); P(i+1,6) = wy(j);
P(i,3) = wz(j); P(i+1,7) = wz(j);
P(i,4) = 1; P(i+1,8) = 1;
P(i,9:12) = P(i,1:4)*-1*ix(j);
P(i+1,9:12) = P(i,1:4)*-1*iy(j);
j = j+1;
end


%% Perform SVD of P
[U S V] = svd(P);
[min_val, min_index] = min(diag(S(1:12,1:12)));
%m is given by right singular vector of min. singular value
m = V(1:12, min_index);
%normalize M to make the norm of third rotation vecto unity
norm_31 = norm(m(9:11));
m_canonical = m / norm_31;
M(1,1:4) = m_canonical(1:4);
M(2,1:4) = m_canonical(5:8);
M(3,1:4) = m_canonical(9:12);
5
a1 = M(1,1:3);
a2 = M(2,1:3);
a3 = M(3,1:3);
b = M(1:3, 4);
r3 = a3;


%% compute the intrinsic parameters
u_0 = a1*a3';
v_0 = a2*a3';
cross_a1a3 = cross(a1,a3);
cross_a2a3 = cross(a2,a3);
theta = acos (-1*cross_a1a3*cross_a2a3'/(norm(cross_a1a3)*norm(cross_a2a3)));
alpha = norm(cross_a1a3) * sin(theta);
beta = norm(cross_a2a3) * sin(theta);


%% compute the extrinsic parameters
r1 = cross_a2a3/norm(cross_a2a3);
r2 = cross(r3, r1);
K = [alpha -1*alpha*cot(theta) u_0
0 beta/sin(theta) v_0
0 0 1];
t = inv(K) * b; %translation vector


%% rotation matrix
R(1,1:3) = r1;
R(2,1:3) = r2;
R(3,1:3) = r3;
%%Test of calibration estimates: reconstruction error of 4 new points
%image coordinates (measured)
ix_test_measured = 1.0e+003 *[0.5177 0.7365 1.2485 1.3517];
iy_test_measured = 1.0e+003 *[0.2740 0.2794 1.0609 1.0797];
%world coordinates (measured)
wx_test_measured = [0.9+5*2.8 0.9+3*2.8 0 0];
wy_test_measured = [0 0 0.9+2*2.8 0.9+3*2.8];
wz_test_measured = [8*2.8 8*2.8 2*2.8 2*2.8];


%% reconstruct image coordinates and calculate estimation error
for i=1:4
temp(1:4) = [wx_test_measured(i) wy_test_measured(i) wz_test_measured(i) 1];
image_reconstructed = M * temp';
image_recons_x(i) = image_reconstructed(1)/image_reconstructed(3);
image_recons_y(i) = image_reconstructed(2)/image_reconstructed(3);
error(i) = norm([image_recons_x(i)-ix_test_measured(i) image_recons_y(i)-iy_test_measured(i)]);
end
%print all the values
%intrinsic
theta
u_0
v_0
alpha
beta
%extrinsic
R
t
%reconstruction error
error
6


end

