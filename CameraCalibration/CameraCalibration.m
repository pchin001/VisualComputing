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
    linsys(row,9:12) = -D2(1,i) .* D3(1:4,i);
    % v coordinate row
    linsys(row+1,5:8) = D3(:,i);
    linsys(row+1,9:12) = -D2(2,i) .* D3(1:4,i);
end


%% Solve calibration matrix using SVD
[U,S,V] = svd(linsys,0);
C = reshape(V(:,end),4,3)';

% Scale the camera matrix with a scale factor
C = C ./ sqrt(C(3,1)^2 + C(3,2)^2 + C(3,3)^2); % singularity free
% C = C ./ C(end,end); % singularity if last cell is small
C;

%% Find Position and Orientation of each C-Matrix using Zhang's method
% http://cronos.rutgers.edu/~meer/TEACHTOO/PAPERS/zhang.pdf

% setup
B = C(1:3,1:3);
b = C(:,end);
K = B*B';

% find intrinsic matrix, A
uo = K(1,3);
vo = K(2,3);
ku = K(1,1);
kc = K(2,1);
kv = K(2,2);
beta = sqrt(kv - vo^2);
lambda = (kc - uo*vo) / beta;
alpha = sqrt(ku - uo^2 - lambda^2);

% put A together
zhangA = [
    alpha   lambda  uo;
    0       beta    vo;
    0       0       1;
];

% find extrinsic properties of camera, R and T
zhangR = inv(zhangA) * B;
zhangT = inv(zhangA) * b;
zhangT = -zhangR' * zhangT;
zhangT = inv(B) * b;


%% Find Position and Orientation of each C-Matrix
B = C(1:3,1:3);

% find A and R using RQ-Decomposition
% this implements RQ Decomposition using QR Decomposition
ReverseRows = [0 0 1; 0 1 0 ; 1 0 0];
[wikiR, wikiA] = qr((ReverseRows * B)');
wikiA = ReverseRows * wikiA' * ReverseRows;
wikiR = ReverseRows * wikiR';

% find A and R using QR-Decomposition
% [wikiR,wikiA] = qr(B); % R is Rotation Matrix, K is Intrinsic Properies

% normalize A so that last cell is 1
wikiA = wikiA ./ wikiA(end:end);

% find T
% wikiT = -wikiR' * inv(wikiA) * C(:,end); % A*R*T = C(last column)
wikiT = -inv(B) * C(:,end); % same as line above, since B = A*R

% print A, R, and T for verification
wikiA;
wikiR;
wikiT;
% asin(-R(3,1)) * 180 / pi; % rotation around y-axis



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
translationVector=T;
% The rotation matrix R obtained with this estimation procedure is not guaranteed to be orthogonal.
% Therefore we calculate the rotation matrix that is closest to the estimated matrix (in the
% Frobenius norm sense). Let R = UDV'T then R = UV'T.
[U,D,V] = svd(R);
R = U*V';
rotationMatrix=R;

%% Return the ones we want
% 
% A = zhangA;
% R = zhangR;
% T = zhangT;

A = wikiA;
R = wikiR;
T = wikiT;

end

