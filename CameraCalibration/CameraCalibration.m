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

% Scale the camera matrix with a scale factor
% C = C ./ sqrt(C(3,1)^2 + C(3,2)^2 + C(3,3)^2); % singularity free
C = C ./ C(end,end); % singularity if last cell is small


%% Find A,R,T using Zhang's method
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
zhangR = inv(zhangA) * B; % K-1 B = K-1 K R
zhangT = inv(B) * b;

% zhangT is the position of 3D origin (0,0,0) with respect to the camera's
% coordinate system. We reverse it here to find the camera's position with
% respect to the global 3D coordinate system.
zhangT = -zhangT;
zhangR = -zhangR;

%% Find A,R,T using Majumder and Wiki
B = C(1:3,1:3);

% find A and R using RQ-Decomposition
% this implements RQ Decomposition using QR Decomposition
% http://www.physicsforums.com/showthread.php?t=261739
% ReverseRows = [0 0 1; 0 1 0 ; 1 0 0];
% [wikiR, wikiA] = qr((ReverseRows * B)');
% wikiA = ReverseRows * wikiA' * ReverseRows;
% wikiR = ReverseRows * wikiR';


% find A and R using RQ-Decomposition from online
% http://www.mathworks.com/matlabcentral/fileexchange/24119-dont-let-that-inv-go-past-your-eyes-to-solve-that-system-factorize/content/Factorize/rq.m
[wikiA,wikiR] = rq(B);



% find A and R using QR-Decomposition
% [wikiR,wikiA] = qr(B); % R is Rotation Matrix, A is Intrinsic Properies

% normalize A so that last cell is 1
wikiA = wikiA ./ wikiA(end:end);

% find T
% wikiT = -wikiR' * inv(wikiA) * C(:,end); % A*R*T = C(last column)
wikiT = inv(B) * C(:,end); % same as line above, since B = A*R

% wikiT is the position of 3D origin (0,0,0) with respect to the camera's
% coordinate system. We reverse it here to find the camera's position with
% respect to the global 3D coordinate system.
wikiT = -wikiT;
% wikiR = -wikiR;


%% Return the ones we want

% A = zhangA;
% R = zhangR;
% T = zhangT;

A = wikiA;
R = wikiR;
T = wikiT;

end

