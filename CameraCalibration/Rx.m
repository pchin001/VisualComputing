% http://www.dept.aoe.vt.edu/~cdhall/courses/matlab/
% Row-Major matrix. Transpose to get Column-Major
function Rotmat = Rx(angle)
   c = cos(angle); s = sin(angle);
   Rotmat = [ 1 0 0; 0 c s; 0 -s c ];
