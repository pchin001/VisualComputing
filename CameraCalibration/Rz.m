% http://www.dept.aoe.vt.edu/~cdhall/courses/matlab/
% Row-Major matrix. Transpose to get Column-Major
function Rotmat = Rz(angle)
   c = cos(angle); s = sin(angle);
   Rotmat = [ c s 0; -s c 0; 0 0 1];
