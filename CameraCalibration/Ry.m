% http://www.dept.aoe.vt.edu/~cdhall/courses/matlab/
% Row-Major matrix. Transpose to get Column-Major
function Rotmat = Ry(angle)
   c = cos(angle); s = sin(angle);
   Rotmat = [ c 0 -s; 0 1 0; s 0 c];
