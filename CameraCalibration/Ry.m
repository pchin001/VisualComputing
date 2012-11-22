% http://www.dept.aoe.vt.edu/~cdhall/courses/matlab/
function Rotmat = Ry(angle)
   c = cos(angle); s = sin(angle);
   Rotmat = [ c 0 -s; 0 1 0; s 0 c];
