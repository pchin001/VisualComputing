% http://www.dept.aoe.vt.edu/~cdhall/courses/matlab/
function Rotmat = Rz(angle)
   c = cos(angle); s = sin(angle);
   Rotmat = [ c s 0; -s c 0; 0 0 1];
