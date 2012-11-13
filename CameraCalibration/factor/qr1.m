%QR1     Orthogonal-triangular decomposition.
%   [Q,R] = QR1(A), where A is m-by-n, produces an m-by-n upper triangular
%   matrix R and an m-by-m unitary matrix Q so that A = Q*R.
%
%   [Q,R] = QR1(A,0) produces the "economy size" decomposition.
%   If m>n, only the first n columns of Q and the first n rows of R are
%   computed. If m<=n, this is the same as [Q,R] = QR1(A).
%
%   X = QR1(A) and X = QR1(A,0) return the output of LAPACK's *GEQRF
%   routine. TRIU(X) is the upper triangular factor R.
%
%   See also QR.