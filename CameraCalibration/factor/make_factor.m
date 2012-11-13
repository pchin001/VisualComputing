%MAKE_FACTOR Compilation of QR/QL/LQ/RQ mex-files
MATLAB_PATH = matlabroot;
COMPILE_OPTIONS = '';
v = ver('matlab');
matver = sscanf(v.Version, '%d.%d.%d')';
COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -DMATLAB_VERSION=0x' sprintf('%02d%02d', matver(1), matver(2)) ];
MATLAB_VERSION = matver(2);

if strcmpi('GLNX86', computer) || strcmpi('GLNXA64', computer) ...
        || strcmpi('MACI', computer) || strcmpi('MAC', computer) ...
        || strcmpi('MACI64', computer)
    % GNU/Linux (x86-32 or x86-64) or MacOS (Intel or PPC)
    LAPACK_PATH = ' -lmwlapack';
    if MATLAB_VERSION < 5;
        BLAS_PATH = '';
    else
        BLAS_PATH = ' -lmwblas';
    end
    if strcmpi('GLNX86', computer) || strcmpi('GLNXA64', computer)
        COMPILE_OPTIONS = [COMPILE_OPTIONS,' -DSkip_f2c_Undefs',' -DNON_UNIX_STDIO'];
    else
        COMPILE_OPTIONS = [COMPILE_OPTIONS,' -DSkip_f2c_Undefs'];
    end
elseif strcmpi('PCWIN', computer) || strcmpi('PCWIN64', computer)
    % Windows (x86-32 or x86-64)
    if strcmpi('PCWIN', computer)
        if MATLAB_VERSION < 6;
            MANUFACTURER = 'lcc';
        else
            cc = mex.getCompilerConfigurations('Any','Selected');
            MANUFACTURER = cc.Manufacturer;
        end
        switch lower(MANUFACTURER)
            case {'lcc'}
                BLAS_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win32', 'lcc', 'libmwblas.lib');
                LAPACK_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win32', 'lcc', 'libmwlapack.lib');
                COMPILE_OPTIONS = [COMPILE_OPTIONS,' -DLCCWIN32'];
            case {'microsoft'}
                BLAS_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win32', 'microsoft', 'libmwblas.lib');
                LAPACK_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win32', 'microsoft', 'libmwlapack.lib');
            case {'sybase'}
                BLAS_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win32', 'watcom', 'libmwblas.lib');
                LAPACK_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win32', 'watcom', 'libmwlapack.lib');
                COMPILE_OPTIONS = [COMPILE_OPTIONS,' -DWATCOMWIN32'];
            otherwise
                disp('Try "mex -setup", because BLAS/LAPACK library is not available!')
        end
    else
        BLAS_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win64', 'microsoft', 'libmwblas.lib');
        LAPACK_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win64', 'microsoft', 'libmwlapack.lib');
    end
    if MATLAB_VERSION < 5;
        BLAS_PATH = ''; % On <= 7.4, BLAS in included in LAPACK
    else
        BLAS_PATH = [' "' BLAS_PATH '"'];
    end
    LAPACK_PATH = [' "' LAPACK_PATH '"'];
    COMPILE_OPTIONS = [COMPILE_OPTIONS,' -DMSDOS',' -DUSE_CLOCK',' -DNO_ONEXIT'];
else
    error('Unsupported platform')
end

% Large array dims for 64 bits platforms appeared in Matlab 7.3
if (strcmpi('GLNXA64', computer) || strcmpi('PCWIN64', computer) ...
        || strcmpi('MACI64', computer)) ...
        && ~(MATLAB_VERSION < 3)
    COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -largeArrayDims' ];
end

% Comment next line to suppress optimization
COMPILE_OPTIONS = [ ' -O' COMPILE_OPTIONS ];

% Comment next line to suppress compilation debugging info
%COMPILE_OPTIONS = [ ' -v' COMPILE_OPTIONS ];

disp('Compiling lq...')
eval(['mex ', COMPILE_OPTIONS, ' lq.c', BLAS_PATH, LAPACK_PATH]);
disp('Compiling ql...')
eval(['mex ', COMPILE_OPTIONS, ' ql.c', BLAS_PATH, LAPACK_PATH]);
disp('Compiling rq...')
eval(['mex ', COMPILE_OPTIONS, ' rq.c', BLAS_PATH, LAPACK_PATH]);
disp('Compiling qr1...')
eval(['mex ', COMPILE_OPTIONS, ' qr1.c', BLAS_PATH, LAPACK_PATH]);
