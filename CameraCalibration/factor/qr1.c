/*
 * QR factorization
 *
 * R = qr1(A)
 * [Q,R] = qr1(A)
 * [Q,R] = qr1(A,0)
 *
 * example compile command (see also make_factor.m):
 * mex -O qr1.c libmwlapack.lib
 * or
 * mex -O qr1.c libmwblas.lib libmwlapack.lib (>= R2007B)
 *
 * calls the SGERQF/DGERQF/CGERQF/ZGERQF and
 * SORGRQ/DORGRQ/CUNGRQ/ZUNGRQ named LAPACK functions
 *
 * Ivo Houtzager
 *
 * Delft Center of Systems and Control
 * The Netherlands, 2010
 */

#include "mex.h"
#include "factor.h"
#include "matrix.h"

void qr_double(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSignedIndex lwork, null, info = 1, econ = 0, cplx = 0, dc = 1;
    double *Qpr, *Rpr, *Ipr, *Qpi, *Rpi, *Ipi, *Ap, *ptau, *pwork, *pwork2, *psize, *psize2;
    mwSignedIndex m, n, min_mn, n2, element_size = sizeof(double);
    mwIndex i, j, limit;
    mxClassID classid = mxDOUBLE_CLASS;
    mxComplexity cplxflag = mxREAL;
    
    /* check complex */
    if (mxIsComplex(prhs[0])) {
        cplxflag = mxCOMPLEX;
        cplx = 1;
        dc = 2;
    }
    
    /* get matrix data */
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    if (m == 0 || n == 0) {
        if (nlhs == 1) {
            plhs[0] = mxCreateNumericMatrix(m,n,classid,cplxflag);
        }
        else {
            plhs[1] = mxCreateNumericMatrix(m,n,classid,cplxflag);
            plhs[0] = mxCreateNumericMatrix(m,m,classid,cplxflag);
            if (m != 0) {
                Qpr = mxGetData(plhs[0]);
                for (j=0; j<m; j++) {
                    Qpr[j*m + j] = 1;
                }
            }
        }
        return;
    }
    if (nrhs == 2) {
        null = mxGetScalar(prhs[1]);
        if (null == 0) {
            econ = 1;
        }
    }

    /* allocate tau */
	min_mn = min(m,n);
    ptau = mxCalloc(dc*min_mn,element_size);

    /* allocate A matrix */
    if (m > n && (econ != 1)) {
        Ap = mxCalloc(dc*m*m,element_size);
    }
    else {
        Ap = mxCalloc(dc*m*n,element_size);
    }
    
    /* copy input to A matrix */
    if (cplx) {
        Ipr = mxGetData(prhs[0]);
        Ipi = mxGetImagData(prhs[0]);
        for (i=0; i<m; i++) {
            for (j=0; j<n; j++) {
                Ap[j*2*m+2*i] = Ipr[j*m+i];
                Ap[j*2*m+2*i+1] = Ipi[j*m+i];
            }
        }
    }
    else {
        Ipr = mxGetData(prhs[0]);
        for (i=0; i<m; i++) {
            for (j=0; j<n; j++) {
                Ap[j*m+i] = Ipr[j*m+i];
            }
        }
    }
    
    /* determine blocksize */
    lwork = -1;
    psize = mxCalloc(dc,element_size);
    if (cplx) {
        zgeqrf(&m, &n, Ap, &m, ptau, psize, &lwork, &info);
    }
    else {
        dgeqrf(&m, &n, Ap, &m, ptau, psize, &lwork, &info);
    }
    lwork = psize[0];
    mxFree(psize);

    /* allocate workspace */
    pwork = mxCalloc(dc*lwork,element_size);

    /* calls the DGEQRF function */
    if (cplx) {
        zgeqrf(&m, &n, Ap, &m, ptau, pwork, &lwork, &info);
    }
    else {
        dgeqrf(&m, &n, Ap, &m, ptau, pwork, &lwork, &info);
    }
    mxFree(pwork);
    if (info != 0) {
        mxFree(ptau);
        mxFree(Ap);
        if (cplx) {
            mexErrMsgTxt("ZGEQRF not succesfull");
        }
        else {
            mexErrMsgTxt("DGEQRF not succesfull");
        }
    }

    /* extract lower triangular part */
    if ((econ == 1) && m > n) {
        if (nlhs == 1) {
            plhs[0] = mxCreateNumericMatrix(n,n,classid,cplxflag);
            if (cplx) {
                Rpr = mxGetData(plhs[0]);
                Rpi = mxGetImagData(plhs[0]);
            }
            else {
                Rpr = mxGetData(plhs[0]);
            }
        }
        else {
            plhs[1] = mxCreateNumericMatrix(n,n,classid,cplxflag);
            if (cplx) {
                Rpr = mxGetData(plhs[1]);
                Rpi = mxGetImagData(plhs[1]);
            }
            else {
                Rpr = mxGetData(plhs[1]);
            }
        }
        if (cplx) {
            for (j=0; j<n; j++) {
                limit = j < min_mn-1 ? j : min_mn-1;
                for (i=0; i<=limit; i++) {
                    Rpr[j*n+i] = Ap[2*j*m+2*i];
                    Rpi[j*n+i] = Ap[2*j*m+2*i+1];
                }
            }
        }
        else {
            for (j=0; j<n; j++) {
                limit = j < min_mn-1 ? j : min_mn-1;
                for (i=0; i<=limit; i++) {
                    Rpr[j*n+i] = Ap[j*m+i];
                }
            }
        }
    }
    else {
        if (nlhs == 1) {
            plhs[0] = mxCreateNumericMatrix(m,n,classid,cplxflag);
            if (cplx) {
                Rpr = mxGetData(plhs[0]);
                Rpi = mxGetImagData(plhs[0]);
            }
            else {
                Rpr = mxGetData(plhs[0]);
            }
        }
        else {
            plhs[1] = mxCreateNumericMatrix(m,n,classid,cplxflag);
            if (cplx) {
                Rpr = mxGetData(plhs[1]);
                Rpi = mxGetImagData(plhs[1]);
            }
            else {
                Rpr = mxGetData(plhs[1]);
            }
        }       
        if (cplx) {
            for (j=0; j<n; j++) {
                limit = j < min_mn-1 ? j : min_mn-1;
                for (i=0; i<=limit; i++) {
                    Rpr[j*m+i] = Ap[2*j*m+2*i];
                    Rpi[j*m+i] = Ap[2*j*m+2*i+1];
                }
            }
        }
        else {
            for (j=0; j<n; j++) {
                limit = j < min_mn-1 ? j : min_mn-1;
                for (i=0; i<=limit; i++) {
                    Rpr[j*m+i] = Ap[j*m+i];
                }
            }
        }
    }

    if (nlhs == 2) {
        n2 = (econ == 1) ? min_mn : m;

        /* determine blocksize */
        lwork = -1;
        psize2 = mxCalloc(dc,element_size);
        if (cplx) {
            zungqr(&m, &n2, &min_mn, Ap, &m, ptau, psize2, &lwork, &info);
        }
        else {
            dorgqr(&m, &n2, &min_mn, Ap, &m, ptau, psize2, &lwork, &info);
        }
        lwork = psize2[0];
        mxFree(psize2);

        /* allocate workspace */
        pwork2 = mxCalloc(dc*lwork,element_size);

        /* calls the DORGQR function */
        if (cplx) {
            zungqr(&m, &n2, &min_mn, Ap, &m, ptau, pwork2, &lwork, &info);
        }
        else {
            dorgqr(&m, &n2, &min_mn, Ap, &m, ptau, pwork2, &lwork, &info);
        }
        mxFree(pwork2);
        if (info != 0) {
            mxFree(ptau);
            mxFree(Ap);
            if (cplx) {
                mexErrMsgTxt("DUNGQR not succesfull");
            }
            else {
                mexErrMsgTxt("DORGQR not succesfull");
            }
        }

        /* allocate Q matrix */
        plhs[0] = mxCreateNumericMatrix(m,n2,classid,cplxflag);
        if (cplx) {
            Qpr = mxGetData(plhs[0]);
            Qpi = mxGetImagData(plhs[0]);

            /* copy Ap to Qp */
            for (i=0; i<m; i++) {
                for (j=0; j<n2; j++) {
                    Qpr[j*m+i] = Ap[j*2*m+2*i];
                    Qpi[j*m+i] = Ap[j*2*m+2*i+1];
                }
            }
        }
        else {
            Qpr = mxGetData(plhs[0]);

            /* copy Ap to Qp */
            for (i=0; i<m; i++) {
                for (j=0; j<n2; j++) {
                    Qpr[j*m+i] = Ap[j*m+i];
                }
            }
        }
    }

    mxFree(ptau);
    mxFree(Ap);
}

void qr_single(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSignedIndex lwork, null, info = 1, econ = 0, cplx = 0, dc = 1;
    float *Qpr, *Rpr, *Ipr, *Qpi, *Rpi, *Ipi, *Ap, *ptau, *pwork, *pwork2, *psize, *psize2;
    mwSignedIndex m, n, min_mn, n2, element_size = sizeof(float);
    mwIndex i, j, limit;
    mxClassID classid = mxSINGLE_CLASS;
    mxComplexity cplxflag = mxREAL;
    
    /* check complex */
    if (mxIsComplex(prhs[0])) {
        cplxflag = mxCOMPLEX;
        cplx = 1;
        dc = 2;
    }
    
    /* get matrix data */
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    if (m == 0 || n == 0) {
        if (nlhs == 1) {
            plhs[0] = mxCreateNumericMatrix(m,n,classid,cplxflag);
        }
        else {
            plhs[1] = mxCreateNumericMatrix(m,n,classid,cplxflag);
            plhs[0] = mxCreateNumericMatrix(m,m,classid,cplxflag);
            if (m != 0) {
                Qpr = mxGetData(plhs[0]);
                for (j=0; j<m; j++) {
                    Qpr[j*m + j] = 1;
                }
            }
        }
        return;
    }
    if (nrhs == 2) {
        null = mxGetScalar(prhs[1]);
        if (null == 0) {
            econ = 1;
        }
    }

    /* allocate tau */
	min_mn = min(m,n);
    ptau = mxCalloc(dc*min_mn,element_size);

    /* allocate A matrix */
    if (m > n && (econ != 1)) {
        Ap = mxCalloc(dc*m*m,element_size);
    }
    else {
        Ap = mxCalloc(dc*m*n,element_size);
    }
    
    /* copy input to A matrix */
    if (cplx) {
        Ipr = mxGetData(prhs[0]);
        Ipi = mxGetImagData(prhs[0]);
        for (i=0; i<m; i++) {
            for (j=0; j<n; j++) {
                Ap[j*2*m+2*i] = Ipr[j*m+i];
                Ap[j*2*m+2*i+1] = Ipi[j*m+i];
            }
        }
    }
    else {
        Ipr = mxGetData(prhs[0]);
        for (i=0; i<m; i++) {
            for (j=0; j<n; j++) {
                Ap[j*m+i] = Ipr[j*m+i];
            }
        }
    }
    
    /* determine blocksize */
    lwork = -1;
    psize = mxCalloc(dc,element_size);
    if (cplx) {
        cgeqrf(&m, &n, Ap, &m, ptau, psize, &lwork, &info);
    }
    else {
        sgeqrf(&m, &n, Ap, &m, ptau, psize, &lwork, &info);
    }
    lwork = psize[0];
    mxFree(psize);

    /* allocate workspace */
    pwork = mxCalloc(dc*lwork,element_size);

    /* calls the DGEQRF function */
    if (cplx) {
        cgeqrf(&m, &n, Ap, &m, ptau, pwork, &lwork, &info);
    }
    else {
        sgeqrf(&m, &n, Ap, &m, ptau, pwork, &lwork, &info);
    }
    mxFree(pwork);
    if (info != 0) {
        mxFree(ptau);
        mxFree(Ap);
        if (cplx) {
            mexErrMsgTxt("ZGEQRF not succesfull");
        }
        else {
            mexErrMsgTxt("DGEQRF not succesfull");
        }
    }

    /* extract lower triangular part */
    if ((econ == 1) && m > n) {
        if (nlhs == 1) {
            plhs[0] = mxCreateNumericMatrix(n,n,classid,cplxflag);
            if (cplx) {
                Rpr = mxGetData(plhs[0]);
                Rpi = mxGetImagData(plhs[0]);
            }
            else {
                Rpr = mxGetData(plhs[0]);
            }
        }
        else {
            plhs[1] = mxCreateNumericMatrix(n,n,classid,cplxflag);
            if (cplx) {
                Rpr = mxGetData(plhs[1]);
                Rpi = mxGetImagData(plhs[1]);
            }
            else {
                Rpr = mxGetData(plhs[1]);
            }
        }
        if (cplx) {
            for (j=0; j<n; j++) {
                limit = j < min_mn-1 ? j : min_mn-1;
                for (i=0; i<=limit; i++) {
                    Rpr[j*n+i] = Ap[2*j*m+2*i];
                    Rpi[j*n+i] = Ap[2*j*m+2*i+1];
                }
            }
        }
        else {
            for (j=0; j<n; j++) {
                limit = j < min_mn-1 ? j : min_mn-1;
                for (i=0; i<=limit; i++) {
                    Rpr[j*n+i] = Ap[j*m+i];
                }
            }
        }
    }
    else {
        if (nlhs == 1) {
            plhs[0] = mxCreateNumericMatrix(m,n,classid,cplxflag);
            if (cplx) {
                Rpr = mxGetData(plhs[0]);
                Rpi = mxGetImagData(plhs[0]);
            }
            else {
                Rpr = mxGetData(plhs[0]);
            }
        }
        else {
            plhs[1] = mxCreateNumericMatrix(m,n,classid,cplxflag);
            if (cplx) {
                Rpr = mxGetData(plhs[1]);
                Rpi = mxGetImagData(plhs[1]);
            }
            else {
                Rpr = mxGetData(plhs[1]);
            }
        }       
        if (cplx) {
            for (j=0; j<n; j++) {
                limit = j < min_mn-1 ? j : min_mn-1;
                for (i=0; i<=limit; i++) {
                    Rpr[j*m+i] = Ap[2*j*m+2*i];
                    Rpi[j*m+i] = Ap[2*j*m+2*i+1];
                }
            }
        }
        else {
            for (j=0; j<n; j++) {
                limit = j < min_mn-1 ? j : min_mn-1;
                for (i=0; i<=limit; i++) {
                    Rpr[j*m+i] = Ap[j*m+i];
                }
            }
        }
    }

    if (nlhs == 2) {
        n2 = (econ == 1) ? min_mn : m;

        /* determine blocksize */
        lwork = -1;
        psize2 = mxCalloc(dc,element_size);
        if (cplx) {
            cungqr(&m, &n2, &min_mn, Ap, &m, ptau, psize2, &lwork, &info);
        }
        else {
            sorgqr(&m, &n2, &min_mn, Ap, &m, ptau, psize2, &lwork, &info);
        }
        lwork = psize2[0];
        mxFree(psize2);

        /* allocate workspace */
        pwork2 = mxCalloc(dc*lwork,element_size);

        /* calls the DORGQR function */
        if (cplx) {
            cungqr(&m, &n2, &min_mn, Ap, &m, ptau, pwork2, &lwork, &info);
        }
        else {
            sorgqr(&m, &n2, &min_mn, Ap, &m, ptau, pwork2, &lwork, &info);
        }
        mxFree(pwork2);
        if (info != 0) {
            mxFree(ptau);
            mxFree(Ap);
            if (cplx) {
                mexErrMsgTxt("DUNGQR not succesfull");
            }
            else {
                mexErrMsgTxt("DORGQR not succesfull");
            }
        }

        /* allocate Q matrix */
        plhs[0] = mxCreateNumericMatrix(m,n2,classid,cplxflag);
        if (cplx) {
            Qpr = mxGetData(plhs[0]);
            Qpi = mxGetImagData(plhs[0]);

            /* copy Ap to Qp */
            for (i=0; i<m; i++) {
                for (j=0; j<n2; j++) {
                    Qpr[j*m+i] = Ap[j*2*m+2*i];
                    Qpi[j*m+i] = Ap[j*2*m+2*i+1];
                }
            }
        }
        else {
            Qpr = mxGetData(plhs[0]);

            /* copy Ap to Qp */
            for (i=0; i<m; i++) {
                for (j=0; j<n2; j++) {
                    Qpr[j*m+i] = Ap[j*m+i];
                }
            }
        }
    }

    mxFree(ptau);
    mxFree(Ap);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* check for proper number of arguments */
    if (nrhs != 1 && nrhs != 2) {
        mexErrMsgTxt("QR1 requires one or two input arguments.");
    }
    if (nlhs > 2) {
        mexErrMsgTxt("Too many output arguments.");
    }
    if (!mxIsNumeric(prhs[0]) || mxIsSparse(prhs[0])) {
        mexErrMsgTxt( "Input must be a full matrix." );
    }
    if (mxIsDouble(prhs[0])) {
        qr_double(nlhs, plhs, nrhs, prhs);
    }
    else if (mxIsSingle(prhs[0])) {
        qr_single(nlhs, plhs, nrhs, prhs);
    }
    else {
        mexErrMsgTxt( "Class is not supported." );
    }
}
