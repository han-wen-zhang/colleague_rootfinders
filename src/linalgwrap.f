cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       linalgwrap.f -- simple wrappers for BLAS/LAPACK routines,
c       for both real*8 (d) and complex*16 (z).
c       All wrappers preserve input unless noted otherwise.
c
c       External dependencies: BLAS, LAPACK
c         (e.g. -framework Accelerate on macOS)
c
c   --- BLAS wrappers ---
c     dmatvecblas(m,n,a,x,y)        y = A*x  (real)
c     zmatvecblas(m,n,a,x,y)        y = A*x  (complex)
c     dmatmatblas(m,n,k,a,b,c)      C = A*B  (real)
c     zmatmatblas(m,n,k,a,b,c)      C = A*B  (complex)
c     zmathmatblas(m,n,k,a,b,c)     C = A^* * B  (complex)
c     zmatmathblas(m,n,k,a,b,c)     C = A * B^*  (complex)
c     dmattmatblas(m,n,k,a,b,c)     C = A^T * B  (real)
c     dmatmattblas(m,n,k,a,b,c)     C = A * B^T  (real)
c
c   --- SVD ---
c     dsvd(m,n,a,u,s,vt,ier)        economy SVD (real)
c     zsvd(m,n,a,u,s,vt,ier)        economy SVD (complex)
c     dsvals(m,n,a,s,ier)            singular values only (real)
c     zsvals(m,n,a,s,ier)            singular values only (complex)
c     dcond(m,n,a,cond,ier)          condition number (real)
c     zcond(m,n,a,cond,ier)          condition number (complex)
c
c   --- QR ---
c     dqr(m,n,a,q,r,ier)            full QR (real)
c     dqrecon(m,n,a,q,r,ier)        economy QR (real)
c     dqreconq(m,n,a,q,ier)        economy QR, Q only (real)
c     zqr(m,n,a,q,r,ier)            full QR (complex)
c     zqrecon(m,n,a,q,r,ier)        economy QR (complex)
c     zqreconq(m,n,a,q,ier)        economy QR, Q only (complex)
c
c   --- Column-pivoted QR ---
c     dcpqr(m,n,a,q,r,ipiv,ier)     economy CPQR with Q (real)
c     dcpqrr(m,n,a,r,ipiv,ier)      CPQR, R and pivot only (real)
c     zcpqr(m,n,a,q,r,ipiv,ier)     economy CPQR with Q (complex)
c     zcpqrr(m,n,a,r,ipiv,ier)      CPQR, R and pivot only (complex)
c
c   --- Eigendecomposition ---
c     dsymeig(n,a,w,v,ier)          symmetric (real)
c     zhermeig(n,a,w,v,ier)         Hermitian (complex)
c     dgeig(n,a,w,v,ier)            general (real -> complex w,v)
c     zgeig(n,a,w,v,ier)            general (complex)
c     dsymeigvals(n,a,w,ier)        symmetric eigenvalues only
c     zhermeigvals(n,a,w,ier)       Hermitian eigenvalues only
c     dgeigvals(n,a,w,ier)          general eigenvalues only (real)
c     zgeigvals(n,a,w,ier)          general eigenvalues only (complex)
c
c   --- Linear solve ---
c     dsolve(n,a,b,x,ier)           solve Ax=b (real)
c     zsolve(n,a,b,x,ier)           solve Ax=b (complex)
c
c   --- Least squares ---
c     dlstsq(m,n,a,b,x,ier)         QR-based, full rank (real)
c     zlstsq(m,n,a,b,x,ier)         QR-based, full rank (complex)
c     dlstsqr(m,n,a,b,x,krank,ier)  SVD-based, robust (real)
c     zlstsqr(m,n,a,b,x,krank,ier)  SVD-based, robust (complex)
c     dqrpack(m,n,a,qr,tau,ier)     precompute QR factorization (real)
c     zqrpack(m,n,a,qr,tau,ier)     precompute QR factorization (complex)
c     dqrsolve(m,n,qr,tau,b,x,ier)  solve using precomputed QR (real)
c     zqrsolve(m,n,qr,tau,b,x,ier)  solve using precomputed QR (complex)
c
c   --- Determinant and inverse ---
c     ddet(n,a,det)                  determinant (real)
c     zdet(n,a,det)                  determinant (complex)
c     dinverse(n,a,ainv,ier)         matrix inverse (real)
c     zinverse(n,a,ainv,ier)         matrix inverse (complex)
c
c   --- LU factorization ---
c     dlu(n,a,lu,ipiv,ier)           factorize, preserve a (real)
c     dlui(n,a,ipiv,ier)             factorize in-place (real)
c     dlusolve(n,lu,ipiv,b,x)        solve, single RHS (real)
c     dlusolven(n,nrhs,lu,ipiv,b,x)  solve, multi RHS (real)
c     zlu(n,a,lu,ipiv,ier)           factorize, preserve a (complex)
c     zlui(n,a,ipiv,ier)             factorize in-place (complex)
c     zlusolve(n,lu,ipiv,b,x)        solve, single RHS (complex)
c     zlusolven(n,nrhs,lu,ipiv,b,x)  solve, multi RHS (complex)
c
c   --- Cholesky factorization (SPD/HPD) ---
c     dchol(n,a,l,ier)               factorize, preserve a (real)
c     dcholi(n,a,ier)                factorize in-place (real)
c     dcholsolve(n,l,b,x)            solve, single RHS (real)
c     dcholsolven(n,nrhs,l,b,x)      solve, multi RHS (real)
c     zchol(n,a,l,ier)               factorize, preserve a (complex)
c     zcholi(n,a,ier)                factorize in-place (complex)
c     zcholsolve(n,l,b,x)            solve, single RHS (complex)
c     zcholsolven(n,nrhs,l,b,x)      solve, multi RHS (complex)
c
c   --- Triangular solve ---
c     dtrisolve(n,t,b,x,iuplo)      solve T*x=b (real)
c     ztrisolve(n,t,b,x,iuplo)      solve T*x=b (complex)
c       iuplo = 0 for lower triangular, 1 for upper triangular
c
c   See also: fortutil.f for pure Fortran utilities (no BLAS/LAPACK)
c
c   --- Common LAPACK info values ---
c     All wrappers print a warning and return ier = info on failure.
c     info > 0: algorithmic failure (common causes below)
c     info < 0: illegal argument (bug in calling code)
c
c     SVD (dgesvd/dgesdd): bidiagonal QR did not converge
c     eig (dsyev/zheev):   tridiagonal QR did not converge
c     eig (dgeev/zgeev):   QR algorithm did not converge
c     solve/LU (dgesv/dgetrf): U(i,i)=0, matrix is singular
c     Cholesky (dpotrf):   leading minor not positive definite
c     lstsq (dgels):       A does not have full rank
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c======= BLAS wrappers and naive equivalents ===================
c
c
        subroutine dmatvecblas(m,n,a,x,y)
c
c       y = a * x   using BLAS dgemv
c
c       input:
c       m -- number of rows of a
c       n -- number of columns of a
c       a(m,n) -- matrix
c       x(n) -- input vector
c
c       output:
c       y(m) -- result vector
c
        implicit none
        integer m,n
        real*8 a(m,n),x(n),y(m)
c
c       dgemv: y = alpha*op(A)*x + beta*y
c       'N' = no transpose, alpha=1, beta=0, incx=incy=1
        call dgemv('N',m,n,1.0d0,a,m,x,1,0.0d0,y,1)
c
        return
        end
c
        subroutine dmatmatblas(m,n,k,a,b,c)
c
c       c = a * b   using BLAS dgemm
c
c       input:
c       m -- number of rows of a and c
c       n -- number of columns of a, number of rows of b
c       k -- number of columns of b and c
c       a(m,n) -- left matrix
c       b(n,k) -- right matrix
c
c       output:
c       c(m,k) -- result matrix
c
        implicit none
        integer m,n,k
        real*8 a(m,n),b(n,k),c(m,k)
c
c       dgemm: C = alpha*op(A)*op(B) + beta*C
c       'N','N' = no transpose on A or B, alpha=1, beta=0
c       args: transa,transb, m,k,n, alpha, A,lda, B,ldb, beta, C,ldc
        call dgemm('N','N',m,k,n,1.0d0,a,m,b,n,0.0d0,c,m)
c
        return
        end
c
        subroutine zmatvecblas(m,n,a,x,y)
c
c       y = a * x   using BLAS zgemv
c
c       input:
c       m -- number of rows of a
c       n -- number of columns of a
c       a(m,n) -- complex*16 matrix
c       x(n) -- complex*16 input vector
c
c       output:
c       y(m) -- complex*16 result vector
c
        implicit none
        integer m,n
        complex*16 a(m,n),x(n),y(m),one,zero
        data one/(1.0d0,0.0d0)/,zero/(0.0d0,0.0d0)/
c
c       zgemv: y = alpha*op(A)*x + beta*y
c       'N' = no conjugate transpose, alpha=1, beta=0, incx=incy=1
        call zgemv('N',m,n,one,a,m,x,1,zero,y,1)
c
        return
        end
c
        subroutine zmatmatblas(m,n,k,a,b,c)
c
c       c = a * b   using BLAS zgemm
c
c       input:
c       m -- number of rows of a and c
c       n -- number of columns of a, number of rows of b
c       k -- number of columns of b and c
c       a(m,n) -- complex*16 left matrix
c       b(n,k) -- complex*16 right matrix
c
c       output:
c       c(m,k) -- complex*16 result matrix
c
        implicit none
        integer m,n,k
        complex*16 a(m,n),b(n,k),c(m,k),one,zero
        data one/(1.0d0,0.0d0)/,zero/(0.0d0,0.0d0)/
c
        call zgemm('N','N',m,k,n,one,a,m,b,n,zero,c,m)
c
        return
        end
c
        subroutine zmathmatblas(m,n,k,a,b,c)
c
c       c = a^* * b   using BLAS zgemm
c       (conjugate transpose of a times b)
c
c       input:
c       m -- number of rows of c (= number of columns of a)
c       n -- number of rows of a, number of rows of b
c       k -- number of columns of b and c
c       a(n,m) -- complex*16 left matrix (stored untransposed)
c       b(n,k) -- complex*16 right matrix
c
c       output:
c       c(m,k) -- complex*16 result matrix
c
        implicit none
        integer m,n,k
        complex*16 a(n,m),b(n,k),c(m,k),one,zero
        data one/(1.0d0,0.0d0)/,zero/(0.0d0,0.0d0)/
c
        call zgemm('C','N',m,k,n,one,a,n,b,n,zero,c,m)
c
        return
        end
c
        subroutine zmatmathblas(m,n,k,a,b,c)
c
c       c = a * b^*   using BLAS zgemm
c       (a times conjugate transpose of b)
c
c       input:
c       m -- number of rows of a and c
c       n -- number of columns of a, number of columns of b
c       k -- number of rows of b, number of columns of c
c       a(m,n) -- complex*16 left matrix
c       b(k,n) -- complex*16 right matrix (stored untransposed)
c
c       output:
c       c(m,k) -- complex*16 result matrix
c
        implicit none
        integer m,n,k
        complex*16 a(m,n),b(k,n),c(m,k),one,zero
        data one/(1.0d0,0.0d0)/,zero/(0.0d0,0.0d0)/
c
        call zgemm('N','C',m,k,n,one,a,m,b,k,zero,c,m)
c
        return
        end
c
        subroutine dmattmatblas(m,n,k,a,b,c)
c
c       c = a^T * b   using BLAS dgemm
c
c       input:
c       m -- number of rows of c (= number of columns of a)
c       n -- number of rows of a, number of rows of b
c       k -- number of columns of b and c
c       a(n,m) -- left matrix (stored untransposed)
c       b(n,k) -- right matrix
c
c       output:
c       c(m,k) -- result matrix
c
        implicit none
        integer m,n,k
        real*8 a(n,m),b(n,k),c(m,k)
c
        call dgemm('T','N',m,k,n,1.0d0,a,n,b,n,0.0d0,c,m)
c
        return
        end
c
        subroutine dmatmattblas(m,n,k,a,b,c)
c
c       c = a * b^T   using BLAS dgemm
c
c       input:
c       m -- number of rows of a and c
c       n -- number of columns of a, number of columns of b
c       k -- number of rows of b, number of columns of c
c       a(m,n) -- left matrix
c       b(k,n) -- right matrix (stored untransposed)
c
c       output:
c       c(m,k) -- result matrix
c
        implicit none
        integer m,n,k
        real*8 a(m,n),b(k,n),c(m,k)
c
        call dgemm('N','T',m,k,n,1.0d0,a,m,b,k,0.0d0,c,m)
c
        return
        end
c
        subroutine dcond(m,n,a,cond,ier)
c
c       Computes the 2-norm condition number of a real*8 matrix
c       cond = s(1)/s(min(m,n)), using LAPACK dgesvd.
c       The input matrix a is not modified.
c
c       input:
c       m -- number of rows of a
c       n -- number of columns of a
c       a(m,n) -- matrix
c
c       output:
c       cond -- condition number (huge if matrix is singular)
c       ier -- 0 if successful; >0: SVD did not converge
c
        implicit none
        integer m,n,mn,ier,info,lwork,j,k
        real*8 a(m,n),cond
        real*8, allocatable :: acopy(:,:),s(:),work(:)
        real*8 wopt
c
        mn = min(m,n)
        allocate(acopy(m,n),s(mn))
c
c       Copy a (dgesvd destroys its input).
c
        do k = 1,n
          do j = 1,m
            acopy(j,k) = a(j,k)
          enddo
        enddo
c
c       Workspace query.
c
c       dgesvd: A = U * diag(S) * V^T
c       'N','N' = do not compute U or V
        call dgesvd('N','N',m,n,acopy,m,s,0,1,0,1,
     1    wopt,-1,info)
        lwork = int(wopt)
        allocate(work(lwork))
c
c       Compute singular values only.
c
        call dgesvd('N','N',m,n,acopy,m,s,0,1,0,1,
     1    work,lwork,info)
c
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: dcond failed, info =', info
c
        if(info .eq. 0) then
          if(s(mn) .gt. 0) then
            cond = s(1)/s(mn)
          else
            cond = huge(1.0d0)
          endif
        else
          cond = -1
        endif
c
        deallocate(acopy, s, work)
c
        return
        end
c
c
c
c
        subroutine zcond(m,n,a,cond,ier)
c
c       Computes the 2-norm condition number of a complex*16 matrix
c       cond = s(1)/s(min(m,n)), using LAPACK zgesvd.
c       The input matrix a is not modified.
c
c       input:
c       m -- number of rows of a
c       n -- number of columns of a
c       a(m,n) -- complex*16 matrix
c
c       output:
c       cond -- condition number (huge if matrix is singular)
c       ier -- 0 if successful; >0: SVD did not converge
c
        implicit none
        integer m,n,mn,ier,info,lwork,j,k
        complex*16 a(m,n)
        real*8 cond
        complex*16, allocatable :: acopy(:,:),work(:)
        real*8, allocatable :: s(:),rwork(:)
        complex*16 wopt
c
        mn = min(m,n)
        allocate(acopy(m,n),s(mn),rwork(5*mn))
c
c       Copy a (zgesvd destroys its input).
c
        do k = 1,n
          do j = 1,m
            acopy(j,k) = a(j,k)
          enddo
        enddo
c
c       Workspace query.
c
c       zgesvd: A = U * diag(S) * V^*
c       'N','N' = do not compute U or V
        call zgesvd('N','N',m,n,acopy,m,s,0,1,0,1,
     1    wopt,-1,rwork,info)
        lwork = int(dble(wopt))
        allocate(work(lwork))
c
c       Compute singular values only.
c
        call zgesvd('N','N',m,n,acopy,m,s,0,1,0,1,
     1    work,lwork,rwork,info)
c
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: zcond failed, info =', info
c
        if(info .eq. 0) then
          if(s(mn) .gt. 0) then
            cond = s(1)/s(mn)
          else
            cond = huge(1.0d0)
          endif
        else
          cond = -1
        endif
c
        deallocate(acopy, s, work, rwork)
c
        return
        end
c
c
c
c
        subroutine dsvd(m,n,a,u,s,vt,ier)
c
c       Economy SVD of a real*8 matrix: A = U * diag(S) * Vt
c       using LAPACK dgesdd. The input matrix a is not modified.
c
c       input:
c       m -- number of rows of a
c       n -- number of columns of a
c       a(m,n) -- matrix
c
c       output:
c       u(m,min(m,n)) -- left singular vectors
c       s(min(m,n)) -- singular values (decreasing order)
c       vt(min(m,n),n) -- right singular vectors (transposed)
c       ier -- 0 if successful; >0: SVD did not converge
c
        implicit none
        integer m,n,mn,ier,info,lwork,j,k
        real*8 a(m,n),u(m,*),s(*),vt(min(m,n),n)
        real*8, allocatable :: acopy(:,:),work(:)
        integer, allocatable :: iwork(:)
        real*8 wopt
c
        mn = min(m,n)
        allocate(acopy(m,n),iwork(8*mn))
c
c       Copy a (dgesdd destroys its input).
c
        do k = 1,n
          do j = 1,m
            acopy(j,k) = a(j,k)
          enddo
        enddo
c
c       Workspace query.
c
c       dgesdd: A = U * diag(S) * Vt
c       'S' = economy SVD (min(m,n) columns of U, rows of Vt)
        call dgesdd('S',m,n,acopy,m,s,u,m,vt,mn,
     1    wopt,-1,iwork,info)
        lwork = int(wopt)
        allocate(work(lwork))
c
c       Compute SVD.
c
        call dgesdd('S',m,n,acopy,m,s,u,m,vt,mn,
     1    work,lwork,iwork,info)
c
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: dsvd failed, info =', info
c
        deallocate(acopy, work, iwork)
c
        return
        end
c
c
c
c
        subroutine zsvd(m,n,a,u,s,vt,ier)
c
c       Economy SVD of a complex*16 matrix: A = U * diag(S) * Vt
c       using LAPACK zgesdd. The input matrix a is not modified.
c       Note: Vt here is V^* (conjugate transpose).
c
c       input:
c       m -- number of rows of a
c       n -- number of columns of a
c       a(m,n) -- complex*16 matrix
c
c       output:
c       u(m,min(m,n)) -- left singular vectors (complex)
c       s(min(m,n)) -- singular values (real, decreasing order)
c       vt(min(m,n),n) -- right singular vectors (conjugate transposed)
c       ier -- 0 if successful; >0: SVD did not converge
c
        implicit none
        integer m,n,mn,ier,info,lwork,j,k
        complex*16 a(m,n),u(m,*),vt(min(m,n),n)
        real*8 s(*)
        complex*16, allocatable :: acopy(:,:),work(:)
        integer, allocatable :: iwork(:)
        real*8, allocatable :: rwork(:)
        complex*16 wopt
c
        mn = min(m,n)
        allocate(acopy(m,n),iwork(8*mn))
        allocate(rwork(5*mn*mn + 7*mn))
c
c       Copy a (zgesdd destroys its input).
c
        do k = 1,n
          do j = 1,m
            acopy(j,k) = a(j,k)
          enddo
        enddo
c
c       Workspace query.
c
c       zgesdd: A = U * diag(S) * V^*
c       'S' = economy SVD
        call zgesdd('S',m,n,acopy,m,s,u,m,vt,mn,
     1    wopt,-1,rwork,iwork,info)
        lwork = int(dble(wopt))
        allocate(work(lwork))
c
c       Compute SVD.
c
        call zgesdd('S',m,n,acopy,m,s,u,m,vt,mn,
     1    work,lwork,rwork,iwork,info)
c
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: zsvd failed, info =', info
c
        deallocate(acopy, work, iwork, rwork)
c
        return
        end
c
c
c======= QR factorization ======================================
c
c
        subroutine dqr(m,n,a,q,r,ier)
c
c       Full QR factorization of a real*8 matrix: A = Q * R
c       using LAPACK dgeqrf + dorgqr.
c       The input matrix a is not modified.
c
c       input:
c       m -- number of rows of a
c       n -- number of columns of a
c       a(m,n) -- matrix
c
c       output:
c       q(m,m) -- orthogonal matrix
c       r(m,n) -- upper triangular matrix
c       ier -- 0 if successful; <0: illegal argument
c
        implicit none
        integer m,n,mn,ier,info,lwork,j,k
        real*8 a(m,n),q(m,m),r(m,n)
        real*8, allocatable :: tau(:),work(:)
        real*8 wopt
c
        mn = min(m,n)
        allocate(tau(mn))
c
c       Copy a into r (dgeqrf overwrites its input).
c
        do k = 1,n
          do j = 1,m
            r(j,k) = a(j,k)
          enddo
        enddo
c
c       QR factorization.
c       dgeqrf: A = Q*R via Householder reflectors stored in A and tau
c       lwork=-1 is workspace query
c
        call dgeqrf(m,n,r,m,tau,wopt,-1,info)
        lwork = int(wopt)
        allocate(work(lwork))
        call dgeqrf(m,n,r,m,tau,work,lwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: dqr failed, info =', info
        if(info .ne. 0) then
          deallocate(tau,work)
          return
        endif
c
c       Copy R into q, then form Q from Householder reflectors.
c
        do k = 1,m
          do j = 1,m
            q(j,k) = 0
          enddo
        enddo
        do k = 1,mn
          do j = 1,m
            q(j,k) = r(j,k)
          enddo
        enddo
c
c       dorgqr: forms Q explicitly from Householder reflectors
c       args: m, n_cols_Q, k_reflectors, Q, ldq, tau, work, lwork
c
        deallocate(work)
        call dorgqr(m,m,mn,q,m,tau,wopt,-1,info)
        lwork = int(wopt)
        allocate(work(lwork))
        call dorgqr(m,m,mn,q,m,tau,work,lwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: dqr failed, info =', info
c
c       Zero out the lower triangle of r.
c
        do k = 1,n
          do j = k+1,m
            r(j,k) = 0
          enddo
        enddo
c
        deallocate(tau,work)
c
        return
        end
c
c
c
c
        subroutine dqrecon(m,n,a,q,r,ier)
c
c       Economy QR factorization of a real*8 matrix: A = Q * R
c       using LAPACK dgeqrf + dorgqr.
c       The input matrix a is not modified.
c
c       input:
c       m -- number of rows of a
c       n -- number of columns of a
c       a(m,n) -- matrix
c
c       output:
c       q(m,min(m,n)) -- orthonormal columns
c       r(min(m,n),n) -- upper triangular matrix
c       ier -- 0 if successful; <0: illegal argument
c
        implicit none
        integer m,n,mn,ier,info,lwork,j,k
        real*8 a(m,n),q(m,*),r(min(m,n),n)
        real*8, allocatable :: acopy(:,:),tau(:),work(:)
        real*8 wopt
c
        mn = min(m,n)
        allocate(acopy(m,n),tau(mn))
c
c       Copy a (dgeqrf overwrites its input).
c
        do k = 1,n
          do j = 1,m
            acopy(j,k) = a(j,k)
          enddo
        enddo
c
c       QR factorization (see dqr for parameter explanation).
c
        call dgeqrf(m,n,acopy,m,tau,wopt,-1,info)
        lwork = int(wopt)
        allocate(work(lwork))
        call dgeqrf(m,n,acopy,m,tau,work,lwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: dqrecon failed, info =', info
        if(info .ne. 0) then
          deallocate(acopy,tau,work)
          return
        endif
c
c       Extract R (upper triangle of acopy, first mn rows).
c
        do k = 1,n
          do j = 1,mn
            if(j .le. k) then
              r(j,k) = acopy(j,k)
            else
              r(j,k) = 0
            endif
          enddo
        enddo
c
c       Form Q (first mn columns) from Householder reflectors.
c       dorgqr with n_cols_Q = mn gives economy Q.
c
        do k = 1,mn
          do j = 1,m
            q(j,k) = acopy(j,k)
          enddo
        enddo
c
        deallocate(work)
        call dorgqr(m,mn,mn,q,m,tau,wopt,-1,info)
        lwork = int(wopt)
        allocate(work(lwork))
        call dorgqr(m,mn,mn,q,m,tau,work,lwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: dqrecon failed, info =', info
c
        deallocate(acopy,tau,work)
c
        return
        end
c
c
c
c
        subroutine dqreconq(m,n,a,q,ier)
c
c       Economy QR factorization of a real*8 matrix, returning
c       only Q: A = Q * R, where Q(m,min(m,n)) has orthonormal
c       columns. R is not computed.
c       The input matrix a is not modified.
c
c       input:
c       m -- number of rows of a
c       n -- number of columns of a
c       a(m,n) -- matrix
c
c       output:
c       q(m,min(m,n)) -- orthonormal columns
c       ier -- 0 if successful; <0: illegal argument
c
        implicit none
        integer m,n,mn,ier,info,lwork,j,k
        real*8 a(m,n),q(m,*)
        real*8, allocatable :: acopy(:,:),tau(:),work(:)
        real*8 wopt
c
        mn = min(m,n)
        allocate(acopy(m,n),tau(mn))
c
c       Copy a (dgeqrf overwrites its input).
c
        do k = 1,n
          do j = 1,m
            acopy(j,k) = a(j,k)
          enddo
        enddo
c
c       QR factorization (see dqr for parameter explanation).
c
        call dgeqrf(m,n,acopy,m,tau,wopt,-1,info)
        lwork = int(wopt)
        allocate(work(lwork))
        call dgeqrf(m,n,acopy,m,tau,work,lwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: dqreconq failed, info =', info
        if(info .ne. 0) then
          deallocate(acopy,tau,work)
          return
        endif
c
c       Form Q (first mn columns) from Householder reflectors.
c       dorgqr with n_cols_Q = mn gives economy Q.
c
        do k = 1,mn
          do j = 1,m
            q(j,k) = acopy(j,k)
          enddo
        enddo
c
        deallocate(work)
        call dorgqr(m,mn,mn,q,m,tau,wopt,-1,info)
        lwork = int(wopt)
        allocate(work(lwork))
        call dorgqr(m,mn,mn,q,m,tau,work,lwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: dqreconq failed, info =', info
c
        deallocate(acopy,tau,work)
c
        return
        end
c
c
c
c
        subroutine zqr(m,n,a,q,r,ier)
c
c       Full QR factorization of a complex*16 matrix: A = Q * R
c       using LAPACK zgeqrf + zungqr.
c       The input matrix a is not modified.
c
c       input:
c       m -- number of rows of a
c       n -- number of columns of a
c       a(m,n) -- complex*16 matrix
c
c       output:
c       q(m,m) -- unitary matrix
c       r(m,n) -- upper triangular matrix
c       ier -- 0 if successful; <0: illegal argument
c
        implicit none
        integer m,n,mn,ier,info,lwork,j,k
        complex*16 a(m,n),q(m,m),r(m,n)
        complex*16, allocatable :: tau(:),work(:)
        complex*16 wopt
c
        mn = min(m,n)
        allocate(tau(mn))
c
c       Copy a into r (zgeqrf overwrites its input).
c
        do k = 1,n
          do j = 1,m
            r(j,k) = a(j,k)
          enddo
        enddo
c
c       QR factorization.
c       zgeqrf: complex analogue of dgeqrf
c
        call zgeqrf(m,n,r,m,tau,wopt,-1,info)
        lwork = int(dble(wopt))
        allocate(work(lwork))
        call zgeqrf(m,n,r,m,tau,work,lwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: zqr failed, info =', info
        if(info .ne. 0) then
          deallocate(tau,work)
          return
        endif
c
c       Copy R into q, then form Q from Householder reflectors.
c
        do k = 1,m
          do j = 1,m
            q(j,k) = 0
          enddo
        enddo
        do k = 1,mn
          do j = 1,m
            q(j,k) = r(j,k)
          enddo
        enddo
c
c       zungqr: forms unitary Q from Householder reflectors
c       (complex analogue of dorgqr)
c
        deallocate(work)
        call zungqr(m,m,mn,q,m,tau,wopt,-1,info)
        lwork = int(dble(wopt))
        allocate(work(lwork))
        call zungqr(m,m,mn,q,m,tau,work,lwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: zqr failed, info =', info
c
c       Zero out the lower triangle of r.
c
        do k = 1,n
          do j = k+1,m
            r(j,k) = 0
          enddo
        enddo
c
        deallocate(tau,work)
c
        return
        end
c
c
c
c
        subroutine zqrecon(m,n,a,q,r,ier)
c
c       Economy QR factorization of a complex*16 matrix: A = Q * R
c       using LAPACK zgeqrf + zungqr.
c       The input matrix a is not modified.
c
c       input:
c       m -- number of rows of a
c       n -- number of columns of a
c       a(m,n) -- complex*16 matrix
c
c       output:
c       q(m,min(m,n)) -- orthonormal columns
c       r(min(m,n),n) -- upper triangular matrix
c       ier -- 0 if successful; <0: illegal argument
c
        implicit none
        integer m,n,mn,ier,info,lwork,j,k
        complex*16 a(m,n),q(m,*),r(min(m,n),n)
        complex*16, allocatable :: acopy(:,:),tau(:),work(:)
        complex*16 wopt
c
        mn = min(m,n)
        allocate(acopy(m,n),tau(mn))
c
c       Copy a (zgeqrf overwrites its input).
c
        do k = 1,n
          do j = 1,m
            acopy(j,k) = a(j,k)
          enddo
        enddo
c
c       QR factorization (see zqr for parameter explanation).
c
        call zgeqrf(m,n,acopy,m,tau,wopt,-1,info)
        lwork = int(dble(wopt))
        allocate(work(lwork))
        call zgeqrf(m,n,acopy,m,tau,work,lwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: zqrecon failed, info =', info
        if(info .ne. 0) then
          deallocate(acopy,tau,work)
          return
        endif
c
c       Extract R (upper triangle of acopy, first mn rows).
c
        do k = 1,n
          do j = 1,mn
            if(j .le. k) then
              r(j,k) = acopy(j,k)
            else
              r(j,k) = 0
            endif
          enddo
        enddo
c
c       Form Q (first mn columns) from Householder reflectors.
c       zungqr with n_cols_Q = mn gives economy Q.
c
        do k = 1,mn
          do j = 1,m
            q(j,k) = acopy(j,k)
          enddo
        enddo
c
        deallocate(work)
        call zungqr(m,mn,mn,q,m,tau,wopt,-1,info)
        lwork = int(dble(wopt))
        allocate(work(lwork))
        call zungqr(m,mn,mn,q,m,tau,work,lwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: zqrecon failed, info =', info
c
        deallocate(acopy,tau,work)
c
        return
        end
c
c
c
c
        subroutine zqreconq(m,n,a,q,ier)
c
c       Economy QR factorization of a complex*16 matrix, returning
c       only Q: A = Q * R, where Q(m,min(m,n)) has orthonormal
c       columns. R is not computed.
c       The input matrix a is not modified.
c
c       input:
c       m -- number of rows of a
c       n -- number of columns of a
c       a(m,n) -- complex*16 matrix
c
c       output:
c       q(m,min(m,n)) -- orthonormal columns
c       ier -- 0 if successful; <0: illegal argument
c
        implicit none
        integer m,n,mn,ier,info,lwork,j,k
        complex*16 a(m,n),q(m,*)
        complex*16, allocatable :: acopy(:,:),tau(:),work(:)
        complex*16 wopt
c
        mn = min(m,n)
        allocate(acopy(m,n),tau(mn))
c
c       Copy a (zgeqrf overwrites its input).
c
        do k = 1,n
          do j = 1,m
            acopy(j,k) = a(j,k)
          enddo
        enddo
c
c       QR factorization (see zqr for parameter explanation).
c
        call zgeqrf(m,n,acopy,m,tau,wopt,-1,info)
        lwork = int(dble(wopt))
        allocate(work(lwork))
        call zgeqrf(m,n,acopy,m,tau,work,lwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: zqreconq failed, info =', info
        if(info .ne. 0) then
          deallocate(acopy,tau,work)
          return
        endif
c
c       Form Q (first mn columns) from Householder reflectors.
c       zungqr with n_cols_Q = mn gives economy Q.
c
        do k = 1,mn
          do j = 1,m
            q(j,k) = acopy(j,k)
          enddo
        enddo
c
        deallocate(work)
        call zungqr(m,mn,mn,q,m,tau,wopt,-1,info)
        lwork = int(dble(wopt))
        allocate(work(lwork))
        call zungqr(m,mn,mn,q,m,tau,work,lwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: zqreconq failed, info =', info
c
        deallocate(acopy,tau,work)
c
        return
        end
c
c
c======= Column-pivoted QR ====================================
c
c
        subroutine dcpqr(m,n,a,q,r,ipiv,ier)
c
c       Economy column-pivoted QR: A*P = Q*R, using LAPACK dgeqp3.
c       The input matrix a is not modified.
c
c       input:
c       m -- number of rows of a
c       n -- number of columns of a
c       a(m,n) -- matrix
c
c       output:
c       q(m,min(m,n)) -- orthonormal columns
c       r(min(m,n),n) -- upper triangular
c       ipiv(n) -- column permutation: A(:,ipiv) = Q*R
c       ier -- 0 if successful
c
        implicit none
        integer m,n,mn,ier,info,lwork,j,k,ipiv(n)
        real*8 a(m,n),q(m,*),r(min(m,n),n)
        real*8, allocatable :: acopy(:,:),tau(:),work(:)
        real*8 wopt
c
        mn = min(m,n)
        allocate(acopy(m,n),tau(mn))
c
        do k = 1,n
          do j = 1,m
            acopy(j,k) = a(j,k)
          enddo
          ipiv(k) = 0
        enddo
c
c       dgeqp3: column-pivoted QR
c       ipiv=0 on input means free pivoting
        call dgeqp3(m,n,acopy,m,ipiv,tau,wopt,-1,info)
        lwork = int(wopt)
        allocate(work(lwork))
        call dgeqp3(m,n,acopy,m,ipiv,tau,work,lwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: dcpqr failed, info =', info
        if(info .ne. 0) then
          deallocate(acopy,tau,work)
          return
        endif
c
c       Extract R.
c
        do k = 1,n
          do j = 1,mn
            if(j .le. k) then
              r(j,k) = acopy(j,k)
            else
              r(j,k) = 0
            endif
          enddo
        enddo
c
c       Form Q from Householder reflectors.
c
        do k = 1,mn
          do j = 1,m
            q(j,k) = acopy(j,k)
          enddo
        enddo
c
c       dorgqr with n_cols_Q = mn gives economy Q.
        deallocate(work)
        call dorgqr(m,mn,mn,q,m,tau,wopt,-1,info)
        lwork = int(wopt)
        allocate(work(lwork))
        call dorgqr(m,mn,mn,q,m,tau,work,lwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: dcpqr failed, info =', info
c
        deallocate(acopy,tau,work)
c
        return
        end
c
c
c
c
        subroutine dcpqrr(m,n,a,r,ipiv,ier)
c
c       Column-pivoted QR returning only R and pivot: A*P = Q*R.
c       Does not compute Q (faster). Using LAPACK dgeqp3.
c       The input matrix a is not modified.
c
c       input:
c       m -- number of rows of a
c       n -- number of columns of a
c       a(m,n) -- matrix
c
c       output:
c       r(min(m,n),n) -- upper triangular
c       ipiv(n) -- column permutation: A(:,ipiv) = Q*R
c       ier -- 0 if successful
c
        implicit none
        integer m,n,mn,ier,info,lwork,j,k,ipiv(n)
        real*8 a(m,n),r(min(m,n),n)
        real*8, allocatable :: acopy(:,:),tau(:),work(:)
        real*8 wopt
c
        mn = min(m,n)
        allocate(acopy(m,n),tau(mn))
c
        do k = 1,n
          do j = 1,m
            acopy(j,k) = a(j,k)
          enddo
          ipiv(k) = 0
        enddo
c
        call dgeqp3(m,n,acopy,m,ipiv,tau,wopt,-1,info)
        lwork = int(wopt)
        allocate(work(lwork))
        call dgeqp3(m,n,acopy,m,ipiv,tau,work,lwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: dcpqrr failed, info =', info
c
c       Extract R.
c
        do k = 1,n
          do j = 1,mn
            if(j .le. k) then
              r(j,k) = acopy(j,k)
            else
              r(j,k) = 0
            endif
          enddo
        enddo
c
        deallocate(acopy,tau,work)
c
        return
        end
c
c
c
c
        subroutine zcpqr(m,n,a,q,r,ipiv,ier)
c
c       Economy column-pivoted QR: A*P = Q*R, using LAPACK zgeqp3.
c       The input matrix a is not modified.
c
c       input:
c       m -- number of rows of a
c       n -- number of columns of a
c       a(m,n) -- complex*16 matrix
c
c       output:
c       q(m,min(m,n)) -- unitary columns
c       r(min(m,n),n) -- upper triangular
c       ipiv(n) -- column permutation
c       ier -- 0 if successful
c
        implicit none
        integer m,n,mn,ier,info,lwork,j,k,ipiv(n)
        complex*16 a(m,n),q(m,*),r(min(m,n),n)
        complex*16, allocatable :: acopy(:,:),tau(:),work(:)
        real*8, allocatable :: rwork(:)
        complex*16 wopt
c
        mn = min(m,n)
        allocate(acopy(m,n),tau(mn),rwork(2*n))
c
        do k = 1,n
          do j = 1,m
            acopy(j,k) = a(j,k)
          enddo
          ipiv(k) = 0
        enddo
c
c       zgeqp3: complex column-pivoted QR
        call zgeqp3(m,n,acopy,m,ipiv,tau,wopt,-1,rwork,info)
        lwork = int(dble(wopt))
        allocate(work(lwork))
        call zgeqp3(m,n,acopy,m,ipiv,tau,work,lwork,rwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: zcpqr failed, info =', info
        if(info .ne. 0) then
          deallocate(acopy,tau,work,rwork)
          return
        endif
c
c       Extract R.
c
        do k = 1,n
          do j = 1,mn
            if(j .le. k) then
              r(j,k) = acopy(j,k)
            else
              r(j,k) = 0
            endif
          enddo
        enddo
c
c       Form Q.
c       zungqr with n_cols_Q = mn gives economy Q.
c
        do k = 1,mn
          do j = 1,m
            q(j,k) = acopy(j,k)
          enddo
        enddo
c
        deallocate(work)
        call zungqr(m,mn,mn,q,m,tau,wopt,-1,info)
        lwork = int(dble(wopt))
        allocate(work(lwork))
        call zungqr(m,mn,mn,q,m,tau,work,lwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: zcpqr failed, info =', info
c
        deallocate(acopy,tau,work,rwork)
c
        return
        end
c
c
c
c
        subroutine zcpqrr(m,n,a,r,ipiv,ier)
c
c       Column-pivoted QR returning only R and pivot: A*P = Q*R.
c       Does not compute Q (faster). Using LAPACK zgeqp3.
c       The input matrix a is not modified.
c
c       input:
c       m -- number of rows of a
c       n -- number of columns of a
c       a(m,n) -- complex*16 matrix
c
c       output:
c       r(min(m,n),n) -- upper triangular
c       ipiv(n) -- column permutation
c       ier -- 0 if successful
c
        implicit none
        integer m,n,mn,ier,info,lwork,j,k,ipiv(n)
        complex*16 a(m,n),r(min(m,n),n)
        complex*16, allocatable :: acopy(:,:),tau(:),work(:)
        real*8, allocatable :: rwork(:)
        complex*16 wopt
c
        mn = min(m,n)
        allocate(acopy(m,n),tau(mn),rwork(2*n))
c
        do k = 1,n
          do j = 1,m
            acopy(j,k) = a(j,k)
          enddo
          ipiv(k) = 0
        enddo
c
        call zgeqp3(m,n,acopy,m,ipiv,tau,wopt,-1,rwork,info)
        lwork = int(dble(wopt))
        allocate(work(lwork))
        call zgeqp3(m,n,acopy,m,ipiv,tau,work,lwork,rwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: zcpqrr failed, info =', info
c
c       Extract R.
c
        do k = 1,n
          do j = 1,mn
            if(j .le. k) then
              r(j,k) = acopy(j,k)
            else
              r(j,k) = 0
            endif
          enddo
        enddo
c
        deallocate(acopy,tau,work,rwork)
c
        return
        end
c
c
c======= Eigendecomposition ====================================
c
c
        subroutine dsymeig(n,a,w,v,ier)
c
c       Eigendecomposition of a real*8 symmetric matrix: A*v = w*v
c       using LAPACK dsyev. The input matrix a is not modified.
c
c       input:
c       n -- dimension of a (n x n)
c       a(n,n) -- symmetric matrix (only upper triangle used)
c
c       output:
c       w(n) -- eigenvalues in ascending order (real)
c       v(n,n) -- orthonormal eigenvectors (columns)
c       ier -- 0 if successful; >0: convergence failure
c
        implicit none
        integer n,ier,info,lwork,j,k
        real*8 a(n,n),w(n),v(n,n)
        real*8, allocatable :: work(:)
        real*8 wopt
c
c       Copy a into v (dsyev overwrites its input with eigenvectors).
c
        do k = 1,n
          do j = 1,n
            v(j,k) = a(j,k)
          enddo
        enddo
c
c       dsyev: eigenvalues and eigenvectors of symmetric matrix
c       'V' = compute eigenvectors, 'U' = upper triangle stored
        call dsyev('V','U',n,v,n,w,wopt,-1,info)
        lwork = int(wopt)
        allocate(work(lwork))
        call dsyev('V','U',n,v,n,w,work,lwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: dsymeig failed, info =', info
c
        deallocate(work)
c
        return
        end
c
c
c
c
        subroutine zhermeig(n,a,w,v,ier)
c
c       Eigendecomposition of a complex*16 Hermitian matrix: A*v = w*v
c       using LAPACK zheev. The input matrix a is not modified.
c
c       input:
c       n -- dimension of a (n x n)
c       a(n,n) -- Hermitian matrix (only upper triangle used)
c
c       output:
c       w(n) -- eigenvalues in ascending order (real)
c       v(n,n) -- unitary eigenvectors (columns, complex)
c       ier -- 0 if successful; >0: convergence failure
c
        implicit none
        integer n,ier,info,lwork,j,k
        complex*16 a(n,n),v(n,n)
        real*8 w(n)
        complex*16, allocatable :: work(:)
        real*8, allocatable :: rwork(:)
        complex*16 wopt
c
c       Copy a into v (zheev overwrites its input with eigenvectors).
c
        do k = 1,n
          do j = 1,n
            v(j,k) = a(j,k)
          enddo
        enddo
c
c       zheev: eigenvalues and eigenvectors of Hermitian matrix
c       'V' = compute eigenvectors, 'U' = upper triangle stored
        allocate(rwork(3*n-2))
        call zheev('V','U',n,v,n,w,wopt,-1,rwork,info)
        lwork = int(dble(wopt))
        allocate(work(lwork))
        call zheev('V','U',n,v,n,w,work,lwork,rwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: zhermeig failed, info =', info
c
        deallocate(work,rwork)
c
        return
        end
c
c
c
c
        subroutine dgeig(n,a,w,v,ier)
c
c       Eigendecomposition of a real*8 general matrix: A*v = w*v
c       using LAPACK dgeev. Returns complex eigenvalues and
c       eigenvectors. The input matrix a is not modified.
c
c       input:
c       n -- dimension of a (n x n)
c       a(n,n) -- general real matrix
c
c       output:
c       w(n) -- eigenvalues (complex)
c       v(n,n) -- right eigenvectors (complex, columns)
c       ier -- 0 if successful; >0: QR did not converge
c
        implicit none
        integer n,ier,info,lwork,j,k
        real*8 a(n,n)
        complex*16 w(n),v(n,n)
        real*8, allocatable :: acopy(:,:),wr(:),wi(:),vr(:,:)
        real*8, allocatable :: work(:)
        real*8 wopt,vl(1,1)
c
        allocate(acopy(n,n),wr(n),wi(n),vr(n,n))
c
c       Copy a (dgeev destroys its input).
c
        do k = 1,n
          do j = 1,n
            acopy(j,k) = a(j,k)
          enddo
        enddo
c
c       dgeev: eigenvalues and right eigenvectors of general matrix
c       'N' = do not compute left eigenvectors
c       'V' = compute right eigenvectors
c       wr,wi = real and imaginary parts of eigenvalues
c       vr = right eigenvectors (real, with complex pairs encoded)
        call dgeev('N','V',n,acopy,n,wr,wi,vl,1,vr,n,
     1    wopt,-1,info)
        lwork = int(wopt)
        allocate(work(lwork))
        call dgeev('N','V',n,acopy,n,wr,wi,vl,1,vr,n,
     1    work,lwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: dgeig failed, info =', info
c
c       Combine wr,wi into complex w and decode vr into complex v.
c
        if(info .eq. 0) then
          j = 1
 100      if(j .gt. n) goto 200
c
          if(wi(j) .eq. 0) then
c           Real eigenvalue.
            w(j) = dcmplx(wr(j),0.0d0)
            do k = 1,n
              v(k,j) = dcmplx(vr(k,j),0.0d0)
            enddo
            j = j+1
          else
c           Complex conjugate pair: columns j and j+1 of vr
c           encode v_j = vr(:,j) + i*vr(:,j+1)
c           and   v_{j+1} = vr(:,j) - i*vr(:,j+1)
            w(j) = dcmplx(wr(j),wi(j))
            w(j+1) = dcmplx(wr(j+1),wi(j+1))
            do k = 1,n
              v(k,j) = dcmplx(vr(k,j),vr(k,j+1))
              v(k,j+1) = dcmplx(vr(k,j),-vr(k,j+1))
            enddo
            j = j+2
          endif
          goto 100
 200      continue
        endif
c
        deallocate(acopy,wr,wi,vr,work)
c
        return
        end
c
c
c
c
        subroutine zgeig(n,a,w,v,ier)
c
c       Eigendecomposition of a complex*16 general matrix: A*v = w*v
c       using LAPACK zgeev. The input matrix a is not modified.
c
c       input:
c       n -- dimension of a (n x n)
c       a(n,n) -- general complex matrix
c
c       output:
c       w(n) -- eigenvalues (complex)
c       v(n,n) -- right eigenvectors (complex, columns)
c       ier -- 0 if successful; >0: QR did not converge
c
        implicit none
        integer n,ier,info,lwork,j,k
        complex*16 a(n,n),w(n),v(n,n)
        complex*16, allocatable :: acopy(:,:),work(:)
        complex*16 vl(1,1)
        real*8, allocatable :: rwork(:)
        complex*16 wopt
c
        allocate(acopy(n,n),rwork(2*n))
c
c       Copy a (zgeev destroys its input).
c
        do k = 1,n
          do j = 1,n
            acopy(j,k) = a(j,k)
          enddo
        enddo
c
c       zgeev: eigenvalues and right eigenvectors of general matrix
c       'N' = do not compute left eigenvectors
c       'V' = compute right eigenvectors
        call zgeev('N','V',n,acopy,n,w,vl,1,v,n,
     1    wopt,-1,rwork,info)
        lwork = int(dble(wopt))
        allocate(work(lwork))
        call zgeev('N','V',n,acopy,n,w,vl,1,v,n,
     1    work,lwork,rwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: zgeig failed, info =', info
c
        deallocate(acopy,work,rwork)
c
        return
        end
c
c
c======= Linear solve ==========================================
c
c
        subroutine dsolve(n,a,b,x,ier)
c
c       Solves Ax = b for a real*8 n x n matrix using LAPACK dgesv
c       (LU factorization with partial pivoting).
c       The input matrix a and vector b are not modified.
c
c       input:
c       n -- dimension of a (n x n)
c       a(n,n) -- coefficient matrix
c       b(n) -- right-hand side vector
c
c       output:
c       x(n) -- solution vector
c       ier -- 0 if successful; >0: matrix is singular
c              (ier > 0 means singular matrix)
c
        implicit none
        integer n,ier,info,j,k
        real*8 a(n,n),b(n),x(n)
        real*8, allocatable :: acopy(:,:)
        integer, allocatable :: ipiv(:)
c
        allocate(acopy(n,n),ipiv(n))
c
c       Copy a and b (dgesv overwrites both).
c
        do k = 1,n
          do j = 1,n
            acopy(j,k) = a(j,k)
          enddo
        enddo
        do j = 1,n
          x(j) = b(j)
        enddo
c
c       dgesv: solves A*x = b via LU with partial pivoting
c       nrhs=1, overwrites acopy with LU factors, x with solution
        call dgesv(n,1,acopy,n,ipiv,x,n,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: dsolve failed, info =', info
c
        deallocate(acopy,ipiv)
c
        return
        end
c
c
c
c
        subroutine zsolve(n,a,b,x,ier)
c
c       Solves Ax = b for a complex*16 n x n matrix using LAPACK zgesv
c       (LU factorization with partial pivoting).
c       The input matrix a and vector b are not modified.
c
c       input:
c       n -- dimension of a (n x n)
c       a(n,n) -- complex*16 coefficient matrix
c       b(n) -- complex*16 right-hand side vector
c
c       output:
c       x(n) -- complex*16 solution vector
c       ier -- 0 if successful; >0: matrix is singular
c              (ier > 0 means singular matrix)
c
        implicit none
        integer n,ier,info,j,k
        complex*16 a(n,n),b(n),x(n)
        complex*16, allocatable :: acopy(:,:)
        integer, allocatable :: ipiv(:)
c
        allocate(acopy(n,n),ipiv(n))
c
c       Copy a and b (zgesv overwrites both).
c
        do k = 1,n
          do j = 1,n
            acopy(j,k) = a(j,k)
          enddo
        enddo
        do j = 1,n
          x(j) = b(j)
        enddo
c
c       zgesv: solves A*x = b via LU with partial pivoting
c       nrhs=1, overwrites acopy with LU factors, x with solution
        call zgesv(n,1,acopy,n,ipiv,x,n,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: zsolve failed, info =', info
c
        deallocate(acopy,ipiv)
c
        return
        end
c
c
c
c
        subroutine dsvals(m,n,a,s,ier)
c
c       Computes singular values only of a real*8 m x n matrix
c       using LAPACK dgesvd. The input matrix a is not modified.
c
c       input:
c       m -- number of rows of a
c       n -- number of columns of a
c       a(m,n) -- matrix
c
c       output:
c       s(min(m,n)) -- singular values in decreasing order
c       ier -- 0 if successful; >0: SVD did not converge
c
        implicit none
        integer m,n,mn,ier,info,lwork,j,k
        real*8 a(m,n),s(*)
        real*8, allocatable :: acopy(:,:),work(:)
        real*8 wopt
c
        mn = min(m,n)
        allocate(acopy(m,n))
c
        do k = 1,n
          do j = 1,m
            acopy(j,k) = a(j,k)
          enddo
        enddo
c
c       dgesvd: 'N','N' = singular values only, no U or V
        call dgesvd('N','N',m,n,acopy,m,s,0,1,0,1,
     1    wopt,-1,info)
        lwork = int(wopt)
        allocate(work(lwork))
        call dgesvd('N','N',m,n,acopy,m,s,0,1,0,1,
     1    work,lwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: dsvals failed, info =', info
c
        deallocate(acopy,work)
c
        return
        end
c
c
c
c
        subroutine zsvals(m,n,a,s,ier)
c
c       Computes singular values only of a complex*16 m x n matrix
c       using LAPACK zgesvd. The input matrix a is not modified.
c
c       input:
c       m -- number of rows of a
c       n -- number of columns of a
c       a(m,n) -- complex*16 matrix
c
c       output:
c       s(min(m,n)) -- singular values in decreasing order (real)
c       ier -- 0 if successful; >0: SVD did not converge
c
        implicit none
        integer m,n,mn,ier,info,lwork,j,k
        complex*16 a(m,n)
        real*8 s(*)
        complex*16, allocatable :: acopy(:,:),work(:)
        real*8, allocatable :: rwork(:)
        complex*16 wopt
c
        mn = min(m,n)
        allocate(acopy(m,n),rwork(5*mn))
c
        do k = 1,n
          do j = 1,m
            acopy(j,k) = a(j,k)
          enddo
        enddo
c
c       zgesvd: 'N','N' = singular values only, no U or V
        call zgesvd('N','N',m,n,acopy,m,s,0,1,0,1,
     1    wopt,-1,rwork,info)
        lwork = int(dble(wopt))
        allocate(work(lwork))
        call zgesvd('N','N',m,n,acopy,m,s,0,1,0,1,
     1    work,lwork,rwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: zsvals failed, info =', info
c
        deallocate(acopy,work,rwork)
c
        return
        end
c
c
c
c
        subroutine dsymeigvals(n,a,w,ier)
c
c       Eigenvalues only of a real*8 symmetric matrix using
c       LAPACK dsyev. The input matrix a is not modified.
c
c       input:
c       n -- dimension of a (n x n)
c       a(n,n) -- symmetric matrix (only upper triangle used)
c
c       output:
c       w(n) -- eigenvalues in ascending order (real)
c       ier -- 0 if successful; >0: convergence failure
c
        implicit none
        integer n,ier,info,lwork,j,k
        real*8 a(n,n),w(n)
        real*8, allocatable :: acopy(:,:),work(:)
        real*8 wopt
c
        allocate(acopy(n,n))
c
        do k = 1,n
          do j = 1,n
            acopy(j,k) = a(j,k)
          enddo
        enddo
c
c       dsyev: 'N' = eigenvalues only, 'U' = upper triangle
        call dsyev('N','U',n,acopy,n,w,wopt,-1,info)
        lwork = int(wopt)
        allocate(work(lwork))
        call dsyev('N','U',n,acopy,n,w,work,lwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: dsymeigvals failed, info =', info
c
        deallocate(acopy,work)
c
        return
        end
c
c
c
c
        subroutine zhermeigvals(n,a,w,ier)
c
c       Eigenvalues only of a complex*16 Hermitian matrix using
c       LAPACK zheev. The input matrix a is not modified.
c
c       input:
c       n -- dimension of a (n x n)
c       a(n,n) -- Hermitian matrix (only upper triangle used)
c
c       output:
c       w(n) -- eigenvalues in ascending order (real)
c       ier -- 0 if successful; >0: convergence failure
c
        implicit none
        integer n,ier,info,lwork,j,k
        complex*16 a(n,n)
        real*8 w(n)
        complex*16, allocatable :: acopy(:,:),work(:)
        real*8, allocatable :: rwork(:)
        complex*16 wopt
c
        allocate(acopy(n,n),rwork(3*n-2))
c
        do k = 1,n
          do j = 1,n
            acopy(j,k) = a(j,k)
          enddo
        enddo
c
c       zheev: 'N' = eigenvalues only, 'U' = upper triangle
        call zheev('N','U',n,acopy,n,w,wopt,-1,rwork,info)
        lwork = int(dble(wopt))
        allocate(work(lwork))
        call zheev('N','U',n,acopy,n,w,work,lwork,rwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: zhermeigvals failed, info =', info
c
        deallocate(acopy,work,rwork)
c
        return
        end
c
c
c
c
        subroutine dgeigvals(n,a,w,ier)
c
c       Eigenvalues only of a real*8 general matrix using
c       LAPACK dgeev. Returns complex eigenvalues.
c       The input matrix a is not modified.
c
c       input:
c       n -- dimension of a (n x n)
c       a(n,n) -- general real matrix
c
c       output:
c       w(n) -- eigenvalues (complex)
c       ier -- 0 if successful; >0: QR did not converge
c
        implicit none
        integer n,ier,info,lwork,j,k
        real*8 a(n,n)
        complex*16 w(n)
        real*8, allocatable :: acopy(:,:),wr(:),wi(:),work(:)
        real*8 wopt,vl(1,1),vr(1,1)
c
        allocate(acopy(n,n),wr(n),wi(n))
c
        do k = 1,n
          do j = 1,n
            acopy(j,k) = a(j,k)
          enddo
        enddo
c
c       dgeev: 'N','N' = no left or right eigenvectors
        call dgeev('N','N',n,acopy,n,wr,wi,vl,1,vr,1,
     1    wopt,-1,info)
        lwork = int(wopt)
        allocate(work(lwork))
        call dgeev('N','N',n,acopy,n,wr,wi,vl,1,vr,1,
     1    work,lwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: dgeigvals failed, info =', info
c
c       Combine wr,wi into complex w.
c
        if(info .eq. 0) then
          do j = 1,n
            w(j) = dcmplx(wr(j),wi(j))
          enddo
        endif
c
        deallocate(acopy,wr,wi,work)
c
        return
        end
c
c
c
c
        subroutine zgeigvals(n,a,w,ier)
c
c       Eigenvalues only of a complex*16 general matrix using
c       LAPACK zgeev. The input matrix a is not modified.
c
c       input:
c       n -- dimension of a (n x n)
c       a(n,n) -- general complex matrix
c
c       output:
c       w(n) -- eigenvalues (complex)
c       ier -- 0 if successful; >0: QR did not converge
c
        implicit none
        integer n,ier,info,lwork,j,k
        complex*16 a(n,n),w(n)
        complex*16, allocatable :: acopy(:,:),work(:)
        complex*16 vl(1,1),vr(1,1)
        real*8, allocatable :: rwork(:)
        complex*16 wopt
c
        allocate(acopy(n,n),rwork(2*n))
c
        do k = 1,n
          do j = 1,n
            acopy(j,k) = a(j,k)
          enddo
        enddo
c
c       zgeev: 'N','N' = no left or right eigenvectors
        call zgeev('N','N',n,acopy,n,w,vl,1,vr,1,
     1    wopt,-1,rwork,info)
        lwork = int(dble(wopt))
        allocate(work(lwork))
        call zgeev('N','N',n,acopy,n,w,vl,1,vr,1,
     1    work,lwork,rwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: zgeigvals failed, info =', info
c
        deallocate(acopy,work,rwork)
c
        return
        end
c
c
c======= Determinant and inverse ===============================
c
c
        subroutine ddet(n,a,det)
c
c       Computes the determinant of a real*8 n x n matrix via
c       LU factorization. The input matrix a is not modified.
c
c       input:
c       n -- dimension of a
c       a(n,n) -- matrix
c
c       output:
c       det -- determinant of a
c
        implicit none
        integer n,j,k,info
        real*8 a(n,n),det
        real*8, allocatable :: acopy(:,:)
        integer, allocatable :: ipiv(:)
c
        allocate(acopy(n,n),ipiv(n))
c
        do k = 1,n
          do j = 1,n
            acopy(j,k) = a(j,k)
          enddo
        enddo
c
        call dgetrf(n,n,acopy,n,ipiv,info)
c
        if(info .gt. 0) then
c         Singular matrix — determinant is zero.
          det = 0
          deallocate(acopy,ipiv)
          return
        endif
c
c       det = product of diagonal of U, with sign from pivots.
c
        det = 1.0d0
        do j = 1,n
          det = det * acopy(j,j)
          if(ipiv(j) .ne. j) det = -det
        enddo
c
        deallocate(acopy,ipiv)
c
        return
        end
c
c
c
c
        subroutine zdet(n,a,det)
c
c       Computes the determinant of a complex*16 n x n matrix via
c       LU factorization. The input matrix a is not modified.
c
c       input:
c       n -- dimension of a
c       a(n,n) -- complex*16 matrix
c
c       output:
c       det -- determinant of a (complex)
c
        implicit none
        integer n,j,k,info
        complex*16 a(n,n),det
        complex*16, allocatable :: acopy(:,:)
        integer, allocatable :: ipiv(:)
c
        allocate(acopy(n,n),ipiv(n))
c
        do k = 1,n
          do j = 1,n
            acopy(j,k) = a(j,k)
          enddo
        enddo
c
        call zgetrf(n,n,acopy,n,ipiv,info)
c
        if(info .gt. 0) then
c         Singular matrix — determinant is zero.
          det = 0
          deallocate(acopy,ipiv)
          return
        endif
c
c       det = product of diagonal of U, with sign from pivots.
c
        det = 1.0d0
        do j = 1,n
          det = det * acopy(j,j)
          if(ipiv(j) .ne. j) det = -det
        enddo
c
        deallocate(acopy,ipiv)
c
        return
        end
c
c
c
c
        subroutine dinverse(n,a,ainv,ier)
c
c       Computes the inverse of a real*8 n x n matrix using
c       LAPACK dgetrf + dgetri. The input matrix a is not modified.
c
c       input:
c       n -- dimension of a
c       a(n,n) -- matrix to invert
c
c       output:
c       ainv(n,n) -- inverse of a
c       ier -- 0 if successful, >0 if singular
c
        implicit none
        integer n,ier,info,j,k,lwork
        real*8 a(n,n),ainv(n,n)
        integer, allocatable :: ipiv(:)
        real*8, allocatable :: work(:)
        real*8 wopt
c
        allocate(ipiv(n))
c
c       Copy a into ainv.
c
        do k = 1,n
          do j = 1,n
            ainv(j,k) = a(j,k)
          enddo
        enddo
c
c       dgetrf: LU factorize ainv in place.
        call dgetrf(n,n,ainv,n,ipiv,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: dinverse failed, info =', info
        if(info .ne. 0) then
          deallocate(ipiv)
          return
        endif
c
c       dgetri: compute inverse from LU factors.
        call dgetri(n,ainv,n,ipiv,wopt,-1,info)
        lwork = int(wopt)
        allocate(work(lwork))
        call dgetri(n,ainv,n,ipiv,work,lwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: dinverse failed, info =', info
c
        deallocate(ipiv,work)
c
        return
        end
c
c
c
c
        subroutine zinverse(n,a,ainv,ier)
c
c       Computes the inverse of a complex*16 n x n matrix using
c       LAPACK zgetrf + zgetri. The input matrix a is not modified.
c
c       input:
c       n -- dimension of a
c       a(n,n) -- complex*16 matrix to invert
c
c       output:
c       ainv(n,n) -- inverse of a
c       ier -- 0 if successful, >0 if singular
c
        implicit none
        integer n,ier,info,j,k,lwork
        complex*16 a(n,n),ainv(n,n)
        integer, allocatable :: ipiv(:)
        complex*16, allocatable :: work(:)
        complex*16 wopt
c
        allocate(ipiv(n))
c
c       Copy a into ainv.
c
        do k = 1,n
          do j = 1,n
            ainv(j,k) = a(j,k)
          enddo
        enddo
c
c       zgetrf: LU factorize ainv in place.
        call zgetrf(n,n,ainv,n,ipiv,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: zinverse failed, info =', info
        if(info .ne. 0) then
          deallocate(ipiv)
          return
        endif
c
c       zgetri: compute inverse from LU factors.
        call zgetri(n,ainv,n,ipiv,wopt,-1,info)
        lwork = int(dble(wopt))
        allocate(work(lwork))
        call zgetri(n,ainv,n,ipiv,work,lwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: zinverse failed, info =', info
c
        deallocate(ipiv,work)
c
        return
        end
c
c
c======= Least squares =========================================
c
c
        subroutine dlstsq(m,n,a,b,x,ier)
c
c       Least squares solve using LAPACK dgels (QR-based).
c       Assumes a has full rank.
c       If m >= n (overdetermined): minimizes ||Ax - b||_2.
c       If m < n (underdetermined): minimum norm solution.
c       The input matrix a and vector b are not modified.
c
c       input:
c       m -- number of rows of a
c       n -- number of columns of a
c       a(m,n) -- matrix
c       b(m) -- right-hand side vector
c
c       output:
c       x(n) -- solution vector
c       ier -- 0 if successful; >0: not full rank
c
        implicit none
        integer m,n,mn,ier,info,lwork,j,k
        real*8 a(m,n),b(m),x(n)
        real*8, allocatable :: acopy(:,:),bx(:),work(:)
        real*8 wopt
c
        mn = max(m,n)
        allocate(acopy(m,n),bx(mn))
c
c       Copy a (dgels overwrites it).
c
        do k = 1,n
          do j = 1,m
            acopy(j,k) = a(j,k)
          enddo
        enddo
c
c       Copy b into bx (dgels needs length max(m,n)).
c       Zero-pad if m < n.
c
        do j = 1,m
          bx(j) = b(j)
        enddo
        do j = m+1,mn
          bx(j) = 0
        enddo
c
c       dgels: solve least squares via QR or LQ
c       'N' = no transpose
        call dgels('N',m,n,1,acopy,m,bx,mn,wopt,-1,info)
        lwork = int(wopt)
        allocate(work(lwork))
        call dgels('N',m,n,1,acopy,m,bx,mn,work,lwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: dlstsq failed, info =', info
c
c       Extract solution (first n elements of bx).
c
        do j = 1,n
          x(j) = bx(j)
        enddo
c
        deallocate(acopy,bx,work)
c
        return
        end
c
c
c
c
        subroutine zlstsq(m,n,a,b,x,ier)
c
c       Least squares solve using LAPACK zgels (QR-based).
c       Assumes a has full rank.
c       If m >= n (overdetermined): minimizes ||Ax - b||_2.
c       If m < n (underdetermined): minimum norm solution.
c       The input matrix a and vector b are not modified.
c
c       input:
c       m -- number of rows of a
c       n -- number of columns of a
c       a(m,n) -- complex*16 matrix
c       b(m) -- complex*16 right-hand side vector
c
c       output:
c       x(n) -- complex*16 solution vector
c       ier -- 0 if successful; >0: not full rank
c
        implicit none
        integer m,n,mn,ier,info,lwork,j,k
        complex*16 a(m,n),b(m),x(n)
        complex*16, allocatable :: acopy(:,:),bx(:),work(:)
        complex*16 wopt
c
        mn = max(m,n)
        allocate(acopy(m,n),bx(mn))
c
c       Copy a (zgels overwrites it).
c
        do k = 1,n
          do j = 1,m
            acopy(j,k) = a(j,k)
          enddo
        enddo
c
c       Copy b into bx (zgels needs length max(m,n)).
c
        do j = 1,m
          bx(j) = b(j)
        enddo
        do j = m+1,mn
          bx(j) = 0
        enddo
c
c       zgels: solve least squares via QR or LQ
c       'N' = no conjugate transpose
        call zgels('N',m,n,1,acopy,m,bx,mn,wopt,-1,info)
        lwork = int(dble(wopt))
        allocate(work(lwork))
        call zgels('N',m,n,1,acopy,m,bx,mn,work,lwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: zlstsq failed, info =', info
c
c       Extract solution (first n elements of bx).
c
        do j = 1,n
          x(j) = bx(j)
        enddo
c
        deallocate(acopy,bx,work)
c
        return
        end
c
c
c
c
        subroutine dlstsqr(m,n,a,b,x,krank,ier)
c
c       Robust least squares solve using LAPACK dgelsd (SVD-based).
c       Handles rank-deficient matrices.
c       If m >= n: minimizes ||Ax - b||_2.
c       If m < n: minimum norm solution.
c       The input matrix a and vector b are not modified.
c
c       input:
c       m -- number of rows of a
c       n -- number of columns of a
c       a(m,n) -- matrix
c       b(m) -- right-hand side vector
c
c       output:
c       x(n) -- solution vector
c       krank -- effective rank of a
c       ier -- 0 if successful; >0: SVD did not converge
c
        implicit none
        integer m,n,mn,krank,ier,info,lwork,j,k
        real*8 a(m,n),b(m),x(n),rcond
        real*8, allocatable :: acopy(:,:),bx(:),s(:),work(:)
        integer, allocatable :: iwork(:)
        real*8 wopt
c
        mn = max(m,n)
        allocate(acopy(m,n),bx(mn),s(min(m,n)))
        allocate(iwork(3*min(m,n)*min(11,2+int(log(dble(min(m,n)))
     1    /log(2d0)))+11*min(m,n)))
c
c       Copy a and b.
c
        do k = 1,n
          do j = 1,m
            acopy(j,k) = a(j,k)
          enddo
        enddo
        do j = 1,m
          bx(j) = b(j)
        enddo
        do j = m+1,mn
          bx(j) = 0
        enddo
c
c       rcond < 0 means use machine precision for rank determination.
        rcond = -1.0d0
c
c       dgelsd: least squares via SVD (divide-and-conquer)
        call dgelsd(m,n,1,acopy,m,bx,mn,s,rcond,krank,
     1    wopt,-1,iwork,info)
        lwork = int(wopt)
        allocate(work(lwork))
        call dgelsd(m,n,1,acopy,m,bx,mn,s,rcond,krank,
     1    work,lwork,iwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: dlstsqr failed, info =', info
c
        do j = 1,n
          x(j) = bx(j)
        enddo
c
        deallocate(acopy,bx,s,work,iwork)
c
        return
        end
c
c
c
c
        subroutine zlstsqr(m,n,a,b,x,krank,ier)
c
c       Robust least squares solve using LAPACK zgelsd (SVD-based).
c       Handles rank-deficient matrices.
c       If m >= n: minimizes ||Ax - b||_2.
c       If m < n: minimum norm solution.
c       The input matrix a and vector b are not modified.
c
c       input:
c       m -- number of rows of a
c       n -- number of columns of a
c       a(m,n) -- complex*16 matrix
c       b(m) -- complex*16 right-hand side vector
c
c       output:
c       x(n) -- complex*16 solution vector
c       krank -- effective rank of a
c       ier -- 0 if successful; >0: SVD did not converge
c
        implicit none
        integer m,n,mn,mnn,krank,ier,info,lwork,lrwork,j,k
        real*8 rcond
        complex*16 a(m,n),b(m),x(n)
        complex*16, allocatable :: acopy(:,:),bx(:),work(:)
        real*8, allocatable :: s(:),rwork(:)
        integer, allocatable :: iwork(:)
        complex*16 wopt
        real*8 rwopt
c
        mn = max(m,n)
        mnn = min(m,n)
        allocate(acopy(m,n),bx(mn),s(mnn))
        allocate(iwork(3*mnn*min(11,2+int(log(dble(mnn))
     1    /log(2d0)))+11*mnn))
c
c       Copy a and b.
c
        do k = 1,n
          do j = 1,m
            acopy(j,k) = a(j,k)
          enddo
        enddo
        do j = 1,m
          bx(j) = b(j)
        enddo
        do j = m+1,mn
          bx(j) = 0
        enddo
c
c       rcond < 0 means use machine precision for rank determination.
        rcond = -1.0d0
c
c       zgelsd: least squares via SVD (divide-and-conquer)
c       Workspace query.
        lrwork = 1
        allocate(rwork(lrwork))
        call zgelsd(m,n,1,acopy,m,bx,mn,s,rcond,krank,
     1    wopt,-1,rwork,iwork,info)
        lwork = int(dble(wopt))
        lrwork = int(rwork(1))
        deallocate(rwork)
        allocate(work(lwork),rwork(lrwork))
c
c       Re-copy since workspace query may have touched acopy.
        do k = 1,n
          do j = 1,m
            acopy(j,k) = a(j,k)
          enddo
        enddo
        do j = 1,m
          bx(j) = b(j)
        enddo
        do j = m+1,mn
          bx(j) = 0
        enddo
c
        call zgelsd(m,n,1,acopy,m,bx,mn,s,rcond,krank,
     1    work,lwork,rwork,iwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: zlstsqr failed, info =', info
c
        do j = 1,n
          x(j) = bx(j)
        enddo
c
        deallocate(acopy,bx,s,work,rwork,iwork)
c
        return
        end
c
c
c
c
        subroutine dqrpack(m,n,a,qr,tau,ier)
c
c       Precomputes QR factorization of a real*8 m x n matrix
c       (m >= n) using LAPACK dgeqrf. The packed format stores R
c       in the upper triangle and Householder reflectors below.
c       The input matrix a is not modified.
c
c       input:
c       m -- number of rows (m >= n)
c       n -- number of columns
c       a(m,n) -- matrix
c
c       output:
c       qr(m,n) -- packed QR factorization
c       tau(n) -- Householder scalar factors
c       ier -- 0 if successful
c
        implicit none
        integer m,n,ier,info,lwork,j,k
        real*8 a(m,n),qr(m,n),tau(n)
        real*8, allocatable :: work(:)
        real*8 wopt
c
        do k = 1,n
          do j = 1,m
            qr(j,k) = a(j,k)
          enddo
        enddo
c
        call dgeqrf(m,n,qr,m,tau,wopt,-1,info)
        lwork = int(wopt)
        allocate(work(lwork))
        call dgeqrf(m,n,qr,m,tau,work,lwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: dqrpack failed, info =', info
c
        deallocate(work)
c
        return
        end
c
c
c
c
        subroutine zqrpack(m,n,a,qr,tau,ier)
c
c       Precomputes QR factorization of a complex*16 m x n matrix
c       (m >= n) using LAPACK zgeqrf. The packed format stores R
c       in the upper triangle and Householder reflectors below.
c       The input matrix a is not modified.
c
c       input:
c       m -- number of rows (m >= n)
c       n -- number of columns
c       a(m,n) -- complex*16 matrix
c
c       output:
c       qr(m,n) -- packed QR factorization (complex*16)
c       tau(n) -- Householder scalar factors (complex*16)
c       ier -- 0 if successful
c
        implicit none
        integer m,n,ier,info,lwork,j,k
        complex*16 a(m,n),qr(m,n),tau(n)
        complex*16, allocatable :: work(:)
        complex*16 wopt
c
        do k = 1,n
          do j = 1,m
            qr(j,k) = a(j,k)
          enddo
        enddo
c
        call zgeqrf(m,n,qr,m,tau,wopt,-1,info)
        lwork = int(dble(wopt))
        allocate(work(lwork))
        call zgeqrf(m,n,qr,m,tau,work,lwork,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: zqrpack failed, info =', info
c
        deallocate(work)
c
        return
        end
c
c
c
c
        subroutine dqrsolve(m,n,qr,tau,b,x,ier)
c
c       Solves the least squares problem min ||Ax - b|| using a
c       precomputed QR factorization (from dqrpack).
c       Requires m >= n. Also works as a linear solve when m = n.
c       The inputs qr, tau, and b are not modified.
c
c       input:
c       m -- number of rows (m >= n)
c       n -- number of columns
c       qr(m,n) -- packed QR factorization from dqrpack
c       tau(n) -- Householder scalar factors from dqrpack
c       b(m) -- right-hand side vector
c
c       output:
c       x(n) -- solution vector
c       ier -- 0 if successful
c
        implicit none
        integer m,n,ier,info,lwork,j
        real*8 qr(m,n),tau(n),b(m),x(n)
        real*8, allocatable :: bx(:),work(:)
        real*8 wopt
c
c       Copy b (dormqr overwrites it).
c
        allocate(bx(m))
        do j = 1,m
          bx(j) = b(j)
        enddo
c
c       Apply Q^T to b: bx <- Q^T * b
c
        call dormqr('L','T',m,1,n,qr,m,tau,bx,m,wopt,-1,info)
        lwork = int(wopt)
        allocate(work(lwork))
        call dormqr('L','T',m,1,n,qr,m,tau,bx,m,work,lwork,info)
        if(info .ne. 0) then
          print *, 'linalgwrap: dqrsolve dormqr failed, info =',info
          ier = info
          deallocate(bx,work)
          return
        endif
c
c       Solve R*x = (Q^T*b): first n entries of bx
c
        call dtrtrs('U','N','N',n,1,qr,m,bx,m,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: dqrsolve dtrtrs failed, info =',info
c
        do j = 1,n
          x(j) = bx(j)
        enddo
c
        deallocate(bx,work)
c
        return
        end
c
c
c
c
        subroutine zqrsolve(m,n,qr,tau,b,x,ier)
c
c       Solves the least squares problem min ||Ax - b|| using a
c       precomputed QR factorization (from zqrpack).
c       Requires m >= n. Also works as a linear solve when m = n.
c       The inputs qr, tau, and b are not modified.
c
c       input:
c       m -- number of rows (m >= n)
c       n -- number of columns
c       qr(m,n) -- packed QR factorization from zqrpack (complex*16)
c       tau(n) -- Householder scalar factors from zqrpack (complex*16)
c       b(m) -- right-hand side vector (complex*16)
c
c       output:
c       x(n) -- solution vector (complex*16)
c       ier -- 0 if successful
c
        implicit none
        integer m,n,ier,info,lwork,j
        complex*16 qr(m,n),tau(n),b(m),x(n)
        complex*16, allocatable :: bx(:),work(:)
        complex*16 wopt
c
c       Copy b (zunmqr overwrites it).
c
        allocate(bx(m))
        do j = 1,m
          bx(j) = b(j)
        enddo
c
c       Apply Q^H to b: bx <- Q^H * b
c
        call zunmqr('L','C',m,1,n,qr,m,tau,bx,m,wopt,-1,info)
        lwork = int(dble(wopt))
        allocate(work(lwork))
        call zunmqr('L','C',m,1,n,qr,m,tau,bx,m,work,lwork,info)
        if(info .ne. 0) then
          print *, 'linalgwrap: zqrsolve zunmqr failed, info =',info
          ier = info
          deallocate(bx,work)
          return
        endif
c
c       Solve R*x = (Q^H*b): first n entries of bx
c
        call ztrtrs('U','N','N',n,1,qr,m,bx,m,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: zqrsolve ztrtrs failed, info =',info
c
        do j = 1,n
          x(j) = bx(j)
        enddo
c
        deallocate(bx,work)
c
        return
        end
c
c
c======= LU factorization =====================================
c
c
        subroutine dlu(n,a,lu,ipiv,ier)
c
c       LU factorization of a real*8 n x n matrix: P*A = L*U
c       using LAPACK dgetrf. The input matrix a is not modified.
c
c       input:
c       n -- dimension of a
c       a(n,n) -- matrix to factorize
c
c       output:
c       lu(n,n) -- combined L and U factors (LAPACK packed format)
c       ipiv(n) -- pivot indices
c       ier -- 0 if successful, >0 if singular
c
        implicit none
        integer n,ier,info,j,k,ipiv(n)
        real*8 a(n,n),lu(n,n)
c
        do k = 1,n
          do j = 1,n
            lu(j,k) = a(j,k)
          enddo
        enddo
c
c       dgetrf: A = P*L*U via partial pivoting
        call dgetrf(n,n,lu,n,ipiv,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: dlu failed, info =', info
c
        return
        end
c
c
c
c
        subroutine dlui(n,a,ipiv,ier)
c
c       In-place LU factorization of a real*8 n x n matrix.
c       Overwrites a with the combined L and U factors.
c
c       input:
c       n -- dimension of a
c       a(n,n) -- matrix to factorize
c
c       output:
c       a(n,n) -- combined L and U factors (LAPACK packed format)
c       ipiv(n) -- pivot indices
c       ier -- 0 if successful, >0 if singular
c
        implicit none
        integer n,ier,info,ipiv(n)
        real*8 a(n,n)
c
        call dgetrf(n,n,a,n,ipiv,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: dlui failed, info =', info
c
        return
        end
c
c
c
c
        subroutine dlusolve(n,lu,ipiv,b,x)
c
c       Solves Ax = b using an existing LU factorization (from dlu
c       or dlui). The vector b is not modified.
c
c       input:
c       n -- dimension
c       lu(n,n) -- LU factors from dlu or dlui
c       ipiv(n) -- pivot indices from dlu or dlui
c       b(n) -- right-hand side vector
c
c       output:
c       x(n) -- solution vector
c
        implicit none
        integer n,j,ipiv(n),info
        real*8 lu(n,n),b(n),x(n)
c
        do j = 1,n
          x(j) = b(j)
        enddo
c
c       dgetrs: solve using LU factors, 'N' = no transpose
        call dgetrs('N',n,1,lu,n,ipiv,x,n,info)
c
        return
        end
c
c
c
c
        subroutine dlusolven(n,nrhs,lu,ipiv,b,x)
c
c       Solves AX = B for multiple right-hand sides using an
c       existing LU factorization. B is not modified.
c
c       input:
c       n -- dimension
c       nrhs -- number of right-hand sides
c       lu(n,n) -- LU factors from dlu or dlui
c       ipiv(n) -- pivot indices from dlu or dlui
c       b(n,nrhs) -- right-hand side matrix
c
c       output:
c       x(n,nrhs) -- solution matrix
c
        implicit none
        integer n,nrhs,j,k,ipiv(n),info
        real*8 lu(n,n),b(n,nrhs),x(n,nrhs)
c
        do k = 1,nrhs
          do j = 1,n
            x(j,k) = b(j,k)
          enddo
        enddo
c
c       dgetrs: solve using LU factors, 'N' = no transpose
        call dgetrs('N',n,nrhs,lu,n,ipiv,x,n,info)
c
        return
        end
c
c
c
c
        subroutine zlu(n,a,lu,ipiv,ier)
c
c       LU factorization of a complex*16 n x n matrix: P*A = L*U
c       using LAPACK zgetrf. The input matrix a is not modified.
c
c       input:
c       n -- dimension of a
c       a(n,n) -- complex*16 matrix to factorize
c
c       output:
c       lu(n,n) -- combined L and U factors (LAPACK packed format)
c       ipiv(n) -- pivot indices
c       ier -- 0 if successful, >0 if singular
c
        implicit none
        integer n,ier,info,j,k,ipiv(n)
        complex*16 a(n,n),lu(n,n)
c
        do k = 1,n
          do j = 1,n
            lu(j,k) = a(j,k)
          enddo
        enddo
c
c       zgetrf: A = P*L*U via partial pivoting
        call zgetrf(n,n,lu,n,ipiv,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: zlu failed, info =', info
c
        return
        end
c
c
c
c
        subroutine zlui(n,a,ipiv,ier)
c
c       In-place LU factorization of a complex*16 n x n matrix.
c       Overwrites a with the combined L and U factors.
c
c       input:
c       n -- dimension of a
c       a(n,n) -- complex*16 matrix to factorize
c
c       output:
c       a(n,n) -- combined L and U factors (LAPACK packed format)
c       ipiv(n) -- pivot indices
c       ier -- 0 if successful, >0 if singular
c
        implicit none
        integer n,ier,info,ipiv(n)
        complex*16 a(n,n)
c
        call zgetrf(n,n,a,n,ipiv,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: zlui failed, info =', info
c
        return
        end
c
c
c
c
        subroutine zlusolve(n,lu,ipiv,b,x)
c
c       Solves Ax = b using an existing complex*16 LU factorization.
c       The vector b is not modified.
c
c       input:
c       n -- dimension
c       lu(n,n) -- complex*16 LU factors from zlu or zlui
c       ipiv(n) -- pivot indices from zlu or zlui
c       b(n) -- complex*16 right-hand side vector
c
c       output:
c       x(n) -- complex*16 solution vector
c
        implicit none
        integer n,j,ipiv(n),info
        complex*16 lu(n,n),b(n),x(n)
c
        do j = 1,n
          x(j) = b(j)
        enddo
c
c       zgetrs: solve using LU factors, 'N' = no transpose
        call zgetrs('N',n,1,lu,n,ipiv,x,n,info)
c
        return
        end
c
c
c
c
        subroutine zlusolven(n,nrhs,lu,ipiv,b,x)
c
c       Solves AX = B for multiple right-hand sides using an
c       existing complex*16 LU factorization. B is not modified.
c
c       input:
c       n -- dimension
c       nrhs -- number of right-hand sides
c       lu(n,n) -- complex*16 LU factors from zlu or zlui
c       ipiv(n) -- pivot indices from zlu or zlui
c       b(n,nrhs) -- complex*16 right-hand side matrix
c
c       output:
c       x(n,nrhs) -- complex*16 solution matrix
c
        implicit none
        integer n,nrhs,j,k,ipiv(n),info
        complex*16 lu(n,n),b(n,nrhs),x(n,nrhs)
c
        do k = 1,nrhs
          do j = 1,n
            x(j,k) = b(j,k)
          enddo
        enddo
c
c       zgetrs: solve using LU factors, 'N' = no transpose
        call zgetrs('N',n,nrhs,lu,n,ipiv,x,n,info)
c
        return
        end
c
c
c======= Cholesky factorization (SPD/HPD) =====================
c
c
        subroutine dchol(n,a,l,ier)
c
c       Cholesky factorization of a real*8 symmetric positive definite
c       matrix: A = L * L^T, using LAPACK dpotrf.
c       The input matrix a is not modified.
c
c       input:
c       n -- dimension of a
c       a(n,n) -- SPD matrix
c
c       output:
c       l(n,n) -- lower triangular Cholesky factor
c       ier -- 0 if successful, >0 if not SPD
c
        implicit none
        integer n,ier,info,j,k
        real*8 a(n,n),l(n,n)
c
        do k = 1,n
          do j = 1,n
            l(j,k) = a(j,k)
          enddo
        enddo
c
c       dpotrf: 'L' = compute lower triangular factor
        call dpotrf('L',n,l,n,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: dchol failed, info =', info
c
c       Zero upper triangle.
c
        do k = 2,n
          do j = 1,k-1
            l(j,k) = 0
          enddo
        enddo
c
        return
        end
c
c
c
c
        subroutine dcholi(n,a,ier)
c
c       In-place Cholesky factorization of a real*8 SPD matrix.
c       Overwrites a with the lower triangular factor L.
c
c       input:
c       n -- dimension of a
c       a(n,n) -- SPD matrix
c
c       output:
c       a(n,n) -- lower triangular Cholesky factor
c       ier -- 0 if successful, >0 if not SPD
c
        implicit none
        integer n,ier,info,j,k
        real*8 a(n,n)
c
        call dpotrf('L',n,a,n,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: dcholi failed, info =', info
c
c       Zero upper triangle.
c
        do k = 2,n
          do j = 1,k-1
            a(j,k) = 0
          enddo
        enddo
c
        return
        end
c
c
c
c
        subroutine dcholsolve(n,l,b,x)
c
c       Solves Ax = b using an existing Cholesky factorization
c       A = L*L^T. The vector b is not modified.
c
c       input:
c       n -- dimension
c       l(n,n) -- lower triangular Cholesky factor from dchol/dcholi
c       b(n) -- right-hand side vector
c
c       output:
c       x(n) -- solution vector
c
        implicit none
        integer n,j,info
        real*8 l(n,n),b(n),x(n)
c
        do j = 1,n
          x(j) = b(j)
        enddo
c
c       dpotrs: solve using Cholesky factor, 'L' = lower triangular
        call dpotrs('L',n,1,l,n,x,n,info)
c
        return
        end
c
c
c
c
        subroutine dcholsolven(n,nrhs,l,b,x)
c
c       Solves AX = B for multiple RHS using an existing Cholesky
c       factorization. B is not modified.
c
c       input:
c       n -- dimension
c       nrhs -- number of right-hand sides
c       l(n,n) -- lower triangular Cholesky factor
c       b(n,nrhs) -- right-hand side matrix
c
c       output:
c       x(n,nrhs) -- solution matrix
c
        implicit none
        integer n,nrhs,j,k,info
        real*8 l(n,n),b(n,nrhs),x(n,nrhs)
c
        do k = 1,nrhs
          do j = 1,n
            x(j,k) = b(j,k)
          enddo
        enddo
c
        call dpotrs('L',n,nrhs,l,n,x,n,info)
c
        return
        end
c
c
c
c
        subroutine zchol(n,a,l,ier)
c
c       Cholesky factorization of a complex*16 Hermitian positive
c       definite matrix: A = L * L^*, using LAPACK zpotrf.
c       The input matrix a is not modified.
c
c       input:
c       n -- dimension of a
c       a(n,n) -- HPD matrix
c
c       output:
c       l(n,n) -- lower triangular Cholesky factor
c       ier -- 0 if successful, >0 if not HPD
c
        implicit none
        integer n,ier,info,j,k
        complex*16 a(n,n),l(n,n)
c
        do k = 1,n
          do j = 1,n
            l(j,k) = a(j,k)
          enddo
        enddo
c
c       zpotrf: 'L' = compute lower triangular factor
        call zpotrf('L',n,l,n,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: zchol failed, info =', info
c
c       Zero upper triangle.
c
        do k = 2,n
          do j = 1,k-1
            l(j,k) = 0
          enddo
        enddo
c
        return
        end
c
c
c
c
        subroutine zcholi(n,a,ier)
c
c       In-place Cholesky factorization of a complex*16 HPD matrix.
c       Overwrites a with the lower triangular factor L.
c
c       input:
c       n -- dimension of a
c       a(n,n) -- HPD matrix
c
c       output:
c       a(n,n) -- lower triangular Cholesky factor
c       ier -- 0 if successful, >0 if not HPD
c
        implicit none
        integer n,ier,info,j,k
        complex*16 a(n,n)
c
        call zpotrf('L',n,a,n,info)
        ier = info
        if(info .ne. 0)
     1    print *, 'linalgwrap: zcholi failed, info =', info
c
c       Zero upper triangle.
c
        do k = 2,n
          do j = 1,k-1
            a(j,k) = 0
          enddo
        enddo
c
        return
        end
c
c
c
c
        subroutine zcholsolve(n,l,b,x)
c
c       Solves Ax = b using an existing Cholesky factorization
c       A = L*L^*. The vector b is not modified.
c
c       input:
c       n -- dimension
c       l(n,n) -- lower triangular Cholesky factor from zchol/zcholi
c       b(n) -- complex*16 right-hand side vector
c
c       output:
c       x(n) -- complex*16 solution vector
c
        implicit none
        integer n,j,info
        complex*16 l(n,n),b(n),x(n)
c
        do j = 1,n
          x(j) = b(j)
        enddo
c
c       zpotrs: solve using Cholesky factor, 'L' = lower triangular
        call zpotrs('L',n,1,l,n,x,n,info)
c
        return
        end
c
c
c
c
        subroutine zcholsolven(n,nrhs,l,b,x)
c
c       Solves AX = B for multiple RHS using an existing complex*16
c       Cholesky factorization. B is not modified.
c
c       input:
c       n -- dimension
c       nrhs -- number of right-hand sides
c       l(n,n) -- complex*16 lower triangular Cholesky factor
c       b(n,nrhs) -- complex*16 right-hand side matrix
c
c       output:
c       x(n,nrhs) -- complex*16 solution matrix
c
        implicit none
        integer n,nrhs,j,k,info
        complex*16 l(n,n),b(n,nrhs),x(n,nrhs)
c
        do k = 1,nrhs
          do j = 1,n
            x(j,k) = b(j,k)
          enddo
        enddo
c
        call zpotrs('L',n,nrhs,l,n,x,n,info)
c
        return
        end
c
c
c======= Triangular solve =====================================
c
c
        subroutine dtrisolve(n,t,b,x,iuplo)
c
c       Solves T*x = b where T is a real*8 triangular matrix,
c       using BLAS dtrsv. The vector b is not modified.
c
c       input:
c       n -- dimension
c       t(n,n) -- triangular matrix
c       b(n) -- right-hand side vector
c       iuplo -- 0 for lower triangular, 1 for upper triangular
c
c       output:
c       x(n) -- solution vector
c
        implicit none
        integer n,j,iuplo
        character*1 uplo
        real*8 t(n,n),b(n),x(n)
c
        if(iuplo .eq. 0) then
          uplo = 'L'
        else
          uplo = 'U'
        endif
c
        do j = 1,n
          x(j) = b(j)
        enddo
c
c       dtrsv: solve T*x = b
c       uplo = 'L' or 'U', 'N' = no transpose, 'N' = non-unit diagonal
        call dtrsv(uplo,'N','N',n,t,n,x,1)
c
        return
        end
c
c
c
c
        subroutine ztrisolve(n,t,b,x,iuplo)
c
c       Solves T*x = b where T is a complex*16 triangular matrix,
c       using BLAS ztrsv. The vector b is not modified.
c
c       input:
c       n -- dimension
c       t(n,n) -- complex*16 triangular matrix
c       b(n) -- complex*16 right-hand side vector
c       iuplo -- 0 for lower triangular, 1 for upper triangular
c
c       output:
c       x(n) -- complex*16 solution vector
c
        implicit none
        integer n,j,iuplo
        character*1 uplo
        complex*16 t(n,n),b(n),x(n)
c
        if(iuplo .eq. 0) then
          uplo = 'L'
        else
          uplo = 'U'
        endif
c
        do j = 1,n
          x(j) = b(j)
        enddo
c
c       ztrsv: solve T*x = b
c       uplo = 'L' or 'U', 'N' = no transpose, 'N' = non-unit diagonal
        call ztrsv(uplo,'N','N',n,t,n,x,1)
c
        return
        end
