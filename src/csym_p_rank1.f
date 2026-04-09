cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     QR algorithm for complex symmetric tridiagonal + rank-1 update.
c
c     Computes the eigenvalues of
c
c              C = A + pv * qv^T
c
c     where A is complex symmetric tridiagonal (diagonal vd,
c     superdiagonal voffd) and pv, qv are vectors.
c     Uses complex orthogonal (not unitary) Givens rotations
c     to preserve the structure.  Cost is O(n^2).
c
c     Dependencies: none
c
c     Entry point:
c       sparse_root_struct  - sets up exiteps, calls cqr_eigenvals_wilk
c
c     Internal subroutines:
c       cqr_eigenvals_wilk        - shifted QR with deflation stack
c       qr_one_step               - one Wilkinson-shifted QR step
c       cqr_one_iter_dumb         - one O(n) QR sweep
c       find_rot_dumb             - complex orthogonal 2x2 rotation
c       eigenvals_quadr_sol2_elem - eigenvalues of 2x2 matrix
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine sparse_root_struct(croot,vd,voffd,norder,coefs,
     1   ierqr)
        implicit real *8 (a-h,o-z)
        complex *16 croot(norder)
        complex *16 vd(*),voffd(*),coefs(*)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   This subroutine computes the roots of a polynomial
c                     norder
c           p(z)= \Sum       coefs_(i+1) P_i            (1)
c                     i=0
c   with a user specified coefs.
c
c   The vectors vd, voffd specify the three-term recurrence  that define
c   polynomials {P_i}_{i=0}^norder.
c
c   It does so by finding the eigenvalues of the so-called
c   colleague matrix, constructed from the expansion coefficients.
c   This subroutine does not assume {P_i} to classical polynomials;
c   it only assume they satisfy a three-term recurrence, so
c   the colleague matrix is in general complex symmetric, as opposed
c   to being real symmetric. As a result, complex orthogonal transforms
c   are used to diagonalize the colleague matrix.
c
c   The claim to fame of this subroutine is that it has cost O(n^2).
c   Besides, numerical experiments strongly suggest that it is
c   also componenetwise backward stable in the sense that the computed
c   eigenvalues are the exact eigenvalues of a perturbed matrix, where
c   the perturbation is proportional to the size of each entry * eps.
c
c
c                          WARNINGS:
c   1. Detailed perturbation analysis has not been done.
c   2. This subroutine uses complex orthogonal transforms for
c      diagonalizing the colleague matrix, so the numerical stability
c      is not guaranteed!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                               inputs:
c   norder: order of the polynomial expansion
c   vd: diagonal of the colleague matrix, contains norder elements,
c       not touched by the subroutine
c   voffd: superdiagonal of the colleague matrix, contains norder
c          elements not touched by the subroutine
c   coefs: polynomial expansion coefficients, contains norder elements,
c          not touched by the subroutine
c
c                               outputs:
c   croot: computed eigenvalues of the colleague matrix, thus
c          the roots of the polynomial expansion (1), contains
c          norder elements
c   ierqr:
c        0 qr iteration successful
c        1028 qr fails to converge to acc. set by exiteps

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        exiteps=epsilon(1d0)

        csize=0
        do i=1,norder-1
          csize=abs(vd(i))**2 + 2*abs(voffd(i))**2 + csize
        enddo

        csize=csize+abs(vd(norder))**2

        csize=sqrt(csize)

        exiteps=exiteps*csize

        call cqr_eigenvals_wilk(croot,norder,exiteps,ierqr,
     1              vd,voffd,coefs)


        return
        end



        subroutine cqr_one_iter_dumb(n,vd,voffd,pv,qv)
        implicit real *8 (a-h,o-z)
        complex *16 vd(*),voffd(*),qv(*),pv(*)
        complex *16 cpqv(10 000),gv(10 000)
        complex *16 sprot(2,n),c,s,cc,cs
        complex *16 d1,d2

c
c sparse qr starts here
c
        do i=1,n
          cpqv(i)=qv(i)
          gv(i)=voffd(i)
        enddo

        do k=n,2,-1

          d1=voffd(k-1)+pv(k-1)*qv(k)
          d2=vd(k)+pv(k)*qv(k)

          call find_rot_dumb(d1,d2,c,s)
          sprot(1,k-1)=c
          sprot(2,k-1)=s

          c=sprot(1,k-1)
          s=sprot(2,k-1)
          cc=c
          cs=s

c
c sparse rotate sub and subsub
c
          if(k.gt.2) then
            d1=c*gv(k-2)-s*(-cpqv(k)*pv(k-2))
            gv(k-2)=d1
          endif
c
c sparse rotate sub and diag
c
          d1=c*vd(k-1)-s*gv(k-1)
          d2=cs*vd(k-1)+cc*gv(k-1)
          vd(k-1)=d1
          gv(k-1)=d2
c
c sparse rotate diag and sup
c
          d1=c*voffd(k-1)-s*vd(k)
          d2=cs*voffd(k-1)+cc*vd(k)
          voffd(k-1)=d1
          vd(k)=d2
c
c rotate pv
c
          d1=c*pv(k-1)-s*pv(k)
          d2=cs*pv(k-1)+cc*pv(k)
          pv(k-1)=d1
          pv(k)=d2
c
c apply correction
c
          ddd1=abs(pv(k-1)*qv(k))**2+abs(pv(k)*qv(k))**2
          ddd2=abs(voffd(k-1))**2+abs(vd(k))**2
          if(ddd1.ge.ddd2) then
            pv(k-1)=-voffd(k-1)/qv(k)
          endif
c
c rotate cpqv
c
          d1=c*cpqv(k-1)-s*cpqv(k)
          d2=cs*cpqv(k-1)+cc*cpqv(k)
          cpqv(k-1)=d1
          cpqv(k)=d2

        enddo


c
c sparse rotate column
c
        do k=n,2,-1

          c=sprot(1,k-1)
          s=sprot(2,k-1)
          cc=c
          cs=s
c
c rotate diag and sup column
c
          d1=cc*vd(k-1)-cs*(-pv(k-1)*qv(k))
          d2=s*vd(k-1)+c*(-pv(k-1)*qv(k))
          vd(k-1)=d1
          voffd(k-1)=d2
c
c rotate sub and diag column
c
          d1=cc*gv(k-1)-cs*vd(k)
          d2=s*gv(k-1)+c*vd(k)
          gv(k-1)=d1
          vd(k)=d2
c
c rotate qv
c
          d1=c*qv(k-1)-s*qv(k)
          d2=cs*qv(k-1)+cc*qv(k)
          qv(k-1)=d1
          qv(k)=d2

        enddo

        return
        end


        subroutine find_rot_dumb(x1,x2,c,s)
        implicit real *8 (a-h,o-z)
        complex *16 x1,x2,c,s,v,vv
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   This subroutine finds the 2*2 complex orthogonal transform Q,
c   so that
c           [c -s] [x1]  =  [      0          ]
c           [s  c] [x2]     [sqrt(x1**2+x2**2)]
c
c                       WARNING:
c   Size of c and s are not bounded! They can be large.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        v=sqrt((x1)**2+(x2)**2)
        vv=sqrt(abs(x1)**2+abs(x2)**2)
        if(abs(vv).eq.0) then
         c=1
         s=0
         goto 1000
        endif

        c=x2/v

        s=x1/v

 1000 continue

        return
        end


        subroutine eigenvals_quadr_sol2_elem(a11,a21,a12,a22,
     1    clam1,clam2)
        implicit real *8 (a-h,o-z)
        complex *16 aa(2,2),clam1,clam2,discr,a11,a21,a12,a22
        data four/4.0d0/,half/0.5d0/

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       This subroutine computes the two eigenvalues of a complex
c       2 by 2 matrix whose elements are a11,a21,a12,a22.
c
c       clam1 and clam2 are the two eigenvalues.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        aa(1,1)=a11
        aa(1,2)=a12
        aa(2,1)=a21
        aa(2,2)=a22


        discr=sqrt( (aa(1,1)+aa(2,2))**2+
     1      four*( aa(1,2)*aa(2,1) -aa(1,1)*aa(2,2)) )
c
        clam1=aa(1,1)+aa(2,2)+discr

        clam1=clam1*half
c
        clam2=aa(1,1)+aa(2,2)-discr
        clam2=clam2*half
c
        return
        end


        subroutine cqr_eigenvals_wilk(clam,n,exiteps,ier,
     1              vd,voffd,coefs)
        implicit real *8 (a-h,o-z)
        dimension nstack(10000,3)
        complex *16 clam(*),coefs(*)
        complex *16 vd(*),voffd(*),pv(10 000),qv(10 000),clam1,clam2
        complex *16 a11,a12,a21,a22,offelem
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       This subroutine computes all n eigenvalues of the matrix
c
c                    C = A + pv * qv^T,                                 (1)
c
c       where A is n-by-n complex symmetric tridiagonal (diagonal vd,
c       superdiagonal voffd), and pv, qv are n-vectors.  The vectors
c       pv and qv are constructed internally from voffd and coefs:
c
c            pv = (0, ..., 0, -1),
c            qv(i) = voffd(n) * coefs(i) / coefs(n+1).
c
c       The matrix C is lower Hessenberg, i.e. C(i,j) = 0 for j > i+1.
c       Together with the rank-1 structure, C is fully determined by
c       vd, voffd, and one of pv or qv.
c
c       The algorithm is a shifted QR iteration with Wilkinson shifts,
c       using complex orthogonal (not unitary) Givens rotations to
c       preserve the complex symmetry of A and the rank-1 structure.
c       Each QR sweep costs O(n), and the total cost is O(n^2).
c
c       Numerical experiments suggest componentwise backward stability:
c       the computed eigenvalues are the exact eigenvalues of
c
c                   (A+dA) + (pv+dpv) * (qv+dqv)^T                     (2)
c
c       with dA ~ eps*||A||, dpv ~ eps*||pv||, dqv ~ eps*||qv||.
c       Detailed perturbation analysis has not been done.
c
c                               WARNING:
c       Complex orthogonal rotations can have unbounded norms, so
c       numerical stability is not guaranteed and breakdowns can
c       occur in principle.  However, with Wilkinson shifts the QR
c       converges fast enough that breakdowns have not been observed.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                      input parameters:
c
c   n: size of the matrix
c   vd: diagonal of A, n elements.  Destroyed by the subroutine.
c   voffd: n elements. voffd(1),...,voffd(n-1) are the superdiagonal
c          of A; voffd(n) is used to form qv with coefs.
c          Destroyed by the subroutine.
c   coefs: n+1 elements, used to form qv via
c          qv(i) = voffd(n) * coefs(i) / coefs(n+1).
c   exiteps: ABSOLUTE precision for convergence of off-diagonal
c            elements.
c
c                      output parameters:
c
c   clam: n computed eigenvalues of C
c   ier:
c        0    successful
c        1028 QR fails to converge within 10000 iterations
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        do i=1,n
          qv(i)=voffd(n)*coefs(i)/coefs(n+1)
          pv(i)=0
        enddo

        pv(n)=-1

        iternum=0
        is=1
        maxstack=is
        nstack(is,1)=1
        nstack(is,2)=n
        nstack(is,3)=nstack(is,2)-nstack(is,1)+1


        ier=0
        do 5000 ijk=1,10000
        if(is.eq.0) goto 6000

        nn=nstack(is,3)
        if((nn.eq.1).or.(nn.eq.2)) goto 2000

        do 1000 i=nstack(is,1),(nstack(is,2)-1)
        offelem=voffd(i)+pv(i)*qv(i+1)

        if(abs(offelem).le.exiteps) then

c
c go down the stack
c
        is=is+1
        nstack(is,1)=nstack(is-1,1)
        nstack(is,2)=i
        nstack(is,3)=nstack(is,2)-nstack(is,1)+1

        nstack(is-1,1)=i+1
        nstack(is-1,3)=nstack(is-1,2)-nstack(is-1,1)+1

        if(maxstack.lt.is) maxstack=is

        goto 2000

        endif

 1000 continue
 2000 continue

        ii=nstack(is,1)
        nn=nstack(is,3)

        if(nn.eq.1) then
        a11=vd(ii)+pv(ii)*qv(ii)
        clam(ii)=a11

c
c go up the stack
c
        is=is-1

        goto 3000
        endif

        if(nn.eq.2) then

        a11=vd(ii)+pv(ii)*qv(ii)
        a12=voffd(ii)+pv(ii)*qv(ii+1)
        a21=voffd(ii)+pv(ii+1)*qv(ii)
        a22=vd(ii+1)+pv(ii+1)*qv(ii+1)

        call eigenvals_quadr_sol2_elem(a11,a21,a12,a22,
     1    clam1,clam2)
        clam(ii)=clam1
        clam(ii+1)=clam2

c
c go up the stack
c
        is=is-1

        goto 3000
        endif

        call qr_one_step(nn,vd(ii),voffd(ii),pv(ii),qv(ii))
        iternum=iternum+1

 3000 continue

 5000 continue
        ier=1028
 6000 continue


        if(ier.eq.1028) then
        write(6,*) 'cqr wilk fails to converge after 10000 iterations'
        write(13,*) 'cqr wilk fails to converge after 10000 iterations'
        endif


        return
        end

        subroutine qr_one_step(n,vd,voffd,pv,qv)
        implicit real *8 (a-h,o-z)
        complex *16 vd(*),voffd(*),pv(*),qv(*)
        complex *16 clam1,clam2,shift,a11,a12,a21,a22

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       This subroutine performs one (complex orthogonal) QR iteration
c       to the matrix represented by generators vd, voffd, pv, qv,
c       all containing n complex numbers. The generators are updated
c       by the subroutine.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        a11=vd(1)+pv(1)*qv(1)
        a12=voffd(1)+pv(1)*qv(2)
        a21=voffd(1)+pv(2)*qv(1)
        a22=vd(2)+pv(2)*qv(2)
        call eigenvals_quadr_sol2_elem(a11,a21,a12,a22,
     1    clam1,clam2)

        shift=clam1
        if(abs(clam1-a11) .gt. abs(clam2-a11) ) shift=clam2

c
c apply shift
c
        do i=1,n
          vd(i)=vd(i)-shift
        enddo

        call cqr_one_iter_dumb(n,vd,voffd,pv,qv)

        do i=1,n
          vd(i)=vd(i)+shift
        enddo


        return
        end
