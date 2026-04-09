cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Real rootfinder via Chebyshev expansion and colleague matrix.
c
c     Given a user-supplied real function fun(x,par1,par2,val),
c     finds all real roots on the interval [a,b]. The algorithm
c     expands fun in a Chebyshev basis, constructs the colleague
c     matrix, and finds roots as eigenvalues using the O(n^2)
c     Hermitian + rank-1 QR eigensolver (herm_p_rank1).
c     Optional Newton refinement via Chebyshev derivative.
c
c     Dependencies:
c       herm_p_rank1.f - Hermitian + rank-1 eigenvalue solver
c       No LAPACK dependency.
c
c     Entry points:
c       droots_cheb      - real roots on [a,b] (single interval)
c       droots_cheb_adap - real roots on [a,b] (adaptive subdivision)
c       (for complex roots near [a,b], use zroots_cheb.f)
c
c     Internal subroutines:
c       droots_cheb0         - per-interval core (no Newton, no setup)
c       droots_cheb_dedup    - remove duplicate roots
c       colleague_roots_m1p1 - real roots of Chebyshev expansion on [-1,1]
c       colleague_roots      - all complex roots of Chebyshev expansion
c       droots_chebexps      - Chebyshev nodes, weights, matrices
c       droots_chebexev      - evaluate Chebyshev expansion
c       droots_chfunder      - evaluate expansion and derivative
c       droots_matvec        - matrix-vector multiply
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Find real roots of fun on [a,b] via Chebyshev expansion.
c
c     Maps [a,b] to [-1,1], expands fun in Chebyshev basis,
c     finds roots using colleague matrix, maps back to [a,b].
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine droots_cheb(ifnewton,fun,par1,par2,
     1      a,b,n,eps,roots,nroots,errest,ier)
        implicit real *8 (a-h,o-z)
        real *8 roots(*),errest(*),par1(*),par2(*)
        external fun
c
c       local workspace on the stack
c
        real *8 x(n+5),u(n+5,n+5),v(n+5,n+5),whts(n+5)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       input:
c       ifnewton - 1 to refine roots with Newton
c       fun(x,par1,par2,val,dval) - user function on [a,b],
c            x is real*8, val=f(x), dval=f'(x).
c            dval is only used when ifnewton=1.
c       par1, par2 - user parameters passed to fun
c       a, b - interval endpoints
c       n - expansion order (number of Chebyshev terms is n+1)
c       eps - machine precision for eigenvalue computation
c
c       output:
c       roots - real roots found on [a,b]
c       nroots - number of roots
c       errest - error estimates for each root
c            if ifnewton=0: |f(root)|
c            if ifnewton=1: |f(root)/f'(root)| (Newton step)
c       ier - 0=success, 10=not converged, 128=QR failure
c
c       Internally, the function is expanded in Chebyshev basis
c       on [a,b] and roots are found as eigenvalues of the
c       colleague matrix. Complex eigenvalues near the real axis
c       are filtered by a strip of half-width delta=1d-3
c       relative to (b-a)/2 (i.e. imag part < 1d-3*(b-a)/2).
c       Spurious roots are discarded using the Chebyshev
c       polynomial residual |p(root)| > coff*||coefs|| with
c       coff=1000*eps.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c precompute Chebyshev nodes and matrices
c
        itype=2
        call droots_chebexps(itype,n+1,x,u,v,whts)

c
c find roots on [a,b]
c
        call droots_cheb0(fun,par1,par2,a,b,n,eps,
     1      x,u,roots,nroots,ier)

c
c compute error estimates
c
        do i=1,nroots
          call fun(roots(i),par1,par2,val,dval)
          errest(i)=abs(val)
        enddo

c
c Newton refinement using true function (2 iterations)
c overwrite errest with |f/f'|
c
        if(ifnewton.eq.1) then
          do i=1,nroots
            do ijk=1,2
              call fun(roots(i),par1,par2,val,dval)
              if(abs(dval).gt.0) roots(i)=roots(i)-val/dval
            enddo
            call fun(roots(i),par1,par2,val,dval)
            errest(i)=abs(val/dval)
          enddo
        endif

        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Adaptive real rootfinder via Chebyshev expansion.
c
c     Subdivides [a,b] when the expansion of order n
c     doesn't converge. Uses a stack-based binary
c     subdivision. Newton and errest are applied after
c     dedup on the final root set.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine droots_cheb_adap(ifnewton,fun,par1,par2,
     1      a,b,n,eps,nexdp,roots,nroots,errest,ier)
        implicit real *8 (a-h,o-z)
        real *8 roots(*),errest(*),par1(*),par2(*)
        external fun
c
c       precomputed Chebyshev data (same for all subintervals)
c
        real *8 x(n+5),u(n+5,n+5),v(n+5,n+5),whts(n+5)
c
c       stack for binary subdivision
c
        parameter (maxstack=20 000,maxlevel=25)
        real *8, allocatable :: stack_a(:),stack_b(:)
        integer, allocatable :: stack_ifd(:)
c
c       temporary root storage before dedup
c
        real *8, allocatable :: rootslist(:)
        real *8 rts(n+5)
        parameter (maxroots=10 000)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       input:
c       ifnewton - 1 to refine roots with Newton
c       fun(x,par1,par2,val,dval) - user function
c       par1, par2 - user parameters
c       a, b - interval endpoints
c       n - Chebyshev order per subinterval
c       eps - machine precision
c       nexdp - extra subdivision levels after convergence
c            (recommended 0)
c
c       output:
c       roots - real roots found on [a,b]
c       nroots - number of roots
c       errest - error estimates
c       ier - 0=success, 100=QR failure, 1024=max level,
c            2048=max stack, 4096=max roots
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        allocate(rootslist(maxroots))
        allocate(stack_a(maxstack),stack_b(maxstack))
        allocate(stack_ifd(maxstack))

c
c precompute Chebyshev nodes and matrices
c
        itype=2
        call droots_chebexps(itype,n+1,x,u,v,whts)

c
c initialize BFS interval list
c
        stack_a(1)=a
        stack_b(1)=b
        stack_ifd(1)=-nexdp

        ifirst=1
        ilast=1
        ilastnew=1
        nrtot=0
        nboxes=0
        ier=0

c
c BFS: process all intervals level by level
c
        do level=0,maxlevel

        iffini=1

        do icount=ifirst,ilast
          ai=stack_a(icount)
          bi=stack_b(icount)
          ifdone=stack_ifd(icount)

          call droots_cheb0(fun,par1,par2,ai,bi,n,eps,
     1        x,u,rts,nr,ier_sub)

          nboxes=nboxes+1

          if(ier_sub.eq.0) ifdone=ifdone+1

          if(ier_sub.eq.128) then
c
c QR failure: skip this interval, report
c
            ier=100
            write(6,*) 'QR failed on [',ai,',',bi,']'
            write(13,*) 'QR failed on [',ai,',',bi,']'

          elseif(ifdone.ge.1) then
c
c converged: collect roots
c
            if(nrtot+nr.gt.maxroots) then
              ier=4096
              goto 2000
            endif
            do i=1,nr
              rootslist(i+nrtot)=rts(i)
            enddo
            nrtot=nrtot+nr

          else
c
c not done: subdivide into two halves
c
            iffini=0
            amid=(ai+bi)/2

            ilastnew=ilastnew+1
            if(ilastnew.gt.maxstack) then
              ier=2048
              goto 2000
            endif
            stack_a(ilastnew)=ai
            stack_b(ilastnew)=amid
            stack_ifd(ilastnew)=ifdone

            ilastnew=ilastnew+1
            if(ilastnew.gt.maxstack) then
              ier=2048
              goto 2000
            endif
            stack_a(ilastnew)=amid
            stack_b(ilastnew)=bi
            stack_ifd(ilastnew)=ifdone
          endif

        enddo

        if(iffini.eq.1) goto 2000

        ifirst=ilast+1
        ilast=ilastnew

        enddo

        ier=1024

 2000   continue

        write(6,*) 'total intervals is',nboxes
        write(13,*) 'total intervals is',nboxes
        write(6,*) 'max subdivisions is',level
        write(13,*) 'max subdivisions is',level

c
c Newton on rootslist before dedup (2 iterations)
c pulls duplicates from adjacent subintervals together
c
        if(ifnewton.eq.1) then
          do i=1,nrtot
            do ijk=1,2
              call fun(rootslist(i),par1,par2,val,dval)
              if(abs(dval).gt.0)
     1            rootslist(i)=rootslist(i)-val/dval
            enddo
          enddo
        endif

c
c dedup roots
c
        if(nboxes.gt.1) then
          call droots_cheb_dedup(eps,rootslist,nrtot,roots,nroots)
        else
          nroots=nrtot
          do i=1,nroots
            roots(i)=rootslist(i)
          enddo
        endif

        deallocate(rootslist,stack_a,stack_b,stack_ifd)

c
c compute error estimates
c
        do i=1,nroots
          call fun(roots(i),par1,par2,val,dval)
          errest(i)=abs(val)
          if(ifnewton.eq.1) errest(i)=abs(val/dval)
        enddo

        write(6,*) 'total roots found',nroots
        write(13,*) 'total roots found',nroots

        return
        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Remove duplicate roots from adjacent intervals.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine droots_cheb_dedup(eps,roots1,n1,roots2,n2)
        implicit real *8 (a-h,o-z)
        real *8 roots1(*),roots2(*)

        roots2(1)=roots1(1)
        n2=1

        tol=sqrt(eps)*1d-2

        do i=2,n1

        ifput=1
        do j=1,n2
          d=abs(roots2(j)-roots1(i))
          if(d.lt.tol*max(1d0,abs(roots1(i)))) ifput=0
        enddo

        if(ifput.eq.1) then
          n2=n2+1
          roots2(n2)=roots1(i)
        endif

        enddo

        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Core per-interval rootfinder. Precomputed x, u are
c     passed in. No Newton, no errest — just finds roots.
c     Returns ier=10 if expansion doesn't converge.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine droots_cheb0(fun,par1,par2,a,b,n,eps,
     1      x,u,roots,nroots,ier)
        implicit real *8 (a-h,o-z)
        real *8 roots(*),par1(*),par2(*)
        real *8 x(*),u(n+1,n+1)
        external fun
c
c       local workspace on the stack
c
        real *8 y(n+5),coefs(n+5)

        delta=1d-3
        coff=sqrt(eps)*1d-2

c
c evaluate fun at mapped nodes
c
        uu=(b+a)/2
        vv=(b-a)/2

        do i=1,n+1
          call fun(uu+vv*x(i),par1,par2,val,dval)
          y(i)=val
        enddo

c
c check for Inf/NaN in function values (e.g. pole on a
c Chebyshev node). y.ne.y detects NaN (IEEE 754: NaN != NaN).
c abs(y)>1d300 detects Inf and near-overflow (double
c precision overflows at ~1.8d308). return ier=10 so
c the adaptive rootfinder subdivides around the pole.
c
        do i=1,n+1
          if(y(i).ne.y(i) .or. abs(y(i)).gt.1d300) then
            ier=10
            nroots=0
            return
          endif
        enddo

c
c compute Chebyshev coefficients
c
        call droots_matvec(u,y,n+1,coefs)

c
c check convergence
c
        coefsmax=0
        do i=1,n+1
          if(abs(coefs(i)).gt.coefsmax) coefsmax=abs(coefs(i))
        enddo

        coefscut=0
        istart=n*(4d0/5d0)
        do i=istart,n+1
          if(abs(coefs(i)).gt.coefscut) coefscut=abs(coefs(i))
        enddo

        nroots=0
        if(coefscut.gt.eps**(3d0/4d0)*coefsmax) then
          ier=10
          return
        endif

c
c find real roots on [-1,1], map back to [a,b]
c
        call colleague_roots_m1p1(ier,coefs,n,eps,delta,coff,
     1      roots,nroots)

        do i=1,nroots
          roots(i)=uu+vv*roots(i)
        enddo

        return
        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Colleague matrix rootfinding subroutines.
c
c     colleague_roots_m1p1 - real roots of Chebyshev expansion on [-1,1]
c     colleague_roots      - all complex roots of Chebyshev expansion
c
c     Both use the O(n^2) Hermitian + rank-1 QR eigensolver
c     (herm_p_rank1) and are backward stable.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine colleague_roots_m1p1(ier,coefs,n,eps,delta,coff,
     1      roots,nroots)
        implicit real *8 (a-h,o-z)
        real *8 coefs(*),roots(*)
        complex *16 croots(n+5),ima
        real *8 roots2(n+5)
        data ima/(0.0d0,1.0d0)/
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Finds the real roots of the Chebyshev expansion
c
c           f(x) = sum_{i=0}^{n} coefs(i+1) * T_i(x)
c
c       on the interval [-1,1]. Finds all n complex roots
c       via colleague_roots, then filters for roots in the
c       strip [-1-delta,1+delta] x [-delta,delta] and
c       discards those where |f(root)| > coff*||coefs||.
c
c       input:
c       coefs - Chebyshev expansion coefficients (n+1 entries)
c       n - order of the expansion (n+1 terms)
c       eps - absolute precision for eigenvalue computation
c       delta - half-width of complex strip for filtering
c       coff - cutoff for discarding spurious roots
c
c       output:
c       roots - real roots found on [-1,1]
c       nroots - number of roots
c       ier - 0=success, 128=QR failure
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        call colleague_roots(ier,coefs,n,eps,croots)
        if (ier .ne. 0) return
c
c get roots near [-1,1]
c
        ii=0
        do i=1,n
          xx=croots(i)
          yy=-ima*croots(i)

          if (xx .lt. -1-delta .or. xx .gt. 1+delta) goto 1400
          if (yy .lt. -delta .or. yy .gt. delta) goto 1400

          ii=ii+1
          roots2(ii)=croots(i)
 1400   continue
        enddo

        nroots2=ii
c
c compute norm of coefs
c
        dd=0
        do i=1,n+1
          dd=dd+coefs(i)**2
        enddo
        dd=sqrt(dd)
c
c discard roots with large error
c
        ii=0
        do i=1,nroots2
          call droots_chebexev(roots2(i),val,coefs,n)
          if (abs(val) .gt. coff*dd) goto 1800
          ii=ii+1
          roots(ii)=roots2(i)
 1800   continue
        enddo

        nroots=ii

        return
        end


        subroutine colleague_roots(ier,coefs,n,eps,roots)
        implicit real *8 (a-h,o-z)
        complex *16 roots(*)
        real *8 coefs(*)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Finds all n complex roots of the Chebyshev expansion
c
c           f(z) = sum_{i=0}^{n} coefs(i+1) * T_i(z)
c
c       by computing the eigenvalues of the colleague matrix
c       using the O(n^2) Hermitian + rank-1 QR eigensolver.
c       Backward stable: computed roots are exact roots of
c       a perturbed expansion with ||dcoefs|| < eps*||coefs||.
c
c       input:
c       coefs - Chebyshev expansion coefficients (n+1 entries)
c       n - order of the expansion (n+1 terms)
c       eps - absolute precision for eigenvalue computation
c
c       output:
c       roots - n complex roots (complex*16)
c       ier - 0=success, 128=QR failure
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       local workspace on the stack
c
        complex *16 supdiag(n+5),p(n+5),q(n+5),w(15*n+100)

        done=1
c
c construct the colleague matrix
c
        do i=2,n-1
          supdiag(i)=done/2
        enddo
        supdiag(1)=1/sqrt(2*done)

        do i=1,n
          p(i)=0
          roots(i)=0
        enddo
        p(n)=1

        do i=1,n
          q(i)=-coefs(i)/2/coefs(n+1)
        enddo
        q(1)=q(1)*sqrt(2*done)
c
c find the eigenvalues
c
        call herm_p_rank1(ier,roots,supdiag,p,q,n,eps,
     1      w,niter_max,aver_niter)

        return
        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Chebyshev utility subroutines (from chebexps.f, renamed
c     with droots_ prefix to avoid name collisions).
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine droots_chebexps(itype,n,x,u,v,whts)
        implicit real *8 (a-h,o-z)
        real *8 x(*),whts(*),u(n,n),v(n,n)
c
c       Constructs Chebyshev nodes on [-1,1], quadrature weights,
c       and the matrices u (values->coefficients) and
c       v (coefficients->values).
c
c       input:
c       itype - 0: nodes only, 1: nodes+weights, 2: nodes+weights+u+v
c       n - number of nodes
c
c       output:
c       x - Chebyshev nodes
c       u - values to coefficients matrix (n x n)
c       v - coefficients to values matrix (n x n)
c       whts - quadrature weights
c
        done=1
        pi=atan(done)*4
        h=pi/(2*n)
        do i=1,n
          t=(2*i-1)*h
          x(n-i+1)=cos(t)
        enddo

        if(itype .eq. 0) return
c
c construct weights and first two rows of u
c
        if(itype .ge. 2) then
          do i=1,n
            u(1,i)=1
            u(2,i)=x(i)
          enddo
        endif

        do i=1,n
          tjm2=1
          tjm1=x(i)
          whts(i)=2

          ic=-1
          do j=2,n-1
            tj=2*x(i)*tjm1-tjm2
            if(itype .eq. 2) u(j+1,i)=tj
            tjm2=tjm1
            tjm1=tj

            ic=-ic
            if(ic .lt. 0) goto 1400
            rint=-2*(done/(j+1)-done/(j-1))
            whts(i)=whts(i)-rint*tj
 1400       continue
          enddo
          whts(i)=whts(i)/n
        enddo

        if(itype .ne. 2) return
c
c normalize u
c
        do i=1,n
          d=0
          do j=1,n
            d=d+u(i,j)**2
          enddo
          d=done/dsqrt(d)
          do j=1,n
            u(i,j)=u(i,j)*d
          enddo
        enddo
c
c rescale u
c
        ddd=dsqrt(2*done)
        dd=done/dsqrt(n/2d0)
        do i=1,n
          do j=1,n
            u(j,i)=u(j,i)*dd
          enddo
          u(1,i)=u(1,i)/ddd
        enddo
c
c construct v and swap u/v
c
        do i=1,n
          do j=1,n
            v(j,i)=u(i,j)*n/2d0
          enddo
        enddo

        do i=1,n
          do j=1,n
            d=v(j,i)
            v(j,i)=u(j,i)*n/2d0
            u(j,i)=d/n*2d0
          enddo
        enddo

        do i=1,n
          v(1,i)=v(1,i)*2
        enddo

        do i=1,n
          do j=1,i-1
            d=u(j,i)
            u(j,i)=u(i,j)
            u(i,j)=d

            d=v(j,i)
            v(j,i)=v(i,j)
            v(i,j)=d
          enddo
        enddo

        return
        end


        subroutine droots_chebexev(x,val,texp,n)
        implicit real *8 (a-h,o-z)
        real *8 texp(*)
c
c       Evaluates a Chebyshev expansion at point x in [-1,1].
c
c       input: x, texp (n+1 coefficients), n (order)
c       output: val
c
        pjm2=1
        pjm1=x

        val=texp(1)*pjm2+texp(2)*pjm1

        do j=2,n
          pj=2*x*pjm1-pjm2
          val=val+texp(j+1)*pj
          pjm2=pjm1
          pjm1=pj
        enddo

        return
        end


        subroutine droots_chfunder(x,val,der,texp,n)
        implicit real *8 (a-h,o-z)
        real *8 texp(*)
c
c       Evaluates a Chebyshev expansion and its derivative
c       at point x in [-1,1].
c
c       input: x, texp (n+1 coefficients), n (order)
c       output: val, der
c
        tjm2=1
        tjm1=x
        derjm2=0
        derjm1=1

        val=texp(1)*tjm2+texp(2)*tjm1
        der=texp(2)

        do j=2,n
          tj=2*x*tjm1-tjm2
          val=val+texp(j+1)*tj

          derj=2*tjm1+2*x*derjm1-derjm2
          der=der+texp(j+1)*derj

          tjm2=tjm1
          tjm1=tj
          derjm2=derjm1
          derjm1=derj
        enddo

        return
        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Matrix-vector multiply y = A*x (n x n)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine droots_matvec(a,x,n,y)
        implicit real *8 (a-h,o-z)
        real *8 a(n,n),x(n),y(n)

        do i=1,n
          y(i)=0
          do j=1,n
            y(i)=y(i)+a(i,j)*x(j)
          enddo
        enddo

        return
        end
