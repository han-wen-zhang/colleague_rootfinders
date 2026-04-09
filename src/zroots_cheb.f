cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Complex rootfinder via Chebyshev expansion on a real interval.
c
c     Given a complex analytic function fun(z,par1,par2,val,dval),
c     finds complex roots near [a,b] in a strip of half-width delta.
c     The function is expanded in Chebyshev basis on real nodes.
c     The derivative dval=f'(z) is used for Newton refinement.
c
c     Dependencies: herm_p_rank1.f only, no LAPACK.
c
c     Entry point:
c       zroots_cheb        - complex roots near [a,b] (single interval)
c       zroots_cheb_adap   - complex roots near [a,b] (adaptive)
c
c     Internal subroutines:
c       zroots_cheb0       - per-interval core (no Newton, no setup)
c       zroots_cheb_dedup  - remove duplicate complex roots
c       zcolleague_roots_m1p1 - roots near [-1,1] with Bernstein+coff filter
c       zcolleague_roots   - all complex roots of Chebyshev expansion
c       zroots_chebexev    - evaluate complex Chebyshev expansion
c       zroots_chfunder    - expansion and derivative at complex point
c       zroots_chebexps    - Chebyshev nodes, weights, matrices
c       zroots_matvec      - real matrix * complex vector
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine zroots_cheb(ifnewton,fun,par1,par2,
     1      a,b,n,eps,delta,croots,nroots,errest,ier)
        implicit real *8 (a-h,o-z)
        complex *16 croots(*)
        real *8 errest(*),par1(*),par2(*)
        complex *16 zval,zdval
        external fun
c
c       local workspace on the stack
c
        real *8 x(n+5),u(n+5,n+5),v(n+5,n+5),whts(n+5)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       input:
c       ifnewton - 1 to refine roots with Newton
c       fun(z,par1,par2,val,dval) - complex analytic function,
c            z is complex*16, val=f(z), dval=f'(z).
c            dval is only used when ifnewton=1. if the
c            derivative is not available, set ifnewton=0
c            and dval will not be accessed.
c       par1, par2 - user parameters passed to fun
c       a, b - interval endpoints
c       n - expansion order (n+1 Chebyshev terms)
c       eps - absolute precision for eigenvalue computation
c       delta - controls how far from [a,b] complex roots
c            are sought. Roots are kept inside a Bernstein
c            ellipse with foci at the endpoints of [a,b].
c            On [-1,1] (mapped space), the ellipse has:
c              semi-minor axis = delta
c              semi-major axis - 1 ~ delta^2/2 (for small delta)
c            The ellipse narrows near the endpoints,
c            matching the Chebyshev expansion's natural
c            region of accuracy. e.g. delta=0.1 finds
c            roots within ~0.1 of the real axis at the
c            interval center. Keep delta small — accuracy
c            degrades for roots far from [a,b].
c
c       output:
c       croots - complex roots near [a,b] (complex*16)
c       nroots - number of roots
c       errest - error estimates for each root
c            if ifnewton=0: |fun(root)| (residual)
c            if ifnewton=1: |fun(root)/fun'(root)|
c                 (Newton step size)
c       ier - 0=success, 10=expansion not converged,
c            128=QR failure
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c precompute Chebyshev nodes and matrices
c
        itype=2
        call zroots_chebexps(itype,n+1,x,u,v,whts)

c
c convert physical delta to Bernstein ellipse parameter rho.
c the Bernstein ellipse has foci at +-1 in mapped space:
c   semi-minor axis = (rho - 1/rho) / 2
c   semi-major axis = (rho + 1/rho) / 2
c setting semi-minor = dd = delta/((b-a)/2) and solving:
c   rho = dd + sqrt(dd^2 + 1)
c in physical space, the semi-minor axis is delta.
c
        vv=(b-a)/2
        dd=delta/vv
        rho=dd+sqrt(dd**2+1)

c
c find roots on [a,b]
c
        call zroots_cheb0(fun,par1,par2,a,b,n,eps,rho,
     1      x,u,croots,nroots,ier)
        if(ier.ne.0) return

c
c Newton refinement using true function (2 iterations)
c
        if(ifnewton.eq.1) then
          do i=1,nroots
            do ijk=1,2
              call fun(croots(i),par1,par2,zval,zdval)
              if(abs(zdval).gt.0)
     1            croots(i)=croots(i)-zval/zdval
            enddo
          enddo
        endif

c
c compute error estimates
c
        do i=1,nroots
          call fun(croots(i),par1,par2,zval,zdval)
          errest(i)=abs(zval)
          if(ifnewton.eq.1) errest(i)=abs(zval/zdval)
        enddo

        return
        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Core per-interval complex rootfinder. Precomputed
c     x, u are passed in. No Newton, no errest.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine zroots_cheb0(fun,par1,par2,a,b,n,eps,delta,
     1      x,u,croots,nroots,ier)
        implicit real *8 (a-h,o-z)
        complex *16 croots(*)
        real *8 par1(*),par2(*)
        real *8 x(*),u(n+1,n+1)
        external fun
c
c       local workspace on the stack
c
        complex *16 y(n+5),coefs(n+5),zz,zdval

c
c evaluate fun at mapped nodes
c
        uu=(b+a)/2
        vv=(b-a)/2

        do i=1,n+1
          zz=uu+vv*x(i)
          call fun(zz,par1,par2,y(i),zdval)
        enddo

c
c check for Inf/NaN in function values
c
        do i=1,n+1
          if(abs(y(i)).ne.abs(y(i)) .or.
     1       abs(y(i)).gt.1d300) then
            ier=10
            nroots=0
            return
          endif
        enddo

c
c compute complex Chebyshev coefficients
c
        call zroots_matvec(u,y,n+1,coefs)

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
c find roots near [-1,1]: Bernstein ellipse + coff filter
c delta here is rho (Bernstein parameter), passed from caller
c
        coff=sqrt(eps)
        call zcolleague_roots_m1p1(ier,coefs,n,eps,delta,coff,
     1      croots,nroots)

c
c map roots back to [a,b]
c
        do i=1,nroots
          croots(i)=croots(i)*vv+uu
        enddo

        return
        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Complex roots of a Chebyshev expansion near [-1,1].
c     Finds all roots via zcolleague_roots, then filters
c     by delta strip and coff residual (via chebexev).
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine zcolleague_roots_m1p1(ier,coefs,n,eps,rho,coff,
     1      roots,nroots)
        implicit real *8 (a-h,o-z)
        complex *16 coefs(*),roots(*),ima
        complex *16 croots(n+5),zval,zw
        data ima/(0d0,1d0)/
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Finds complex roots of the Chebyshev expansion
c       near [-1,1]. Calls zcolleague_roots for all roots,
c       then filters by Bernstein ellipse and coff residual.
c
c       input:
c       coefs - complex Chebyshev coefficients (n+1 entries)
c       n - order of the expansion
c       eps - absolute precision for eigenvalue computation
c       rho - Bernstein ellipse parameter: keep roots with
c            |w| < rho where w = t + sqrt(t^2-1) is the
c            Joukowski preimage. rho=1 is [-1,1] itself.
c            rho=1.01 keeps roots very near [-1,1].
c       coff - cutoff: roots with |p(root)| > coff*||coefs||
c            are discarded.
c
c       output:
c       roots - complex roots near [-1,1]
c       nroots - number of roots
c       ier - 0=success, 128=QR failure
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        call zcolleague_roots(ier,coefs,n,eps,croots)
        if(ier.ne.0) return
c
c Bernstein ellipse filter: keep roots with |w| < rho
c w = t + sqrt(t^2-1) is the Joukowski preimage
c
        ii=0
        do i=1,n
          call zroots_joukowski(croots(i),zw)
          if(abs(zw).lt.rho) then
            ii=ii+1
            roots(ii)=croots(i)
          endif
        enddo
        nroots2=ii
c
c coff filter: discard roots with large |p(root)|
c
        dd=0
        do i=1,n+1
          dd=dd+abs(coefs(i))**2
        enddo
        dd=sqrt(dd)

        nroots=0
        do i=1,nroots2
          call zroots_chebexev(roots(i),zval,coefs,n)
          if(abs(zval).le.coff*dd) then
            nroots=nroots+1
            roots(nroots)=roots(i)
          endif
        enddo

        return
        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     All complex roots of a Chebyshev expansion with
c     complex coefficients, via the colleague matrix.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine zcolleague_roots(ier,coefs,n,eps,roots)
        implicit real *8 (a-h,o-z)
        complex *16 roots(*),coefs(*)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Finds all n complex roots of the Chebyshev expansion
c
c           f(z) = sum_{i=0}^{n} coefs(i+1) * T_i(z)
c
c       with complex coefficients, by computing the eigenvalues
c       of the colleague matrix using the O(n^2) Hermitian +
c       rank-1 QR eigensolver.
c
c       input:
c       coefs - complex Chebyshev coefficients (n+1 entries)
c       n - order of the expansion
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
          q(i)=conjg(-coefs(i)/2/coefs(n+1))
        enddo
        q(1)=q(1)*sqrt(2*done)
c
c find the eigenvalues
c
        call herm_p_rank1(ier,roots,supdiag,p,q,n,eps,
     1      w,niter_max,aver_niter)

        return
        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Chebyshev utilities for complex expansions.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine zroots_chebexev(z,val,texp,n)
        implicit real *8 (a-h,o-z)
        complex *16 z,val,texp(*),pjm2,pjm1,pj
c
c       Evaluates complex Chebyshev expansion at complex point z.
c
        pjm2=1
        pjm1=z

        val=texp(1)*pjm2+texp(2)*pjm1

        do j=2,n
          pj=2*z*pjm1-pjm2
          val=val+texp(j+1)*pj
          pjm2=pjm1
          pjm1=pj
        enddo

        return
        end


        subroutine zroots_chfunder(z,val,der,texp,n)
        implicit real *8 (a-h,o-z)
        complex *16 z,val,der,texp(*)
        complex *16 tjm2,tjm1,tj,derjm2,derjm1,derj
c
c       Evaluates complex Chebyshev expansion and its derivative
c       at complex point z.
c
        tjm2=1
        tjm1=z
        derjm2=0
        derjm1=1

        val=texp(1)*tjm2+texp(2)*tjm1
        der=texp(2)

        do j=2,n
          tj=2*z*tjm1-tjm2
          val=val+texp(j+1)*tj

          derj=2*tjm1+2*z*derjm1-derjm2
          der=der+texp(j+1)*derj

          tjm2=tjm1
          tjm1=tj
          derjm2=derjm1
          derjm1=derj
        enddo

        return
        end


        subroutine zroots_chebexps(itype,n,x,u,v,whts)
        implicit real *8 (a-h,o-z)
        real *8 x(*),whts(*),u(n,n),v(n,n)
c
c       Constructs Chebyshev nodes, weights, and matrices on [-1,1].
c       Same as droots_chebexps (real nodes and matrices).
c
        done=1
        pi=atan(done)*4
        h=pi/(2*n)
        do i=1,n
          t=(2*i-1)*h
          x(n-i+1)=cos(t)
        enddo

        if(itype .eq. 0) return

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

        ddd=dsqrt(2*done)
        dd=done/dsqrt(n/2d0)
        do i=1,n
          do j=1,n
            u(j,i)=u(j,i)*dd
          enddo
          u(1,i)=u(1,i)/ddd
        enddo

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


        subroutine zroots_matvec(a,x,n,y)
        implicit real *8 (a-h,o-z)
        real *8 a(n,n)
        complex *16 x(n),y(n)
c
c       Real matrix times complex vector: y = A * x
c
        do i=1,n
          y(i)=0
          do j=1,n
            y(i)=y(i)+a(i,j)*x(j)
          enddo
        enddo

        return
        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Joukowski preimage: w = t + sqrt(t^2 - 1)
c     |w| = 1 on [-1,1], |w| > 1 outside.
c     Takes the branch with |w| >= 1.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine zroots_joukowski(t,w)
        implicit real *8 (a-h,o-z)
        complex *16 t,w,w2

        w=t+sqrt(t**2-1)
        w2=t-sqrt(t**2-1)
        if(abs(w2).gt.abs(w)) w=w2

        return
        end
        subroutine zroots_cheb_adap(ifnewton,fun,par1,par2,
     1      a,b,n,eps,delta,nexdp,
     2      croots,nroots,errest,ier)
        implicit real *8 (a-h,o-z)
        complex *16 croots(*)
        real *8 errest(*),par1(*),par2(*)
        external fun
c
c       precomputed Chebyshev data
c
        real *8 x(n+5),u(n+5,n+5),v(n+5,n+5),whts(n+5)
c
c       BFS interval list
c
        parameter (maxstack=20 000,maxlevel=25)
        real *8, allocatable :: stack_a(:),stack_b(:)
        integer, allocatable :: stack_ifd(:)
c
c       temporary root storage
c
        complex *16 rts(n+5)
        complex *16, allocatable :: rootslist(:)
        complex *16 zval,zdval
        parameter (maxroots=10 000)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       input:
c       ifnewton - 1 to refine roots with Newton
c       fun(z,par1,par2,val,dval) - complex analytic function
c       par1, par2 - user parameters
c       a, b - real interval endpoints
c       n - Chebyshev order per subinterval
c       eps - machine precision
c       delta - controls Bernstein ellipse for root filtering.
c            On the top-level [a,b], delta is the semi-minor
c            axis (same as zroots_cheb). Internally, delta is
c            converted to a Bernstein parameter rho which is
c            fixed for all subintervals. As subintervals shrink,
c            the physical ellipse shrinks proportionally.
c            NOTE: there is no guarantee that all roots inside
c            the top-level ellipse are found. Roots far from
c            the real axis may be outside the subinterval
c            ellipses. Increasing delta improves coverage but
c            may admit more spurious roots.
c       nexdp - extra subdivision levels (recommended 0)
c
c       output:
c       croots - complex roots near [a,b]
c       nroots - number of roots
c       errest - error estimates
c       ier - 0=success, 100=QR failure, 1024=max level,
c            2048=max intervals, 4096=max roots
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        allocate(rootslist(maxroots))
        allocate(stack_a(maxstack),stack_b(maxstack))
        allocate(stack_ifd(maxstack))

c
c precompute Chebyshev nodes and matrices
c
        itype=2
        call zroots_chebexps(itype,n+1,x,u,v,whts)

c
c convert physical delta to Bernstein rho (once for all subintervals)
c
        vv0=(b-a)/2
        dd=delta/vv0
        rho=dd+sqrt(dd**2+1)

c
c initialize BFS
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

          call zroots_cheb0(fun,par1,par2,ai,bi,n,eps,rho,
     1        x,u,rts,nr,ier_sub)

          nboxes=nboxes+1

          if(ier_sub.eq.0) ifdone=ifdone+1

          if(ier_sub.eq.128) then
c
c QR failure: skip, report
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
c not done: subdivide
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
c Post-processing pipeline for the complex adaptive case:
c
c 1. Filter by |fun(root)| < sqrt(eps).
c 2. Newton refinement.
c 3. Dedup.
c 4. Error estimates.
c

c
c step 1: Newton refinement
c 2. Dedup with tol=sqrt(eps)*1e-2.
c 3. Error estimates.
c

c
c step 1: Newton refinement
c
        if(ifnewton.eq.1) then
          do i=1,nrtot
            do ijk=1,2
              call fun(rootslist(i),par1,par2,zval,zdval)
              if(abs(zdval).gt.0)
     1            rootslist(i)=rootslist(i)-zval/zdval
            enddo
          enddo
        endif

c
c step 3: dedup roots (tight tol after Newton)
c
        if(nboxes.gt.1) then
          call zroots_cheb_dedup(eps,rootslist,nrtot,
     1        croots,nroots)
        else
          nroots=nrtot
          do i=1,nroots
            croots(i)=rootslist(i)
          enddo
        endif

        deallocate(rootslist,stack_a,stack_b,stack_ifd)

c
c step 4: compute error estimates
c
        do i=1,nroots
          call fun(croots(i),par1,par2,zval,zdval)
          errest(i)=abs(zval)
          if(ifnewton.eq.1) errest(i)=abs(zval/zdval)
        enddo

        write(6,*) 'total roots found',nroots
        write(13,*) 'total roots found',nroots

        return
        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Remove duplicate complex roots.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine zroots_cheb_dedup(eps,roots1,n1,roots2,n2)
        implicit real *8 (a-h,o-z)
        complex *16 roots1(*),roots2(*)

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
