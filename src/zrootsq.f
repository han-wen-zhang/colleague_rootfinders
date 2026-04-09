cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Complex rootfinder on a square region of the complex plane.
c
c     Given a user-supplied analytic function fun(val,dval,z) and
c     its derivative, finds all roots inside a square of side sqw
c     centered at center.
c
c     The algorithm expands fun on the boundary using polynomials
c     orthogonal on the boundary (evaluated via three-term
c     recurrence), checks convergence, and finds roots as
c     eigenvalues of the generalized colleague matrix via an
c     O(n^2) QR algorithm with complex orthogonal rotations.
c     Optional Newton refinement and residue-based root counting.
c
c     Dependencies:
c       csym_p_rank1.f  - colleague matrix eigenvalue solver
c       LAPACK zgeqrf,zunmqr,ztrtrs - precomputed QR least squares
c
c     Entry points:
c       zrootsq           - single-box driver
c       zrootsq_adap      - adaptive driver (quadtree subdivision)
c
c     Internal subroutines:
c       zrootsq0           - single-box expansion + root extraction
c       zrootsq_divide     - subdivide square into 4 children
c       zrootsq_dedup      - remove duplicate roots
c       zrootsq_pts_wts    - Gauss-Legendre nodes on the square
c       zrootsq_poly_load  - load recurrence and evaluate basis
c       zrootsq_factor     - precompute QR of weighted basis
c       zrootsq_leastsq1   - weighted least squares via precomputed QR
c       zrootsq_residue    - root counting via contour integral
c       zrootsq_residue_int - integrand f'/f for residue formula
c       zroot_adapgau      - adaptive Gaussian quadrature
c       zroot_adinrec      - recursive subdivision for adapgau
c       zroot_oneint       - single-interval quadrature
c       zroot_legewhts     - Gauss-Legendre nodes and weights
c       zroot_legepol_sum  - Legendre polynomial evaluation
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine zrootsq(ifprint,ifnewton,ifres,
     1      fun,par1,par2,norder,eps,center,sqw,
     2      nrtot,errest,roots,ier)
        implicit real *8 (a-h,o-z)
        complex *16 center,roots(*)
        dimension errest(*),par1(*),par2(*)
c
c       max norder = 150, npt = 4*150 = 600
c
        complex *16 z(600),cw(600)
        complex *16 q(600,151),aqr(600,151),tau(151)
        complex *16 vd(150),voffd(150)
        complex *16 coefs(151),rootsall(150),rts(150)
        complex *16 ima,z0,val,dval
        data ima/(0.0d0,1.0d0)/
c
c       Finds all complex roots of fun in a square of side sqw
c       centered at center.
c
c       The algorithm expands fun on the boundary using orthogonal
c       polynomials of order norder (evaluated via three-term
c       recurrence). When the expansion converges to accuracy
c       eps, roots are found as eigenvalues of the colleague
c       matrix. Optional Newton refinement and residue-based root
c       counting are available.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                   input parameters:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       ifprint: 1 to print expansion coefficients for diagnostics
c       ifnewton: 1 to refine roots with Newton iterations
c       ifres: 1 to use residue formula for root counting
c       fun: the function whose roots are to be found, of the form
c            fun(z,par1,par2,val,dval), where val=f(z), dval=f'(z).
c            The derivative dval MUST be provided.
c       par1, par2: user parameters passed through to fun.
c            Can be arrays or scalars, used as needed by fun.
c       norder: order of the expansion (max 150,
c               recommended 40 for double precision)
c       eps: relative accuracy of polynomial expansion.
c               Sets the lower bound on root precision.
c               Recommended: 1d-4 to 1d-12 (double precision)
c       center: center of the square (complex*16)
c       sqw: side length of the square
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                   output parameters:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       nrtot: total number of roots found
c       errest: error estimates for each root (length nrtot)
c               if ifnewton=0: |f(root_i)|  (residual)
c               if ifnewton=1: |f(root_i)/f'(root_i)|  (Newton step)
c       roots: all roots found (length nrtot)
c       ier: error flag
c            0   success
c            10  expansion did not converge (norder too small)
c            100 QR eigenvalue solver failed
c

c
c points on the square and weights
c
        n = norder
        call zrootsq_pts_wts(z,cw,n)

c
c load recurrence coefficients and evaluate polynomial basis
c
        call zrootsq_poly_load(q,vd,voffd,norder,z,4*n)

c
c precompute QR factorization of weighted basis matrix
c
        call zrootsq_factor(q,cw,n,norder,aqr,tau)

        call zrootsq0(eps,ifdone,fun,par1,par2,center,sqw,
     1      vd,voffd,q,aqr,tau,z,cw,n,rootsall,norder,rts,nr,
     2      ierqr,ifprint)

        ier=0
        nrtot=0

c       expansion did not converge
        if(ifdone.eq.0) then
          ier=10
          goto 4100
        endif

c       QR eigenvalue solver failed
        if(ierqr.ne.0) then
          ier=100
          goto 4100
        endif
        do i=1,nr
          call fun(rts(i),par1,par2,val,dval)
          roots(i)=rts(i)
          nrtot=nrtot+1
          errest(i)=abs(val)
        enddo

        write(6,*) 'total root number found by the root finder',nrtot
        write(13,*) 'total root number found by the root finder',nrtot

c
c residue for total root num
c
        if(ifres.eq.1) then
         call zrootsq_residue(fun,par1,par2,nres,maxrec,center,sqw)
         write(6,*) 'total number of roots from residue is',nres
         write(13,*) 'total number of roots from residue is',nres

         ntor=5+log10(sqw)/log10(2d0)

         if(maxrec.gt.ntor) then
           write(6,*) 'Roots may be near the boundary, maxrec=',maxrec
           write(13,*) 'Roots may be near the boundary, maxrec=',maxrec
         endif

        endif

c
c refine with newton
c
        if(ifnewton.eq.1) then
         dzmax=-1
         do i=1,nrtot
           z0=roots(i)
           do ijk=1,2
             call fun(z0,par1,par2,val,dval)
             dzmag=abs(val/dval)
             if(dzmag.gt.dzmax) dzmax=dzmag
             z0=z0-val/dval
           enddo
           roots(i)=z0
         enddo
         write(6,*) 'max newton step size is',dzmax
         write(13,*) 'max newton step size is',dzmax

c
c update error after newton
c
         do i=1,nrtot
           call fun(roots(i),par1,par2,val,dval)
           errest(i)=abs(val/dval)
         enddo
        endif

 4100 continue

        return
        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Adaptive complex rootfinder on a square region.
c
c     Finds all roots of fun inside a square by adaptive
c     quadtree subdivision. Each box is expanded using
c     zrootsq0; boxes that don't converge are subdivided
c     into 4 children. Uses a depth-first stack traversal.
c
c     No dependency on prini/prinf/prin2.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine zrootsq_adap(ifprint,ifnewton,ifres,
     1      fun,par1,par2,norder,eps,nexdp,center,sqw,
     2      nrtot,errest,roots,centers,nc,ier)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*),errest(*)
        complex *16 center,roots(*),centers(*)
c
c       max norder = 150, npt = 4*150 = 600
c
        complex *16 z(600),cw(600)
        complex *16 q(600,151),aqr(600,151),tau(151)
        complex *16 vd(150),voffd(150)
        complex *16 rootsall(150),rts(150)
        complex *16, allocatable :: rootslist(:)
c
c       BFS interval list (allocatable)
c
        parameter (maxstack=100 000,maxlevel=25,maxroots=10 000)
        complex *16, allocatable :: stack_c0(:)
        real *8, allocatable :: stack_sqw(:)
        integer, allocatable :: stack_ifd(:)
c
        complex *16 c0i,c0div(4),z0,val,dval
        complex *16 ima
        data ima/(0.0d0,1.0d0)/

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       input:
c       ifprint - 1 to print expansion coefficients
c       ifnewton - 1 to refine roots with Newton
c       ifres - 1 to count roots via residue
c       fun(z,par1,par2,val,dval) - user function,
c            must return val=f(z) and dval=f'(z)
c       par1, par2 - user parameters passed to fun
c       norder - expansion order (max 150, recommended 40)
c       eps - expansion accuracy (recommended 1d-4 to 1d-12)
c       nexdp - extra subdivision levels after convergence.
c            recommended 0: the convergence test is already
c            conservative (last 20% of coefficients < eps^(3/4)).
c            nexdp > 0 forces additional subdivision and can
c            help for difficult problems (e.g. roots near box
c            boundaries or high multiplicity).
c       center - center of the square (complex*16)
c       sqw - side length of the square
c
c       output:
c       nrtot - number of roots found
c       errest - error estimates for each root
c            if ifnewton=0: |f(root)|
c            if ifnewton=1: |f(root)/f'(root)|
c       roots - all roots found (complex*16, length nrtot)
c       centers - centers of all leaf boxes (length nc)
c       nc - number of leaf boxes
c       ier - error flag:
c            0    success
c            100  QR failure in one or more boxes
c            1024 max subdivision level reached
c            2048 max stack size exceeded
c            4096 max roots exceeded
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        allocate(rootslist(maxroots))
        allocate(stack_c0(maxstack),stack_sqw(maxstack))
        allocate(stack_ifd(maxstack))

c
c points on the square and weights
c
        n = norder
        call zrootsq_pts_wts(z,cw,n)

c
c load recurrence coefficients and evaluate polynomial basis
c
        call zrootsq_poly_load(q,vd,voffd,norder,z,4*n)

c
c precompute QR factorization of weighted basis matrix
c
        call zrootsq_factor(q,cw,n,norder,aqr,tau)

c
c initialize BFS
c
        stack_c0(1)=center
        stack_sqw(1)=sqw
        stack_ifd(1)=-nexdp

        ifirst=1
        ilast=1
        ilastnew=1
        nrtot=0
        nc=0
        nboxes=0
        ier=0

c
c BFS: process all boxes level by level
c
        do level=0,maxlevel

        iffini=1

        do icount=ifirst,ilast
          c0i=stack_c0(icount)
          sqwi=stack_sqw(icount)
          ifdone=stack_ifd(icount)

          call zrootsq0(eps,ifconv,fun,par1,par2,c0i,
     1        sqwi,vd,voffd,q,aqr,tau,z,cw,n,rootsall,norder,
     2        rts,nr,ierqr,ifprint)

          nboxes=nboxes+1

          if(ifconv.eq.1) ifdone=ifdone+1

          if(ierqr.gt.0) then
c
c QR failure: skip box, report
c
            ier=100
            write(6,*) 'QR failed at center=',c0i,' sqw=',sqwi
            write(13,*) 'QR failed at center=',c0i,' sqw=',sqwi

          elseif(ifdone.ge.1) then
c
c converged and done: collect roots
c
            if(nrtot+nr.gt.maxroots) then
              ier=4096
              goto 2000
            endif
            do i=1,nr
              rootslist(i+nrtot)=rts(i)
            enddo
            nrtot=nrtot+nr
            nc=nc+1
            centers(nc)=c0i

          else
c
c not done: subdivide into 4 children
c
            iffini=0
            call zrootsq_divide(c0i,sqwi,c0div,sqwdiv)

            do k=1,4
              ilastnew=ilastnew+1
              if(ilastnew.gt.maxstack) then
                ier=2048
                goto 2000
              endif
              stack_c0(ilastnew)=c0div(k)
              stack_sqw(ilastnew)=sqwdiv
              stack_ifd(ilastnew)=ifdone
            enddo
          endif

        enddo

        if(iffini.eq.1) goto 2000

        ifirst=ilast+1
        ilast=ilastnew

        enddo

        ier=1024

 2000   continue

        write(6,*) 'total box number is',nboxes
        write(13,*) 'total box number is',nboxes
        write(6,*) 'max subdivisions is',level
        write(13,*) 'max subdivisions is',level

c
c remove duplicated roots if subdivided
c
        if(nboxes.gt.1) then
          call zrootsq_dedup(eps,rootslist,nrtot,roots,nr)
        else
          nr=nrtot
          do i=1,nr
            roots(i)=rootslist(i)
          enddo
        endif

        nrtot=0
        do i=1,nr
          call fun(roots(i),par1,par2,val,dval)
          rootslist(i)=roots(i)
          nrtot=nrtot+1
          errest(i)=abs(val)
        enddo

        write(6,*) 'total root number found by the root finder',nrtot
        write(13,*) 'total root number found by the root finder',nrtot

c
c residue for total root count
c
        if(ifres.eq.1) then
          call zrootsq_residue(fun,par1,par2,nres,maxrec,
     1        center,sqw)
          write(6,*) 'total number of roots from residue is',nres
          write(13,*) 'total number of roots from residue is',nres

          ntor=5+log10(sqw)/log10(2d0)
          if(maxrec.gt.ntor) then
            write(6,*) 'Roots may be near boundary, maxrec=',maxrec
            write(13,*) 'Roots may be near boundary, maxrec=',maxrec
          endif
        endif

c
c refine with newton
c
        if(ifnewton.eq.1) then
          dzmax=-1
          do i=1,nrtot
            z0=rootslist(i)
            do ijk=1,2
              call fun(z0,par1,par2,val,dval)
              dzmag=abs(val/dval)
              if(dzmag.gt.dzmax) dzmax=dzmag
              z0=z0-val/dval
            enddo
            roots(i)=z0
          enddo
          write(6,*) 'max newton step size is',dzmax
          write(13,*) 'max newton step size is',dzmax

c
c update error after newton
c
          do i=1,nrtot
            call fun(roots(i),par1,par2,val,dval)
            errest(i)=abs(val/dval)
          enddo
        endif

        deallocate(rootslist,stack_c0,stack_sqw,stack_ifd)

        return
        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Divide a square into 4 children.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine zrootsq_divide(c0,sqw,c0div,sqwdiv)
        implicit real *8 (a-h,o-z)
        complex *16 c0,c0div(4),ima
        data ima/(0.0d0,1.0d0)/

        sqwdiv=sqw/2d0
        d=sqwdiv/2d0
        c0div(1)=c0+(-1d0+ima)*d
        c0div(2)=c0+(1d0+ima)*d
        c0div(3)=c0+(-1d0-ima)*d
        c0div(4)=c0+(1d0-ima)*d

        return
        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Remove duplicate roots from adjacent squares.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine zrootsq_dedup(eps,roots1,n1,roots2,n2)
        implicit real *8 (a-h,o-z)
        complex *16 roots1(*),roots2(*)

        roots2(1)=roots1(1)
        n2=1

        tol=sqrt(eps)*1d-2

        do i=1,n1

        ifput=1
        do j=1,n2
          d=abs(roots2(j)-roots1(i))
          if(d.lt.tol) ifput=0
        enddo

        if(ifput.eq.1) then
          n2=n2+1
          roots2(n2)=roots1(i)
        endif

        enddo

        return
        end


        subroutine zrootsq0(eps,ifcv,fun,par1,par2,center,
     1     sqw,vvd,vvoffd,q,aqr,tau,z,cw,n,rootsall,norder,
     2     roots,nr,ierqr,ifprint)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        complex *16 vvd(*),vvoffd(*),q(*)
        complex *16 aqr(*),tau(*)
        complex *16 z(*),cw(*),center
        complex *16 val,rootsall(*),roots(*)
c
c       local workspace on the stack
c
        complex *16 coefs(norder+5)
        complex *16 vd(norder+5),voffd(norder+5)
        complex *16 zz(4*n+5),cww(4*n+5)
        complex *16 rhs(4*n+5)
        real *8 coefabs(norder+5)
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Computes polynomial expansion of fun on a single square,
c       checks convergence, and if converged finds roots as
c       eigenvalues of the colleague matrix. Roots outside the
c       square (with tolerance dl) are discarded.
c
c       input:
c       eps, fun, norder, z, cw, n, center, sqw,
c       vvd, vvoffd (colleague matrix diagonal/off-diagonal),
c       q (basis matrix), aqr, tau (precomputed QR factorization)
c
c       output:
c       ifcv: 1 if converged, 0 otherwise
c       rootsall: all roots of the expansion polynomial
c       roots: roots inside the square
c       nr: number of roots inside the square
c       ierqr: -1 (not run), 0 (success), 1028 (QR failure)
cccccccccccccccccccccccccccccccccccccccccccccccccccc

        
        npt=4*n

        do i=1,norder
          vd(i)=vvd(i)
          voffd(i)=vvoffd(i)
        enddo

        ierqr=-1

c
c translate pts into correct domain
c
        do i=1,npt
          zz(i)=z(i)*sqw/(2d0)+center
          cww(i)=cw(i)*sqw/(2d0)
        enddo

c
c evaluate fun at quadrature points
c
c$omp parallel do private(val)
        do i=1,npt
          call fun(zz(i),par1,par2,rhs(i),val)
        enddo
c$omp end parallel do

c
c check for Inf/NaN in function values
c
        do i=1,npt
          if(abs(rhs(i)).ne.abs(rhs(i)) .or.
     1       abs(rhs(i)).gt.1d300) then
            ifcv=0
            goto 5000
          endif
        enddo

c
c find expansion coefs
c
        call zrootsq_leastsq1(rhs,cww,n,aqr,tau,
     1      norder,coefs)

c
c print expansion coefficients (absolute values)
c
        if(ifprint.eq.1) then
        do i=1,norder+1
          coefabs(i)=abs(coefs(i))
        enddo
        write(6,*) 'expansion coefs (absolute values):'
        write(6,'(5(1x,E12.5))') (coefabs(i),i=1,norder+1)
        write(13,*) 'expansion coefs (absolute values):'
        write(13,'(5(1x,E12.5))') (coefabs(i),i=1,norder+1)
        endif

c
c determine the largest |coefs(i)|
c
        coefsmax=-1
        do i=1,norder
          if(coefsmax.lt.abs(coefs(i))) coefsmax=abs(coefs(i))
        enddo

c
c check convergence: last 20% of coefficients should be below
c eps^(3/4) * coefsmax
c
        coefscut=0
        istart=norder*(4d0/5d0)
        do i=istart,norder+1
          if(abs(coefs(i)).gt.coefscut) coefscut=abs(coefs(i))
        enddo

        ifcv=1
        if(coefscut.gt.eps**(3d0/4d0)*coefsmax) then
          ifcv=0
          goto 5000
        endif

        call sparse_root_struct(rootsall,vd,voffd,norder,coefs,
     1    ierqr)

        do i=1,norder
          rootsall(i)=rootsall(i)*sqw/(2d0)+center
        enddo

c
c keep only roots inside the square (with tolerance dl)
c
        dl=1d-6
        nr=0
        do i=1,norder
          ifin=1
          if(abs(dreal(rootsall(i)-center)).ge.(sqw/(2d0)+dl*sqw))
     1        ifin=0
          if(abs(dimag(rootsall(i)-center)).ge.(sqw/(2d0)+dl*sqw))
     1        ifin=0

          if(ifin.eq.1) then
            nr=nr+1
            roots(nr)=rootsall(i)
          endif
        enddo

 5000 continue

        return
        end


        subroutine zrootsq_pts_wts(z,cw,n)
        implicit real *8 (a-h,o-z)
        complex *16 z(*),cw(*)
        dimension ts(10 000),whts(10 000)
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      Constructs Gaussian nodes z and weights cw on the four
c      sides of the unit square [-1,1]^2.  z and cw have length
c      4*n, where n is the number of points per side.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        

c
c Legendre nodes and weights on [-1,1]
c
        ifwhts=1
        call zroot_legewhts(n,ts,whts,ifwhts)

c
c scale weights for each side
c
        do i=1,n
          cw(i)=whts(i)
          cw(n+i)=whts(i)
          cw(2*n+i)=whts(i)
          cw(3*n+i)=whts(i)
        enddo

c
c nodes on the four sides
c
        do i=1,n
          z(i)=ts(i)*ima+1d0
          z(n+i)=-ts(i)+1d0*ima
          z(2*n+i)=-ts(i)*ima-1d0
          z(3*n+i)=ts(i)-1d0*ima
        enddo

        return
        end


        subroutine zrootsq_poly_load(q,vd,voffd,norder,z,npt)
        implicit real *8 (a-h,o-z)
        complex *16 vd(*),voffd(*)
        complex *16 q(npt,norder+1),z(*)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Loads three-term recurrence coefficients (vd, voffd)
c       and evaluates the orthogonal polynomial basis at the
c       quadrature points z(1:npt) via the recurrence:
c
c       P_0(z) = 1
c       voffd(k)*P_k(z) = (z-vd(k))*P_{k-1}(z) - voffd(k-1)*P_{k-2}(z)
c
c       Output: q(i,j) = P_{j-1}(z(i)) for i=1..npt, j=1..norder+1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        call diag_elem(vd,norder)
        call subdiag_elem(voffd,norder)

c
c P_0(z) = 1
c
        do i=1,npt
          q(i,1)=1
        enddo

c
c P_1(z) = (z - vd(1)) / voffd(1)
c
        do i=1,npt
          q(i,2)=(z(i)-vd(1))/voffd(1)
        enddo

c
c three-term recurrence for P_2 through P_norder
c
        do k=2,norder
          do i=1,npt
            q(i,k+1)=((z(i)-vd(k))*q(i,k)-voffd(k-1)*q(i,k-1))
     1              /voffd(k)
          enddo
        enddo

        return
        end




        subroutine zrootsq_factor(q,cw,n,norder,aqr,tau)
        implicit real *8 (a-h,o-z)
        complex *16 q(4*n,norder+1),cw(*)
        complex *16 aqr(4*n,norder+1),tau(*)
        complex *16 zwork(100 000)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Precomputes the QR factorization of the weighted basis matrix
c       aqr(i,j) = sqrt(cw(i)) * q(i,j)
c using LAPACK zgeqrf. The factored matrix aqr and the
c Householder reflectors tau are stored for later use by
c zrootsq_leastsq1.
c
c This is called once; the factorization is reused for every box.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        npt=4*n
        m=norder+1

c
c form weighted basis matrix aqr(i,j) = sqrt(cw(i)) * q(i,j)
c
        do j=1,m
          do i=1,npt
            aqr(i,j)=q(i,j)*sqrt(cw(i))
          enddo
        enddo

c
c QR factorize: aqr = Q * R (stored in-place)
c
        lwork=100 000
        call zgeqrf(npt,m,aqr,npt,tau,zwork,lwork,info)

        return
        end


        subroutine zrootsq_leastsq1(fvals,cw,n,aqr,tau,
     1      norder,coeff)
        implicit real *8 (a-h,o-z)
        complex *16 aqr(4*n,norder+1),tau(*),coeff(*)
        complex *16 fvals(*),cw(*)
        complex *16 rhs(4*n)
        complex *16 zwork(100 000)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Computes expansion coefficients from function values fvals
c by solving the weighted least squares problem
c       min ||sqrt(cw) * (fvals - q*coeff)||
c using the precomputed QR factorization (aqr, tau) from
c zrootsq_factor. Applies Q^H to the rhs, then solves R*x = Q^H*b
c via back-substitution.
c
c fvals are the function values at the quadrature points
c (evaluated by the caller, not by this subroutine).
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        npt=4*n
        m=norder+1

c
c form weighted rhs
c
        do i=1,npt
          rhs(i)=sqrt(cw(i))*fvals(i)
        enddo

c
c apply Q^H to rhs: rhs <- Q^H * rhs
c
        lwork=100 000
        call zunmqr('L','C',npt,1,m,aqr,npt,tau,rhs,npt,
     1      zwork,lwork,info)

c
c solve R * coeff = rhs(1:m) by back-substitution
c
        call ztrtrs('U','N','N',m,1,aqr,npt,rhs,npt,info)

        do i=1,m
          coeff(i)=rhs(i)
        enddo

        return
        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Residue-based root counting via contour integration of f'/f.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine zrootsq_residue(fun,upar1,upar2,nres,mrec,
     1      center,sqw)
        implicit real *8 (a-h,o-z)
        dimension upar1(*),upar2(*)
        complex *16 rint1,rint2,rint3,rint4,rint
        complex *16 center,gpar(2)
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
        external zrootsq_residue_int
c
c count roots by integrating f'/f on the four sides
c
        
        b=sqw/(2d0)
        a=-b

        gpar(1)=center
        gpar(2)=sqw

        eps=sqrt(epsilon(1d0))
        m=32

        mrec=0
        do iside=1,4
          call zroot_adapgau(ier,a,b,zrootsq_residue_int,fun,
     1        upar1,upar2,gpar,iside,m,eps,rint1,maxrec,numint)
          if(ier.ne.0) write(6,*) 'C',iside,' int failure, ier=',ier
          if(mrec.lt.maxrec) mrec=maxrec

          if(iside.eq.1) rint=rint1*ima
          if(iside.eq.2) rint=rint-rint1
          if(iside.eq.3) rint=rint-rint1*ima
          if(iside.eq.4) rint=rint+rint1
        enddo

        pi=4*atan(1d0)
        nres=dimag(rint/(2*pi))+0.1

        return
        end


        subroutine zrootsq_residue_int(x,fun,upar1,upar2,gpar,
     1      iside,val)
        implicit real *8 (a-h,o-z)
        dimension upar1(*),upar2(*)
        complex *16 val,val1,val2,z
        complex *16 gpar(*),c0
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
c       Computes f'(z)/f(z) at a point z on side iside of the square,
c       parameterized by x in [-sqw/2, sqw/2].
c
c       iside=1: right  side, z = sqw/2 + i*x + c0
c       iside=2: top    side, z = x + i*sqw/2 + c0
c       iside=3: left   side, z = -sqw/2 + i*x + c0
c       iside=4: bottom side, z = x - i*sqw/2 + c0
c
        
        sqw=dreal(gpar(2))/(2d0)
        c0=gpar(1)

        if(iside.eq.1) z=(sqw+ima*x)+c0
        if(iside.eq.2) z=(x+ima*sqw)+c0
        if(iside.eq.3) z=(-sqw+ima*x)+c0
        if(iside.eq.4) z=(x-ima*sqw)+c0

        call fun(z,upar1,upar2,val1,val2)
        val=val2/val1

        return
        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Adaptive Gaussian quadrature for complex-valued integrands.
c       Modified from cadapgau_new to pass an integer parameter
c       iside through to the integrand callback.
c
c       Callback signature: fun(x, par1, par2, iside, val)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine zroot_adapgau(ier,a,b,callback,userfun,
     1      upar1,upar2,gpar,iside,m,eps,rint,maxrec,numint)
        implicit real *8 (a-h,o-z)
        complex *16 rint,vals(200)
        dimension t(100),w(100),stack(400),
     1      upar1(*),upar2(*),gpar(*)

        ifwhts=1
        call zroot_legewhts(m,t,w,ifwhts)

        nnmax=100000
        maxdepth=200

        call zroot_adinrec(ier,stack,a,b,callback,userfun,
     1      upar1,upar2,gpar,
     2      iside,t,w,m,vals,nnmax,eps,rint,maxdepth,maxrec,numint)
        return
        end


        subroutine zroot_adinrec(ier,stack,a,b,callback,
     1      userfun,upar1,upar2,gpar,iside,t,w,m,vals,nnmax,eps,
     2      rint,maxdepth,maxrec,numint)
        implicit real *8 (a-h,o-z)
        complex *16 rint,value2,value3,vals(*)
        dimension stack(2,*),t(*),w(*),upar1(*),upar2(*),gpar(*)

        stack(1,1)=a
        stack(2,1)=b
        call zroot_oneint(a,b,callback,userfun,upar1,upar2,gpar,
     1      iside,t,w,m,vals(1))

        j=1
        rint=0
        ier=0
        maxrec=0
        do 3000 i=1,nnmax
        numint=i
        if(j .gt. maxrec) maxrec=j

         c=(stack(1,j)+stack(2,j))/2
        call zroot_oneint(stack(1,j),c,callback,userfun,
     1      upar1,upar2,gpar,iside,t,w,m,value2)

        call zroot_oneint(c,stack(2,j),callback,userfun,
     1      upar1,upar2,gpar,iside,t,w,m,value3)

        dd=cdabs(value2+value3-vals(j))
        ifdone=0
        if(dd .le. eps) ifdone=1

        if(ifdone  .eq. 0) goto 2000

        rint=rint+value2+value3
        j=j-1

        if(j .eq. 0) return
        goto 3000
 2000 continue

        stack(1,j+1)=stack(1,j)
        stack(2,j+1)=(stack(1,j)+stack(2,j))/2
        vals(j+1)=value2

        stack(1,j)=(stack(1,j)+stack(2,j))/2
        vals(j)=value3

        j=j+1

        if(j .le. maxdepth) goto 3000
        ier=8
        return
 3000 continue
        ier=16
        return
        end


        subroutine zroot_oneint(a,b,callback,userfun,upar1,upar2,
     1      gpar,iside,t,w,m,rint)
        implicit real *8 (a-h,o-z)
        dimension t(*),w(*),upar1(*),upar2(*),gpar(*)
        complex *16 rint,val

        rint=0
        u=(b-a)/2
        v=(b+a)/2
        do i=1,m
          tt=u*t(i)+v
          call callback(tt,userfun,upar1,upar2,gpar,iside,val)
          rint=rint+val*w(i)
        enddo
        rint=rint*u
        return
        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Gauss-Legendre quadrature nodes and weights.
c       Extracted from legeexps.f.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine zroot_legewhts(n,ts,whts,ifwhts)
        implicit real *8 (a-h,o-z)
        dimension ts(*),whts(*)
c
c        Constructs the nodes and the weights of the n-point
c        Gaussian quadrature on the interval [-1,1].
c
c        input:  n - number of nodes
c        output: ts - nodes, whts - weights
c        ifwhts: 0 = nodes only, 1 = nodes and weights
c
        eps=1.0d-14
        
        pi=datan(1d0)*4
        h=pi/(2*n)
        do i=1,n
          t=(2*i-1)*h
          ts(n-i+1)=dcos(t)
        enddo

        ts(n/2+1)=0
        do i=1,n/2
          xk=ts(i)
          ifout=0
          do k=1,10
            call zroot_legepol_sum(xk,n,pol,der,sum)
            delta=-pol/der
            xk=xk+delta
            if(abs(delta) .lt. eps) ifout=ifout+1
            if(ifout .eq. 3) goto 1600
          enddo
 1600     continue
          ts(i)=xk
          ts(n-i+1)=-xk
        enddo

        if(ifwhts .eq. 0) return

        do i=1,(n+1)/2
          call zroot_legepol_sum(ts(i),n,pol,der,sum)
          whts(i)=1/sum
          whts(n-i+1)=whts(i)
        enddo

        return
        end


        subroutine zroot_legepol_sum(x,n,pol,der,sum)
        implicit real *8 (a-h,o-z)
        
        sum=0

        pk=1
        pkp1=x
        sum=sum+pk**2/2
        sum=sum+pkp1**2*(1+1d0/2)

        if(n .ge. 2) goto 1200

        sum=0
        pol=1
        der=0
        sum=sum+pol**2/2
        if(n .eq. 0) return

        pol=x
        der=1
        sum=sum+pol**2*(1+1d0/2)
        return
 1200 continue

        do k=1,n-1
          pkm1=pk
          pk=pkp1
          pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
          sum=sum+pkp1**2*(k+1+1d0/2)
        enddo

        pol=pkp1
        der=n*(x*pkp1-pk)/(x**2-1)
        return
        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Precomputed three-term recurrence coefficients for the
c       orthogonal polynomial basis (diagonal and off-diagonal).
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine diag_elem(a,n)
c 
c n is the order of the expansion
c
c a is the vector in the main diagonal and has length n        
c from index=0 to index=n-1
c
        implicit real *8 (a-h,o-z)
        real *8 a(*)
        real *8 b(2*100)
        save 
        data b/
     1    0.3036776984963890719648606804639081D-01,
     2    0.1030196785487212514705985765443933D-02,
     3   -0.2263831926188826289124021679549482D+01,
     4    0.7363266794523236243590459649171935D+00,
     5    0.4274784196977596793656277744964483D+01,
     6   -0.1309141514110488631424179980299076D+01,
     7   -0.1970520023567703072598214907953918D+01,
     8    0.6209215121189803436585337622647937D+00,
     9   -0.1138239117917704292356739019613377D+00,
     *    0.1216471945326044633572477494040142D-01,
     1   -0.1109289680457407418179441520235096D+01,
     2    0.8738259182999032623222451837493250D+01,
     3    0.4521795263543255810612382611407837D+01,
     4   -0.1171480044837908212225188930591947D+02,
     5   -0.3112636421462545873231251870884133D+01,
     6    0.2703021121142773960797733850490684D+01,
     7   -0.1099999912612343372505357102901803D+01,
     8    0.1844286345574250863627077433860568D+00,
     9    0.3073225798711069915325427071047283D+01,
     *    0.1092837699234445739575330983973550D+01,
     1    0.1225605490629735142186991988536894D+01,
     2   -0.1295895495389606652370279775289358D+01,
     3   -0.2409788162918876550868634162352033D+01,
     4   -0.2453872804200626382507427300519134D+01,
     5   -0.7555985663251584151541089890630731D+00,
     6    0.2600689200755379488056526002755205D+01,
     7   -0.2048930238538691565247235996189770D+01,
     8   -0.3292168277772834312237925518654082D+00,
     9    0.2462039989397908806277052710808911D+01,
     *    0.1130552052943225565255058428017123D+01,
     1   -0.1197910341617292552541867855061562D+01,
     2   -0.9123905486276858287956178051513854D+00,
     3    0.1062562478773726318975320478750696D+01,
     4    0.4867560638192171116479637982277501D+00,
     5   -0.1379413990219107470165879852801873D+01,
     6    0.1123868532599553989336245053566958D+01,
     7   -0.3095801227046912474426325131885125D-01,
     8   -0.2916943458135204321976097119714363D+00,
     9    0.1923593664123302519638500331386618D+01,
     *   -0.5583457278897859207843027207900583D+00,
     1   -0.1380483640978667242712477032068018D+01,
     2    0.7742852504604716681514311177772121D+00,
     3   -0.5008883087246050827170090494961545D+00,
     4   -0.2276715672522905105304792982654296D+01,
     5    0.1384539656618468973263308580541730D+01,
     6    0.1302717514246881006871615402743845D+01,
     7   -0.1342339479048226165677317743789184D+01,
     8   -0.4489022173377686908860689859400316D+00,
     9   -0.1313401361707819575198599769581398D+00,
     *   -0.1303056852897112393287446086155597D+01,
     1    0.1182456458961619029110598932043124D+01,
     2    0.5590438321097331764046589335681336D+00,
     3    0.3410515447586137960351244703589851D+00,
     4    0.1096791524441345956341588687472621D+01,
     5   -0.8688618276138566006689761254249336D+00,
     6   -0.8945591837178112955640649710865514D-01,
     7    0.7817267277457608029726594650485875D+00,
     8   -0.2537121891911355292556493725573726D+00,
     9   -0.1306451213032109536825273179010231D+00,
     *   -0.4788887149714040774024380838381513D-01,
     1   -0.1255702884327540569205089342087624D+01,
     2    0.4118460166247315629098104403409304D-02,
     3    0.1398187256113509982165350495475637D+01,
     4    0.1183506688357217851392843206104212D+01,
     5   -0.1264359257183476007874053053287617D+01,
     6   -0.1469244218412277187825721457545656D+00,
     7    0.6744767074455227863341414405869985D+00,
     8   -0.1124118118906685047818266304416692D+01,
     9    0.1213217083497061303825959144155277D-01,
     *    0.1920944459370085807861466252526035D+01,
     1   -0.2015102301093099816831330389373691D+01,
     2   -0.1068279875884921550263762281815137D+01,
     3    0.1723539209364344605853523347503758D+01,
     4   -0.9122742815670968713115765841472152D+00,
     5    0.3038158208168553397499226073288828D+00,
     6    0.6025474996807269043864978091579641D+00,
     7   -0.9511771862035712653403092890019526D+01,
     8   -0.1115415403004397325450416094751441D+01,
     9    0.1064450272664890704298979906605913D+02,
     *   -0.2488575438521785297959483578757373D+00,
     1   -0.6935463520848694096396354095811037D+00,
     2    0.1808665289464023531753898149058460D+01,
     3   -0.4888441687598371868394144070199117D-01,
     4   -0.2128992348443366041066879168919758D+01,
     5   -0.1512965923570476467972748382530732D+00,
     6    0.1369001235173170313249256757367722D+01,
     7    0.2393918595733830524882692709364755D+01,
     8   -0.3960015108115400581065804395655033D+01,
     9   -0.2234515581464662607844576401923207D+01,
     *    0.4495586411431224235425522866844662D+01,
     1   -0.9718402456798038790720777288638035D+00,
     2   -0.1528410138812708664468435714171426D+01,
     3    0.1417392608071106122033045493684528D+01,
     4    0.6529441613781810853468870332998318D+00,
     5   -0.1246742960120980033412133209003466D+01,
     6   -0.3494198724584419930598222584889347D+00,
     7    0.1806558601001270476326315555154864D+01,
     8    0.1132471351107087780920382557876574D+01,
     9   -0.1385967205320723252969419731442502D+01,
     *    0.1067515714529940075553353519105179D+01,
     1    0.2059082442960325120379639426441890D+00,
     2    0.5214759843463801080462922132213995D+00,
     3   -0.1510801488992807916562464675163206D+01,
     4   -0.2685962630412040454186009542041389D+01,
     5    0.7678040228856816973801528082138248D+00,
     6    0.1459323840471671913546427411178064D+01,
     7   -0.3141730129294775918894032155310885D+00,
     8   -0.9362777773754649874852615385110986D+00,
     9    0.1336086952713150652418327676352879D+01,
     *    0.1500022896661031675098411926976107D+01,
     1   -0.5297479896671643007171218878648387D+00,
     2   -0.3352226059238493139984572515856174D+00,
     3    0.6888754524513588686517168231293673D+00,
     4   -0.1021112681321268109009335050762759D+01,
     5   -0.1205275851276851337299899014628440D+00,
     6    0.7536825259689853205455808817748984D+00,
     7   -0.1346849580451283942923330348171103D+01,
     8   -0.1286092107692669303213118320139717D+01,
     9    0.2245540132387460905109843394331260D+00,
     *    0.3302809999905728408917196934345695D+00,
     1    0.1498313747262144279422159739082461D+00,
     2   -0.2884726137386325269250322633861842D+00,
     3    0.4867025478649832949111989805558363D+00,
     4    0.2588327430513351120436823896181191D+01,
     5    0.5011885314913070002392586634543772D-01,
     6   -0.1868733981555842549537928070058611D+01,
     7   -0.2344827553048175917217523577088683D+01,
     8    0.5884352899215413683182806073166353D+00,
     9    0.9130102014848219407037785258338958D+00,
     *   -0.1361484858905212517008346334681233D+01,
     1    0.2219723079897351130505418682592104D+01,
     2    0.4009272272678185086811742960584377D+00,
     3    0.2715286939871629970390446855738645D+00,
     4    0.1636114045238399858977409652323124D+01,
     5   -0.7306359952409286820265016682738024D+00,
     6   -0.2126371299895227026761293072784298D+01,
     7   -0.1093494187214921767118448286803280D+01,
     8    0.1497158578261318495840110707108415D+01,
     9    0.3649452053587675025978085232824911D+00,
     *    0.4309785993912294949234018079996853D+00,
     1    0.3524383131449079230832803876945821D+00,
     2   -0.7654446046629473675917429837232586D+00,
     3   -0.2655387078804165525088029912398489D+00,
     4   -0.8607646642455116594829126414434241D+00,
     5   -0.2245044262538080098594354527638477D+00,
     6   -0.6888795005449932763120893632668534D+00,
     7   -0.5103982779916414169025046217213682D+01,
     8    0.5605688250570435131922099196945735D+01,
     9    0.5453894464785898888941663158876859D+01,
     *   -0.5480479492482228994681407139115602D+01,
     1   -0.3408315770770519042380702001041131D+01,
     2    0.1406485320242628111242125893084713D+01,
     3    0.4164112362677757221155674156376626D+01,
     4   -0.1293015758956267362027320397652979D+01,
     5    0.2217683679959967441865980137982388D+00,
     6   -0.2495301775306532456279050620814569D+01,
     7   -0.8073666194619870599524785897488860D+00,
     8    0.2966031421559124064917503073500114D+01,
     9    0.5525303303447801713544256576393496D+00,
     *    0.3322201772471715278489463565266556D+00,
     1   -0.1066766834824709394309848013074905D+01,
     2    0.3530247673698011410646404325209052D+00,
     3    0.1810054534439419207085979252636124D+01,
     4   -0.1135726084194808501722757882002050D+01,
     5   -0.8751204509997197611798125933629023D+00,
     6   -0.1119433514697289738153584845731421D+01,
     7   -0.1368342880684185610261691245130986D+01,
     8    0.1260464191931636642725505911722078D+01,
     9    0.1392116715094043456476986697683124D+01,
     *    0.2660469823327112570224178730691291D-01,
     1   -0.4038433692015993290360621043849076D+00,
     2    0.1655527362495244311099560939785704D+00,
     3   -0.5664716776282386118822509457432331D+01,
     4   -0.1238070440148386194079330935930206D+01,
     5    0.4956200455191687026649523112214894D+01,
     6    0.1405432830021015537673114747256375D+01,
     7   -0.1805147207552500285517679399970688D+01,
     8    0.2248468656978182600931328144948788D+00,
     9    0.2001842910107993337022983372463778D+01,
     *   -0.3297959822798429296277311972471277D+00,
     1    0.3734960226840081025300645371553058D+01,
     2   -0.3254666100292063453452001718383957D+01,
     3   -0.3705185478843245144583381334076763D+01,
     4    0.3157005522465866659031874423482300D+01,
     5    0.4869897259549685799790843553474425D+00,
     6   -0.3868911409679611780479589821626428D+00,
     7   -0.6576355184234647894770489851479289D+00,
     8    0.2239911848961670875975236308715238D+00,
     9   -0.5393925595737874054478852092621873D+01,
     *    0.7508796617199647510486448925576017D+00,
     1    0.6420589778526103684230005146605683D+01,
     2   -0.8547293734328067561124467285775607D+00,
     3   -0.1345209585750250295380287128805357D+01,
     4    0.4300822629707437552231064306425613D+00,
     5   -0.6723772644395028321566104614033330D-01,
     6    0.1385059281432101324082047826565999D+01,
     7    0.1415026220028163895619400461374435D+01,
     8   -0.2638215640462922448973590603884720D+01,
     9    0.6775630068374320404832165411033024D+00,
     *    0.2696977664675180035186220647489082D+01/





        do i=1,2*n
          a(i)=b(i)
        enddo

        return
        end
        
        subroutine subdiag_elem(a,n)
c 
c n is the order of the expansion
c
c The vector a is in the sup/subdiagonal and has length n        
c from index=0 to index=n-1
c
c From 0 to n-2 goes to the sup/subdiagonal of the colleague matrix
c n-1 is used in the rank one part

        implicit real *8 (a-h,o-z)
        real *8 a(*)
        real *8 b(2*100)
        save 
        data b/
     1    0.1923495204367377238112581268734139D+00,
     2    0.5203057485431622924639683675445002D-01,
     3    0.1705565363979569071320835043663309D+01,
     4    0.6338155740481245583728099264856847D+01,
     5    0.5152198481923434649505716470148896D+01,
     6   -0.1279129040773978969933997449086109D+01,
     7    0.5674229393367071797141091540760772D-01,
     8   -0.1944827566759655014554418972549039D+00,
     9    0.9823003975496773909889618700516142D-01,
     *    0.9277971471521227617414035658077208D-01,
     1    0.9852565132165426638966259786109800D+01,
     2    0.3778099847959403120132755923562080D+01,
     3    0.3484824931667706817341473796477764D+01,
     4   -0.4450675682889525811026629514796543D+00,
     5    0.5111350888807004028455451540982299D-01,
     6   -0.1721716124377019594205149061445525D+00,
     7    0.2560826835244521023347032037422735D+00,
     8   -0.6693961125609358576925468900584187D+00,
     9    0.1064508650550630002440160471388471D+01,
     *   -0.1480736016564448840777254891997209D+01,
     1    0.4504538387296152367555549877477558D+00,
     2   -0.1914244136540372857619850566082410D+01,
     3    0.1932613761257918864659971804394271D+01,
     4   -0.3075082374001163087630008206771555D+00,
     5    0.2563529612817569277993906450180676D+00,
     6   -0.1673630293210639747554562155736408D-01,
     7    0.1387540743958712093538320772246432D+01,
     8   -0.2178252178034106837646802199498505D+01,
     9    0.4204322373694631312076650771055014D+00,
     *    0.1815639751360799894531531475964947D+01,
     1    0.9661702332026339100066871953851090D+00,
     2    0.4757553519403311427893049308997992D+00,
     3    0.3757689796949451290100935213231637D+00,
     4   -0.7837320707272966792906683540869793D+00,
     5    0.9037965977083404029655591209963919D+00,
     6    0.8892790946070312584960343015505839D+00,
     7    0.5057023667244722466212001717261307D+00,
     8    0.8944623788489772149270212854974374D+00,
     9    0.1672820713409693680880234782565098D+00,
     *   -0.1095072396494123656421541817457585D+01,
     1    0.1540857117021665357100366279928892D+01,
     2    0.3634902691921575729412115634633821D+00,
     3    0.8834384727891074928561989862005162D-01,
     4   -0.1157033109157053145634733631671541D+01,
     5    0.1063351962685883088251919884865472D+01,
     6   -0.1246344904133967752057578923325773D+01,
     7    0.3695593092627493513208771193449963D+00,
     8    0.2806450446923304144673256325988179D+00,
     9    0.1369034044429016778499666953035782D+01,
     *   -0.8360892271463443947994804118509994D+00,
     1    0.8750849565779075195111509997143016D+00,
     2    0.8838102000361465201204272872862943D+00,
     3    0.4344086519630630207353278008046807D+00,
     4   -0.1321420194141831021402415386384249D+01,
     5    0.9776036788076114403616667930413652D+00,
     6   -0.6609765706523619971980016360926125D+00,
     7    0.1195319300733472034162429437813271D+01,
     8    0.9684487657584631756763713219991206D+00,
     9    0.5865627065383242165743212539343851D-01,
     *   -0.9032905583050624878615724082394392D+00,
     1    0.8884990828852397391766753012379649D+00,
     2   -0.4984202096705334110956627617785496D+00,
     3    0.5503895548505834190813982636509196D+00,
     4   -0.1467417279478730376430592410909778D+01,
     5    0.7022909624392676719816061150857339D+00,
     6    0.2355457870926016599420037913427247D+00,
     7    0.5277593316275356998779251388623165D-01,
     8    0.7330254101392849356286376845199480D+00,
     9    0.1778195141686376296663816753163242D+01,
     *   -0.9275635731402565611410182070715324D+00,
     1    0.9838505032763365763654931365757430D+00,
     2    0.1464210925123777722878721203549511D+01,
     3    0.2955630763956252413888382830183571D+00,
     4   -0.1186511428195767787088160656072478D+01,
     5    0.2373410992985849820113029810526210D+00,
     6    0.2841330162721890923315891539538513D-01,
     7    0.3937285321290536488748100957397366D+00,
     8   -0.9967258901463986130098821785512026D+01,
     9    0.3773179904830828665141025646704289D+00,
     *    0.7481059356693873952656118433121516D+00,
     1    0.7331073635371181062135250138600962D+00,
     2   -0.7903626517872721847877585895171896D-01,
     3    0.1689413407192417967424516121384873D+01,
     4    0.2579154799836399669939948524241539D+00,
     5    0.1812749168183230912963260196427466D+00,
     6   -0.2856197543553746209292673429281464D+00,
     7    0.4276292687650511861146287253311796D+01,
     8    0.2099590007175576611366379810626569D+01,
     9    0.7503798283984029950545237455235540D+00,
     *    0.6049955215823594751509222781267972D+00,
     1    0.1640393280190247468660540206876447D+00,
     2   -0.9504615533436039322085349454058981D+00,
     3    0.6668395488600469629882273797819978D+00,
     4   -0.6199899189090730386004110075969955D-02,
     5    0.3507997517942377042918010909154722D+00,
     6   -0.2213139928238165031701035200186545D+01,
     7    0.1192357837836846208054951940329545D+01,
     8   -0.7394891766393726380679797180342852D+00,
     9    0.9025262611583402356448236932594012D+00,
     *    0.8707666350680245918354430537603448D-01,
     1    0.1741980980481635586249089114770927D+01,
     2   -0.8242631991381198424925716558795022D+00,
     3    0.1019646544646454484692322700870795D+01,
     4   -0.1698297977681668440631742933822228D+00,
     5    0.1120030247834702896175167807596235D+01,
     6   -0.1067042984177612844782139741859444D+01,
     7    0.6812058412047538607239849990053829D+00,
     8    0.6926678736040007582403375394969102D+00,
     9    0.8343117839275379140682779540610367D+00,
     *   -0.7337235524184659705744529985459219D+00,
     1    0.9424215001997809030417350074607274D+00,
     2    0.4680699185924130475477771205420260D+00,
     3    0.2942543510049305599867362652368230D+00,
     4   -0.7843244790672838958919315492177791D+00,
     5    0.1284323380034938780770125407783293D+01,
     6   -0.1266311306975735964134862411524061D+01,
     7    0.1139069640903432122294432881097467D+01,
     8    0.6404877541927275698699496597224000D+00,
     9    0.6590273562450003710612793409564621D+00,
     *   -0.5573403429802879309711411101832009D+00,
     1    0.1439845928037065786532766656800811D+00,
     2    0.7227010060762646846548402031619742D+00,
     3    0.2325077477751531522476512542909652D+01,
     4   -0.1636877949474514264335818278104618D+00,
     5    0.2250726394436435407247829425413159D+00,
     6    0.5759432103358482107008182737891473D+00,
     7    0.7204092411099942062611650496462265D+00,
     8    0.1779071591494357975008933452458125D+01,
     9    0.3811393362018149143618127931053338D+00,
     *   -0.8081850357578361872870352142938614D+00,
     1    0.1556077835253322069943249301064918D+00,
     2    0.8436934486995830553871571154549411D+00,
     3    0.1640146728608238987847163710601598D+01,
     4   -0.6743770342196251591268772694487992D+00,
     5    0.6247060750309549984170713073563716D+00,
     6    0.5522124973147225458598702957249981D+00,
     7    0.1008067530119974622678633752633641D+01,
     8    0.6062249042992685538890762172131900D+00,
     9    0.8938776046007852870210384807036293D+00,
     *    0.1200202747846392491000687316293717D+00,
     1    0.1400082684779677447480055261192511D+01,
     2   -0.3799406485545469521148016578262382D+00,
     3    0.3903975205809937553197793151127974D+00,
     4    0.1443916964902853238919542525910712D+01,
     5    0.6433724593287331172523935207001767D+00,
     6    0.1250245414157346180687826461791742D+00,
     7    0.5426997583240463952443697182649372D+01,
     8    0.5184658386956210055897003375717946D+01,
     9    0.2970176759404348182381169013588321D+00,
     *   -0.9910683159572319978616203897198209D-02,
     1    0.1380363057140953028028153768664330D+01,
     2    0.3782203443492573709729898176883237D+01,
     3    0.3785105859651342667442488965260882D+00,
     4   -0.1034748099881203899473186831344316D+00,
     5    0.2686925068635927460855427185627015D+01,
     6    0.5131437627913442837822584587390723D+00,
     7    0.5188892441548270273525801145931631D+00,
     8   -0.1364585770325727871028851118168587D+00,
     9    0.7564728419238961588065602816191400D+00,
     *    0.9055453186815995613895350151226593D+00,
     1    0.6036436692904591382208013856165061D+00,
     2   -0.8572636466825099948504114325189405D+00,
     3    0.5333534701588809430858124088541242D+00,
     4    0.1430631643362942278144751267943792D+01,
     5    0.8351406132536766849550802919690700D+00,
     6    0.6274255654581418210947446344471079D+00,
     7    0.2939826027387104968302255843893356D+00,
     8    0.1180525149349258843888910199013807D+01,
     9    0.8407570045947384860009465949324929D+00,
     *    0.6938768225067462363540559444861161D+00,
     1    0.3478351274974918441560276814710451D+00,
     2    0.5450122169470188494083577165706704D+00,
     3    0.1446837084025495808682106917080357D+01,
     4   -0.5266810973676843087983417902190729D+01,
     5    0.1691989033769279824971880240539182D+00,
     6    0.3883333533525291074959839263583391D+00,
     7    0.3890432881470821909866057838680904D+00,
     8    0.2139948730112129684451484942170339D+01,
     9    0.4846377885530350197035589222187354D-01,
     *    0.5483212966432591367655200423335543D+00,
     1    0.3129986535613952910165782719615112D+01,
     2    0.3544795701577519391256122488883616D+01,
     3    0.4420612501281071426339364730901120D-01,
     4    0.7151942124070659374574005980316453D+00,
     5    0.9947844757253360781617425804648552D+00,
     6    0.8499897566020296070025567996418620D+00,
     7    0.3946866836005982293708845088059620D+00,
     8   -0.1296697626882432486345878911477400D+00,
     9    0.7835164219179124998822773566292220D+00,
     *    0.6023886542287286472118000716014290D+01,
     1    0.7209067162227893557409903335308733D+00,
     2   -0.3128311574475158257762956687817291D-01,
     3    0.3385185985192104566019062082377774D+00,
     4    0.4728779346446554383038847021253540D+00,
     5    0.8148972807298324326548144775549758D+00,
     6   -0.6251662559792881577089934878400566D+00,
     7    0.2506562485165803891231484258930403D+01,
     8    0.7849271840879905713156301661682678D+00,
     9    0.9955940555710414273371857493976524D+00,
     *   -0.8883235107109528772870712101986674D+00/





        do i=1,2*n
          a(i)=b(i)
        enddo

        return
        end

