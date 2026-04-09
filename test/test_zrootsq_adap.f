c
c     Test driver for zrootsq_adap adaptive complex rootfinder.
c     Tests several functions requiring quadtree subdivision.
c
        implicit real *8 (a-h,o-z)
        real *8 errest(10 000),par1(1),par2(1)
        complex *16 center,ima
        complex *16 roots(10 000)
        complex *16, allocatable :: centers(:)
        data ima/(0.0d0,1.0d0)/

        external test_fun1
        external test_fun3
        external test_fun6
        external test_fun9
        external test_fun_pole1
        external test_fun_pole2
        external test_fun_bigsin

        allocate(centers(10 000 000))
        call prini(6,13)

        ifnewton=1
c       ifnewton=0

c
c === Test 1: polynomial with multiplicity ===
c 14 roots in [-1.5,1.5]^2, single box should suffice
c tests that adaptive wrapper works on easy cases
c
        print *,'=== test_fun1: polynomial, adap ==='
        norder=40
        center=0d0
        sqw=3d0
        eps=1d-14
        nexdp=1

        call cpu_time(t1)
!$      t1=omp_get_wtime()
        call zrootsq_adap(0,ifnewton,1,
     1      test_fun1,par1,par2,norder,eps,nexdp,center,sqw,
     2      nrtot,errest,roots,centers,nc,ier)
        call cpu_time(t2)
!$      t2=omp_get_wtime()
        call prinf('nrtot *',nrtot,1)
        call prinf('ier *',ier,1)
        write(6,*) 'time (sec): ',t2-t1
        print *,''

c
c === Test 2: sin(3*pi*x/2)/(x-2), 3 roots ===
c search region centered at 1, side 2.13
c expect 3 roots: 0, 2/3, 4/3
c
        print *,'=== test_fun3: sin, 3 roots, adap ==='
        norder=40
        center=1d0
        sqw=2.13d0
        eps=1d-10
        nexdp=1

        call cpu_time(t1)
!$      t1=omp_get_wtime()
        call zrootsq_adap(0,ifnewton,1,
     1      test_fun3,par1,par2,norder,eps,nexdp,center,sqw,
     2      nrtot,errest,roots,centers,nc,ier)
        call cpu_time(t2)
!$      t2=omp_get_wtime()
        call prinf('nrtot *',nrtot,1)
        call prinf('ier *',ier,1)
        write(6,*) 'time (sec): ',t2-t1
        print *,''

c
c === Test 3: sin(100/(x-2)), many roots near x=2 ===
c search region centered at 0, side 2.13
c many roots accumulating near singularity at x=2
c
        print *,'=== test_fun6: sin(100/(x-2)), adap ==='
        norder=40
        center=0d0
        sqw=2.13d0
        eps=1d-10
        nexdp=1

        call cpu_time(t1)
!$      t1=omp_get_wtime()
        call zrootsq_adap(0,ifnewton,1,
     1      test_fun6,par1,par2,norder,eps,nexdp,center,sqw,
     2      nrtot,errest,roots,centers,nc,ier)
        call cpu_time(t2)
!$      t2=omp_get_wtime()
        call prinf('nrtot *',nrtot,1)
        call prinf('nc (leaf boxes) *',nc,1)
        call prinf('ier *',ier,1)
        write(6,*) 'time (sec): ',t2-t1
        print *,''

c
c === Test 4: test_fun9, 104 roots (main adaptive test) ===
c oscillatory sin(d/(x*r0-x0)) with d=50, rotation
c search region centered at 0, side 2.13
c
        print *,'=== test_fun9: 104 roots, adap ==='
        norder=40
        center=0d0
        sqw=2.13d0
        eps=1d-10
        nexdp=1

        call cpu_time(t1)
!$      t1=omp_get_wtime()
        call zrootsq_adap(0,ifnewton,1,
     1      test_fun9,par1,par2,norder,eps,nexdp,center,sqw,
     2      nrtot,errest,roots,centers,nc,ier)
        call cpu_time(t2)
!$      t2=omp_get_wtime()
        call prinf('nrtot *',nrtot,1)
        call prinf('nc (leaf boxes) *',nc,1)
        call prinf('ier *',ier,1)
        call prin2('max errest *',errest,1)
        write(6,*) 'time (sec): ',t2-t1
        print *,''

c
c === Test 5: test_fun9 with nexdp=0 (no extra subdivision) ===
c should still find all roots but with less dedup overhead
c
        print *,'=== test_fun9: nexdp=0 ==='
        norder=40
        center=0d0
        sqw=2.13d0
        eps=1d-10
        nexdp=0

        call cpu_time(t1)
!$      t1=omp_get_wtime()
        call zrootsq_adap(0,ifnewton,1,
     1      test_fun9,par1,par2,norder,eps,nexdp,center,sqw,
     2      nrtot,errest,roots,centers,nc,ier)
        call cpu_time(t2)
!$      t2=omp_get_wtime()
        call prinf('nrtot *',nrtot,1)
        call prinf('nc (leaf boxes) *',nc,1)
        call prinf('ier *',ier,1)
        write(6,*) 'time (sec): ',t2-t1
        print *,''

c
c === Test 6: sin(10z)/(z-0.3-0.2i) on square center=0, sqw=2 ===
c roots of sin(10z) at k*pi/10, pole at 0.3+0.2i
c expect ~12 roots inside the square [-1,1]^2
c
        print *,'=== sin(10z)/(z-0.3-0.2i), adap ==='
        norder=40
        center=0d0
        sqw=2d0
        eps=1d-10
        nexdp=0

        call cpu_time(t1)
!$      t1=omp_get_wtime()
        call zrootsq_adap(0,ifnewton,1,
     1      test_fun_pole1,par1,par2,norder,eps,nexdp,center,sqw,
     2      nrtot,errest,roots,centers,nc,ier)
        call cpu_time(t2)
!$      t2=omp_get_wtime()
        call prinf('nrtot *',nrtot,1)
        call prinf('nc (leaf boxes) *',nc,1)
        call prinf('ier *',ier,1)
        call prin2('max errest *',errest,1)
        write(6,*) 'time (sec): ',t2-t1
        print *,''

c
c === Test 7: test_fun9/(z-0.5-0.3i) — 104 roots + pole ===
c same oscillatory function as test_fun9 but with a pole
c
        print *,'=== test_fun9/(z-0.5-0.3i), adap ==='
        norder=40
        center=0d0
        sqw=2.13d0
        eps=1d-10
        nexdp=0

        call cpu_time(t1)
!$      t1=omp_get_wtime()
        call zrootsq_adap(0,ifnewton,1,
     1      test_fun_pole2,par1,par2,norder,eps,nexdp,center,sqw,
     2      nrtot,errest,roots,centers,nc,ier)
        call cpu_time(t2)
!$      t2=omp_get_wtime()
        call prinf('nrtot *',nrtot,1)
        call prinf('nc (leaf boxes) *',nc,1)
        call prinf('ier *',ier,1)
        errmax=0
        do i=1,nrtot
          if(errest(i).gt.errmax) errmax=errest(i)
        enddo
        call prin2('max errest *',errmax,1)
        write(6,*) 'time (sec): ',t2-t1
        print *,''

c
c === Test 8: sin(z-100-100i) on square center=100+100i, sqw=20 ===
c roots at 100+100i+k*pi, large absolute values
c tests dedup with large roots and subdivision
c
        print *,'=== sin(z-100-100i), center=100+100i, sqw=1000 ==='
        norder=40
        center=100d0+100d0*ima
        sqw=1000d0
        eps=1d-10
        nexdp=0

        call cpu_time(t1)
!$      t1=omp_get_wtime()
        call zrootsq_adap(0,ifnewton,1,
     1      test_fun_bigsin,par1,par2,norder,eps,nexdp,center,sqw,
     2      nrtot,errest,roots,centers,nc,ier)
        call cpu_time(t2)
!$      t2=omp_get_wtime()
        call prinf('nrtot *',nrtot,1)
        call prinf('nc (leaf boxes) *',nc,1)
        call prinf('ier *',ier,1)
        errmax=0
        do i=1,nrtot
          if(errest(i).gt.errmax) errmax=errest(i)
        enddo
        call prin2('max errest *',errmax,1)
        call prin2_long('roots found *',roots,2*nrtot)
        write(6,*) 'time (sec): ',t2-t1
        print *,''

        stop
        end


c
c test functions with poles
c

c       sin(z-100-100i), roots at 100+100i+k*pi
        subroutine test_fun_bigsin(z,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        complex *16 z,val,dval,ima,zs
        data ima/(0d0,1d0)/
        zs=z-(100d0+100d0*ima)
        val=sin(zs)
        dval=cos(zs)
        return
        end

        subroutine test_fun_pole1(z,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        complex *16 z,val,dval,ima,f,df,g,dg
        data ima/(0d0,1d0)/
c       sin(10z)/(z-0.3-0.2i)
        f=sin(10*z)
        df=10*cos(10*z)
        g=z-0.3d0-0.2d0*ima
        dg=1
        val=f/g
        dval=(df*g-f*dg)/g**2
        return
        end

        subroutine test_fun_pole2(z,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        complex *16 z,val,dval,ima,f,df,g,dg,r0
        data ima/(0d0,1d0)/
c       test_fun9(z)/(z-0.5-0.3i)
        pi=atan(1d0)*4
        d=50
        phi=1.1
        x0=2
        r0=exp(phi*ima*z-ima*pi/4)
        f=sin(d/(z*r0-x0*1d0))
        df=-cos(d/(z*r0-x0*1d0))*(d/(z*r0-x0*1d0)**2)*
     1      (r0*ima*phi*z+r0)
        g=z-0.5d0-0.3d0*ima
        dg=1
        val=f/g
        dval=(df*g-f*dg)/g**2
        return
        end
