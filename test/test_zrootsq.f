c
c     Test driver for zrootsq complex rootfinder.
c     Tests several functions with known roots.
c
        implicit real *8 (a-h,o-z)
        dimension errest(10 000),par1(1),par2(1)
        complex *16 center,ima
        complex *16 roots(10 000)
        data ima/(0.0d0,1.0d0)/

        external test_fun1
        external test_fun2
        external test_fun3
        external test_fun4
        external test_fun5
        external test_fun_big

        call prini(6,13)

        ifres=1
        ifnewton=1
        ifprint=1
        eps=1d-14

c
c === Test 1: polynomial with high-multiplicity roots ===
c roots: 0.5 (mult 5), 0.9 (mult 4), -0.8, 0.7i, -0.1i (mult 3)
c 14 roots total (counting multiplicity)
c
        print *,'=== test_fun1: polynomial, high multiplicity ==='
        norder=40
        center=0d0
        sqw=3d0
        call cpu_time(t1)
        call zrootsq(ifprint,ifnewton,ifres,
     1      test_fun1,par1,par2,norder,eps,center,sqw,
     2      nrtot,errest,roots,ier)
        call cpu_time(t2)
        call prin2_long('roots found *',roots,2*nrtot)
        call prin2('root error estimates *',errest,nrtot)
        call prinf('ier *',ier,1)
        call prinf('nrtot *',nrtot,1)
        write(6,*) 'time (sec): ',t2-t1
        print *,''

c
c === Test 2: polynomial with simple roots ===
c roots: 0.5, 0.9, -0.8, 0.7i, -0.1i
c 5 roots total
c
        print *,'=== test_fun2: polynomial, simple roots ==='
        norder=40
        center=0d0
        sqw=3d0
        call cpu_time(t1)
        call zrootsq(ifprint,ifnewton,ifres,
     1      test_fun2,par1,par2,norder,eps,center,sqw,
     2      nrtot,errest,roots,ier)
        call cpu_time(t2)
        call prin2_long('roots found *',roots,2*nrtot)
        call prin2('root error estimates *',errest,nrtot)
        call prinf('ier *',ier,1)
        call prinf('nrtot *',nrtot,1)
        write(6,*) 'time (sec): ',t2-t1
        print *,''

c
c === Test 3: sin(3*pi*x/2)/(x-2) ===
c roots at x = 2k/3 for integer k, pole at x=2
c search region centered at 1-0.1i, side 2.13
c expect 3 roots: 0, 2/3, 4/3
c
        print *,'=== test_fun3: sin(3*pi*x/2)/(x-2) ==='
        norder=40
        center=1d0-0.1d0*ima
        sqw=2.13d0
        call cpu_time(t1)
        call zrootsq(ifprint,ifnewton,ifres,
     1      test_fun3,par1,par2,norder,eps,center,sqw,
     2      nrtot,errest,roots,ier)
        call cpu_time(t2)
        call prin2_long('roots found *',roots,2*nrtot)
        call prin2('root error estimates *',errest,nrtot)
        call prinf('ier *',ier,1)
        call prinf('nrtot *',nrtot,1)
        write(6,*) 'time (sec): ',t2-t1
        print *,''

c
c === Test 4: sinh(3*pi*x/2)/(x-2) ===
c roots at x = 2ki/3 for integer k (imaginary axis), pole at x=2
c search small region on imaginary axis to keep sinh bounded
c expect 3 roots: -2i/3, 0, 2i/3
c
        print *,'=== test_fun4: sinh(3*pi*x/2)/(x-2) ==='
        norder=80
        center=0d0
        sqw=1.8d0
        call cpu_time(t1)
        call zrootsq(ifprint,ifnewton,ifres,
     1      test_fun4,par1,par2,norder,eps,center,sqw,
     2      nrtot,errest,roots,ier)
        call cpu_time(t2)
        call prin2_long('roots found *',roots,2*nrtot)
        call prin2('root error estimates *',errest,nrtot)
        call prinf('ier *',ier,1)
        call prinf('nrtot *',nrtot,1)
        write(6,*) 'time (sec): ',t2-t1
        print *,''

c
c === Test 5: sin(3*pi*x)/(x-2) ===
c roots at x = k/3 for integer k, pole at x=2
c search region centered at 0, side 1.5
c expect 5 roots: -2/3, -1/3, 0, 1/3, 2/3
c
        print *,'=== test_fun5: sin(3*pi*x)/(x-2) ==='
        norder=60
        center=0d0
        sqw=1.5d0
        call cpu_time(t1)
        call zrootsq(ifprint,ifnewton,ifres,
     1      test_fun5,par1,par2,norder,eps,center,sqw,
     2      nrtot,errest,roots,ier)
        call cpu_time(t2)
        call prin2_long('roots found *',roots,2*nrtot)
        call prin2('root error estimates *',errest,nrtot)
        call prinf('ier *',ier,1)
        call prinf('nrtot *',nrtot,1)
        write(6,*) 'time (sec): ',t2-t1
        print *,''

c
c === Test 6: polynomial with large roots ===
c roots at 100+100i, 102+99i, 98+101i
c tests dedup with large absolute values
c
        print *,'=== test_fun_big: roots near 100+100i ==='
        norder=40
        center=100d0+100d0*ima
        sqw=10d0
        call cpu_time(t1)
        call zrootsq(ifprint,ifnewton,ifres,
     1      test_fun_big,par1,par2,norder,eps,center,sqw,
     2      nrtot,errest,roots,ier)
        call cpu_time(t2)
        call prin2_long('roots found *',roots,2*nrtot)
        call prin2('root error estimates *',errest,nrtot)
        call prinf('ier *',ier,1)
        call prinf('nrtot *',nrtot,1)
        write(6,*) 'time (sec): ',t2-t1
        print *,''

        stop
        end
