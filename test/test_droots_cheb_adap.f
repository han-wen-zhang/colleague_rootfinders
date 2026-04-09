c
c     Test droots_cheb_adap: adaptive real rootfinder
c
        implicit real *8 (a-h,o-z)
        real *8 roots(10 000),errest(10 000),par1(1),par2(1)
        external dtest_sin,dtest_sin100
        external dtest_sinpole,dtest_sinpole2,dtest_sinpole3
        external dtest_manysingular,dtest_3poles

        call prini(6,13)
        eps=epsilon(1d0)

c
c test 1: sin(100x) on [-1,1], n=20, ifnewton=1
c
        print *,'=== sin(100x) on [-1,1], n=20, newton=1 ==='
        ifn=1
        call droots_cheb_adap(ifn,dtest_sin100,par1,par2,
     1      -1d0,1d0,20,eps,0,roots,nroots,errest,ier)
        call test_report(ifn,nroots,ier,errest)

c
c test 2: sin(100x) on [-1,1], n=20, ifnewton=0
c
        print *,'=== sin(100x) on [-1,1], n=20, newton=0 ==='
        ifn=0
        call droots_cheb_adap(ifn,dtest_sin100,par1,par2,
     1      -1d0,1d0,20,eps,0,roots,nroots,errest,ier)
        call test_report(ifn,nroots,ier,errest)

c
c test 3: sin(x) on [-100,100], n=40, ifnewton=1
c
        print *,'=== sin(x) on [-100,100], n=40, newton=1 ==='
        ifn=1
        call droots_cheb_adap(ifn,dtest_sin,par1,par2,
     1      -100d0,100d0,40,eps,0,roots,nroots,errest,ier)
        call test_report(ifn,nroots,ier,errest)

c
c test 4: sin(x) on [-100,100], n=40, ifnewton=0
c
        print *,'=== sin(x) on [-100,100], n=40, newton=0 ==='
        ifn=0
        call droots_cheb_adap(ifn,dtest_sin,par1,par2,
     1      -100d0,100d0,40,eps,0,roots,nroots,errest,ier)
        call test_report(ifn,nroots,ier,errest)

c
c test 5: sin(100x) on [-1,1], n=200 (no subdivision)
c
        print *,'=== sin(100x) on [-1,1], n=200 (no subdiv) ==='
        ifn=1
        call droots_cheb_adap(ifn,dtest_sin100,par1,par2,
     1      -1d0,1d0,200,eps,0,roots,nroots,errest,ier)
        call test_report(ifn,nroots,ier,errest)

c
c test 6: sin(100x) on [-1,1], n=20, nexdp=1
c
        print *,'=== sin(100x) on [-1,1], n=20, nexdp=1 ==='
        ifn=1
        call droots_cheb_adap(ifn,dtest_sin100,par1,par2,
     1      -1d0,1d0,20,eps,1,roots,nroots,errest,ier)
        call test_report(ifn,nroots,ier,errest)

c
c test 7: sin(x)/(x-1.01) on [0,2], singularity just outside
c
        print *,'=== sin(x)/(x-1.01) on [0,2], n=40 ==='
        ifn=1
        call droots_cheb_adap(ifn,dtest_sinpole,par1,par2,
     1      0d0,2d0,40,eps,0,roots,nroots,errest,ier)
        call test_report(ifn,nroots,ier,errest)
        call prin2('roots *',roots,nroots)

c
c test 8: sin(x)/(x-1.001) on [0,2], singularity very close
c
        print *,'=== sin(x)/(x-1.001) on [0,2], n=40 ==='
        ifn=1
        call droots_cheb_adap(ifn,dtest_sinpole2,par1,par2,
     1      0d0,2d0,40,eps,0,roots,nroots,errest,ier)
        call test_report(ifn,nroots,ier,errest)
        call prin2('roots *',roots,nroots)

c
c test 9: sin(x)/(x-2) on [0,2], singularity on boundary
c
        print *,'=== sin(x)/(x-2) on [0,2], n=40 ==='
        ifn=1
        call droots_cheb_adap(ifn,dtest_sinpole3,par1,par2,
     1      0d0,2d0,40,eps,0,roots,nroots,errest,ier)
        call test_report(ifn,nroots,ier,errest)
        call prin2('roots *',roots,nroots)

c
c test 10: sin(100x)/(x-0.5) on [-1,1], many roots + pole inside
c expect ~63 roots of sin(100x), minus the one near x=0.5
c
        print *,'=== sin(100x)/(x-0.5) on [-1,1], n=20 ==='
        ifn=1
        call droots_cheb_adap(ifn,dtest_manysingular,par1,par2,
     1      -1d0,1d0,20,eps,0,roots,nroots,errest,ier)
        call test_report(ifn,nroots,ier,errest)
        call prin2('roots *',roots,nroots)

c
c test 11: sin(100x)/((x-0.3)(x+0.7)(x-0.9)) on [-1,1]
c 3 poles inside, many roots, n=20
c
        print *,'=== sin(100x)/((x-0.3)(x+0.7)(x-0.9)), n=20 ==='
        ifn=1
        call droots_cheb_adap(ifn,dtest_3poles,par1,par2,
     1      -1d0,1d0,20,eps,0,roots,nroots,errest,ier)
        call test_report(ifn,nroots,ier,errest)
        call prinf('nroots *',nroots,1)

        stop
        end


        subroutine test_report(ifnewton,nroots,ier,errest)
        implicit real *8 (a-h,o-z)
        real *8 errest(*)

        call prinf('nroots *',nroots,1)
        call prinf('ier *',ier,1)
        if(ifnewton.eq.1) then
          write(6,*) 'Newton: on'
        else
          write(6,*) 'Newton: off'
        endif
        if(nroots.gt.0) then
          errmax=0
          do i=1,nroots
            if(errest(i).gt.errmax) errmax=errest(i)
          enddo
          call prin2('max errest *',errmax,1)
        endif
        print *,''
        print *,''
        return
        end


        subroutine dtest_sin(x,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        val=sin(x)
        dval=cos(x)
        return
        end

        subroutine dtest_sin100(x,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        val=sin(100*x)
        dval=100*cos(100*x)
        return
        end

c       sin(x)/(x-1.01), pole at x=1.01 (just outside [0,2])
        subroutine dtest_sinpole(x,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        val=sin(x)/(x-1.01d0)
        dval=(cos(x)*(x-1.01d0)-sin(x))/(x-1.01d0)**2
        return
        end

c       sin(x)/(x-1.001), pole at x=1.001 (very close)
        subroutine dtest_sinpole2(x,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        val=sin(x)/(x-1.001d0)
        dval=(cos(x)*(x-1.001d0)-sin(x))/(x-1.001d0)**2
        return
        end

c       sin(x)/(x-2), pole at x=2 (on boundary)
        subroutine dtest_sinpole3(x,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        val=sin(x)/(x-2d0)
        dval=(cos(x)*(x-2d0)-sin(x))/(x-2d0)**2
        return
        end

c       sin(100x)/(x-0.5), many roots + pole at x=0.5
        subroutine dtest_manysingular(x,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        val=sin(100*x)/(x-0.5d0)
        dval=(100*cos(100*x)*(x-0.5d0)-sin(100*x))
     1      /(x-0.5d0)**2
        return
        end

c       sin(100x)/((x-0.3)(x+0.7)(x-0.9)), 3 poles
        subroutine dtest_3poles(x,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        real *8 p1,p2,p3
        p1=x-0.3d0
        p2=x+0.7d0
        p3=x-0.9d0
        val=sin(100*x)/(p1*p2*p3)
        dval=(100*cos(100*x)*p1*p2*p3
     1      -sin(100*x)*(p2*p3+p1*p3+p1*p2))
     2      /(p1*p2*p3)**2
        return
        end

