c
c     Test zroots_cheb_adap: adaptive complex rootfinder on [a,b]
c
        implicit real *8 (a-h,o-z)
        complex *16 croots(10 000),ima
        real *8 errest(10 000),par1(1),par2(1)
        data ima/(0d0,1d0)/

        external zfun_expim1,zfun_3poles,zfun_mixed
        external zfun_sin10pole,zfun_sin100_3p,zfun_expipole
        external zfun_shifted3p,zfun_shifted,zfun_shifted2
        external zfun_sin100_3poff
        external zfun_rootpole

        call prini(6,13)
        eps=epsilon(1d0)

c
c test 1: exp(iz)-1 on [-10,10], n=20 (must subdivide)
c roots at 2*k*pi, expect 3 roots: 0, +-2pi
c
        print *,'=== exp(iz)-1 on [-10,10], n=20, delta=0.1 ==='
        ifn=1
        call zroots_cheb_adap(ifn,zfun_expim1,par1,par2,
     1      -10d0,10d0,20,eps,0.1d0,0,
     2      croots,nroots,errest,ier)
        call test_report(ifn,nroots,ier,errest)

c
c test 2: exp(iz)-1 on [-10,10], n=100 (no subdivision)
c expect 3 roots
c
        print *,'=== exp(iz)-1 on [-10,10], n=100, delta=0.1 ==='
        ifn=1
        call zroots_cheb_adap(ifn,zfun_expim1,par1,par2,
     1      -10d0,10d0,100,eps,0.1d0,0,
     2      croots,nroots,errest,ier)
        call test_report(ifn,nroots,ier,errest)

c
c test 3: (z^2+0.01)*(z-0.3+0.05i) on [-5,5], n=20
c expect 3 roots: +-0.1i, 0.3-0.05i
c
        print *,'=== (z^2+0.01)*(z-0.3+0.05i) on [-5,5], n=20 ==='
        ifn=1
        call zroots_cheb_adap(ifn,zfun_mixed,par1,par2,
     1      -5d0,5d0,20,eps,0.2d0,0,
     2      croots,nroots,errest,ier)
        call test_report(ifn,nroots,ier,errest)
        call prin2('croots *',croots,2*nroots)

c
c test 4: (exp(iz)-1) / (z-0.5) on [-10,10], n=20
c expect 3 roots: 0, +-2pi (pole at 0.5 doesn't cancel a root)
c
        print *,'=== exp(iz)-1/(z-0.5) on [-10,10], n=20 ==='
        ifn=1
        call zroots_cheb_adap(ifn,zfun_3poles,par1,par2,
     1      -10d0,10d0,20,eps,0.1d0,0,
     2      croots,nroots,errest,ier)
        call test_report(ifn,nroots,ier,errest)

c
c test 5: same as test 1 but newton off, expect 3 roots
c
        print *,'=== exp(iz)-1 on [-10,10], n=20, newton=0 ==='
        ifn=0
        call zroots_cheb_adap(ifn,zfun_expim1,par1,par2,
     1      -10d0,10d0,20,eps,0.1d0,0,
     2      croots,nroots,errest,ier)
        call test_report(ifn,nroots,ier,errest)

c
c test 6: sin(10z)/(z-0.5) on [-5,5], many roots + pole
c sin(10z)=0 at z=k*pi/10, expect 31 roots (k=-15..15)
c
        print *,'=== sin(10z)/(z-0.5) on [-5,5], n=20, delta=0.2 ==='
        ifn=1
        call zroots_cheb_adap(ifn,zfun_sin10pole,par1,par2,
     1      -5d0,5d0,20,eps,0.2d0,0,
     2      croots,nroots,errest,ier)
        call test_report(ifn,nroots,ier,errest)

c
c test 7: sin(100z)/((z-0.3)(z+1.5)(z-3.7)) on [-5,5]
c expect 319 roots (k=-159..159, k*pi/100), 3 poles on real axis
c
        print *,'=== sin(100z)/3poles on [-5,5], n=20, delta=0.2 ==='
        ifn=1
        call zroots_cheb_adap(ifn,zfun_sin100_3p,par1,par2,
     1      -5d0,5d0,20,eps,0.2d0,0,
     2      croots,nroots,errest,ier)
        call test_report(ifn,nroots,ier,errest)
        nreal=0
        do i=1,nroots
          if(abs(dimag(croots(i))).lt.1d-10) nreal=nreal+1
        enddo
        write(6,*) 'real roots:',nreal,' complex:',nroots-nreal
c       check for near-duplicates
        ndup=0
        dmax=0
        do i=1,nroots
          do j=i+1,nroots
            d=abs(croots(i)-croots(j))
            if(d.lt.1d-6) then
              ndup=ndup+1
              if(d.gt.dmax) dmax=d
            endif
          enddo
        enddo
        write(6,*) 'near-dups:',ndup,' max sep:',dmax

c
c test 8: exp(iz)/(z-1) on [-10,10], expect 3 roots + pole at 1
c
        print *,'=== exp(iz)/(z-1) on [-10,10], n=20, delta=0.1 ==='
        ifn=1
        call zroots_cheb_adap(ifn,zfun_expipole,par1,par2,
     1      -10d0,10d0,20,eps,0.1d0,0,
     2      croots,nroots,errest,ier)
        call test_report(ifn,nroots,ier,errest)
        call prin2('croots *',croots,2*nroots)

c
c test 9: sin(100(z-0.003i)), no poles
c expect 319 roots at k*pi/100 + 0.003i (shifted off real axis)
c needs large delta=1.0 so subinterval ellipses contain roots
c
        print *,'=== sin(100(z-0.003i)), delta=1.0, newton OFF ==='
        ifn=0
        call zroots_cheb_adap(ifn,zfun_shifted2,par1,par2,
     1      -5d0,5d0,40,eps,1.0d0,0,
     2      croots,nroots,errest,ier)
        call test_report(ifn,nroots,ier,errest)

        print *,'=== sin(100(z-0.003i)), delta=1.0, newton ON ==='
        ifn=1
        call zroots_cheb_adap(ifn,zfun_shifted2,par1,par2,
     1      -5d0,5d0,40,eps,1.0d0,0,
     2      croots,nroots,errest,ier)
        call test_report(ifn,nroots,ier,errest)
c       also try without the |fun| filter by using the
c       non-adaptive version on a single interval
        print *,'  (non-adaptive, n=400 for reference:)'
        call zroots_cheb(1,zfun_shifted3p,par1,par2,
     1      -5d0,5d0,400,eps,0.1d0,croots,nroots,errest,ier)
        call prinf('  nroots *',nroots,1)
        nreal=0
        do i=1,nroots
          if(abs(dimag(croots(i))).lt.1d-10) nreal=nreal+1
        enddo
        write(6,*) 'real roots:',nreal,' complex:',nroots-nreal

c
c test 10: sin(100z)/((z-0.3+0.02i)(z+1.5+0.03i)(z-3.7+0.01i))
c expect 319 roots, 3 poles slightly below real axis
c poles don't hit Chebyshev nodes, expansion converges everywhere
c
c test 11: sin(100z)/(z-pi/100) on [-5,5]
c root at z=pi/100 coincides with pole (removable singularity)
c expect 318 roots (319 minus the one at the pole)
c
        print *,'=== sin(100z)/3poles below real, n=20, delta=0.1 ==='
        ifn=1
        call zroots_cheb_adap(ifn,zfun_sin100_3poff,par1,par2,
     1      -5d0,5d0,20,eps,0.1d0,0,
     2      croots,nroots,errest,ier)
        call test_report(ifn,nroots,ier,errest)

        print *,'=== sin(100z)/(z-pi/100) on [-5,5], n=20 ==='
        ifn=1
        call zroots_cheb_adap(ifn,zfun_rootpole,par1,par2,
     1      -5d0,5d0,20,eps,0.1d0,0,
     2      croots,nroots,errest,ier)
        call test_report(ifn,nroots,ier,errest)

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


c
c test functions: fun(z,par1,par2,val,dval)
c

        subroutine zfun_expim1(z,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        complex *16 z,val,dval,ima
        data ima/(0d0,1d0)/
        val=exp(ima*z)-1
        dval=ima*exp(ima*z)
        return
        end

        subroutine zfun_mixed(z,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        complex *16 z,val,dval,ima
        data ima/(0d0,1d0)/
        val=(z**2+0.01d0)*(z-0.3d0+0.05d0*ima)
        dval=2*z*(z-0.3d0+0.05d0*ima)+(z**2+0.01d0)
        return
        end

        subroutine zfun_3poles(z,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        complex *16 z,val,dval,ima,f,df,g,dg
        data ima/(0d0,1d0)/
c       (exp(iz)-1)/(z-0.5), pole at z=0.5
        f=exp(ima*z)-1
        df=ima*exp(ima*z)
        g=z-0.5d0
        dg=1
        val=f/g
        dval=(df*g-f*dg)/g**2
        return
        end

c       sin(10z)/(z-0.5), roots at k*pi/10 + pole at 0.5
        subroutine zfun_sin10pole(z,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        complex *16 z,val,dval,f,df,g,dg
        f=sin(10*z)
        df=10*cos(10*z)
        g=z-0.5d0
        dg=1
        val=f/g
        dval=(df*g-f*dg)/g**2
        return
        end

c       sin(100z)/((z-0.3)(z+1.5)(z-3.7)), 3 poles
        subroutine zfun_sin100_3p(z,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        complex *16 z,val,dval,f,df,g,dg
        complex *16 p1,p2,p3
        f=sin(100*z)
        df=100*cos(100*z)
        p1=z-0.3d0
        p2=z+1.5d0
        p3=z-3.7d0
        g=p1*p2*p3
        dg=p2*p3+p1*p3+p1*p2
        val=f/g
        dval=(df*g-f*dg)/g**2
        return
        end

c       (exp(iz)-1)/(z-1), pole at z=1
        subroutine zfun_expipole(z,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        complex *16 z,val,dval,ima,f,df,g,dg
        data ima/(0d0,1d0)/
        f=exp(ima*z)-1
        df=ima*exp(ima*z)
        g=z-1d0
        dg=1
        val=f/g
        dval=(df*g-f*dg)/g**2
        return
        end

c       sin(100(z-0.03i)), no poles
        subroutine zfun_shifted(z,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        complex *16 z,val,dval,ima,zs
        data ima/(0d0,1d0)/
        zs=z-0.03d0*ima
        val=sin(100*zs)
        dval=100*cos(100*zs)
        return
        end

c       sin(100(z-0.003i)), small shift
        subroutine zfun_shifted2(z,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        complex *16 z,val,dval,ima,zs
        data ima/(0d0,1d0)/
        zs=z-0.003d0*ima
        val=sin(100*zs)
        dval=100*cos(100*zs)
        return
        end

c       sin(100(z-0.03i))/((z-0.3)(z+1.5)(z-3.7))
c       roots shifted 0.03i off real axis
        subroutine zfun_shifted3p(z,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        complex *16 z,val,dval,f,df,g,dg,ima
        complex *16 p1,p2,p3,zs
        data ima/(0d0,1d0)/
        zs=z-0.03d0*ima
        f=sin(100*zs)
        df=100*cos(100*zs)
        p1=z-0.3d0
        p2=z+1.5d0
        p3=z-3.7d0
        g=p1*p2*p3
        dg=p2*p3+p1*p3+p1*p2
        val=f/g
        dval=(df*g-f*dg)/g**2
        return
        end

c       sin(100z)/3poles slightly below real axis
        subroutine zfun_sin100_3poff(z,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        complex *16 z,val,dval,f,df,g,dg,ima
        complex *16 p1,p2,p3
        data ima/(0d0,1d0)/
        f=sin(100*z)
        df=100*cos(100*z)
        p1=z-0.3d0+0.02d0*ima
        p2=z+1.5d0+0.03d0*ima
        p3=z-3.7d0+0.01d0*ima
        g=p1*p2*p3
        dg=p2*p3+p1*p3+p1*p2
        val=f/g
        dval=(df*g-f*dg)/g**2
        return
        end

c       sin(100z)/(z-pi/100), pole coincides with root
        subroutine zfun_rootpole(z,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        complex *16 z,val,dval,f,df,g,dg
        pi=4*atan(1d0)
        f=sin(100*z)
        df=100*cos(100*z)
        g=z-pi/100
        dg=1
        val=f/g
        dval=(df*g-f*dg)/g**2
        return
        end


