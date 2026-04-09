c
c     Test zroots_cheb: complex roots of complex-valued function on [a,b]
c
        implicit real *8 (a-h,o-z)
        complex *16 croots(1000),ima,ztrue
        real *8 errest(1000),par1(10),par2(1)
        data ima/(0d0,1d0)/

        external zfun_x2p1,zfun_shift,zfun_sincos
        external zfun_mixed,zfun_expim1,zfun_3near
        external zfun_dist,zfun_dist2,zfun_layers,zfun_sinz

        call prini(6,13)
        eps=epsilon(1d0)

        ifnewton = 1
        ifnewton = 0
c
c test 1: x^2+1, roots at +/-i (distance 1 from real axis)
c
        print *,'=== x^2+1 on [-2,2], delta=1.1 ==='
        call zroots_cheb(ifnewton,zfun_x2p1,par1,par2,
     1      -2d0,2d0,40,eps,1.1d0,croots,nroots,errest,ier)
        call prinf('nroots *',nroots,1)
        call prinf('ier *',ier,1)
        call prin2('croots *',croots,2*nroots)
        call prin2('errest *',errest,nroots)
        print *,''

c
c test 2: (z-0.5)+0.1i, root at 0.5-0.1i
c
        print *,'=== (z-0.5)+0.1i on [-1,1], delta=0.2 ==='
        call zroots_cheb(ifnewton,zfun_shift,par1,par2,
     1      -1d0,1d0,40,eps,0.2d0,croots,nroots,errest,ier)
        call prinf('nroots *',nroots,1)
        call prinf('ier *',ier,1)
        call prin2('croots *',croots,2*nroots)
        call prin2('errest *',errest,nroots)
        print *,''

c
c test 3: sin(z)+i*cos(z) = i*exp(-iz), no zeros
c
        print *,'=== sin+i*cos on [-5,5], delta=0.5 (expect 0) ==='
        call zroots_cheb(0,zfun_sincos,par1,par2,
     1      -5d0,5d0,60,eps,0.5d0,croots,nroots,errest,ier)
        call prinf('nroots *',nroots,1)
        call prinf('ier *',ier,1)
        print *,''

c
c test 4: (z^2+0.01)*(z-0.3+0.05i), 3 roots near real axis
c
        print *,'=== (z^2+0.01)*(z-0.3+0.05i), delta=0.2 ==='
        call zroots_cheb(ifnewton,zfun_mixed,par1,par2,
     1      -1d0,1d0,40,eps,0.2d0,croots,nroots,errest,ier)
        call prinf('nroots *',nroots,1)
        call prinf('ier *',ier,1)
        call prin2('croots *',croots,2*nroots)
        call prin2('errest *',errest,nroots)
        print *,''

c
c test 5: exp(iz)-1, roots at 2*k*pi, expect 3 on [-10,10]
c
        print *,'=== exp(iz)-1 on [-10,10], delta=0.1 ==='
        call zroots_cheb(ifnewton,zfun_expim1,par1,par2,
     1      -10d0,10d0,100,eps,0.1d0,croots,nroots,errest,ier)
        call prinf('nroots *',nroots,1)
        call prinf('ier *',ier,1)
        call prin2('croots *',croots,2*nroots)
        call prin2('errest *',errest,nroots)
        print *,''

c
c test 6: 3 roots near real axis, delta=0.1
c roots at 1-0.01i, -1+0.01i (inside), 0.5i (outside)
c
        print *,'=== 3 near-real roots, delta=0.2 ==='
        call zroots_cheb(ifnewton,zfun_3near,par1,par2,
     1      -2d0,2d0,40,eps,0.2d0,croots,nroots,errest,ier)
        call prinf('nroots *',nroots,1)
        call prinf('ier *',ier,1)
        call prin2('croots *',croots,2*nroots)
        call prin2('errest *',errest,nroots)
        print *,''

c
c test 7: convergence failure
c
        print *,'=== exp(iz)-1, n=20 (expect ier=10) ==='
        call zroots_cheb(0,zfun_expim1,par1,par2,
     1      -10d0,10d0,20,eps,0.1d0,croots,nroots,errest,ier)
        call prinf('nroots *',nroots,1)
        call prinf('ier *',ier,1)
        print *,''

c
c test 8: accuracy vs distance from real axis (linear)
c
        print *,'=== root at 0.5+d*i, varying d ==='
        do idist=1,7
          d=10d0**(-3+0.5d0*(idist-1))
          par1(1)=d
          call zroots_cheb(ifnewton,zfun_dist,par1,par2,
     1        -1d0,1d0,40,eps,d+0.01d0,croots,nroots,errest,ier)

          ztrue=0.5d0+d*ima
          if(nroots.ge.1) then
            err=abs(croots(1)-ztrue)
            write(6,'(a,e8.1,a,i3,a,e10.3,a,e10.3)')
     1          '    d=',d,' nroots=',nroots,' root err=',err,
     2          ' errest=',errest(1)
          else
            write(6,'(a,e8.1,a,i3,a,i3)')
     1          '    d=',d,' nroots=',nroots,' ier=',ier
          endif
        enddo
        print *,''

c
c test 9: accuracy vs distance (nonlinear, z^2+(0.25+d*i))
c
        print *,'=== z^2+(0.25+d*i), nonlinear, varying d ==='
        do idist=1,7
          d=10d0**(-3+0.5d0*(idist-1))
          par1(1)=d

          ztrue=sqrt(dcmplx(-0.25d0,-d))

          call zroots_cheb(ifnewton,zfun_dist2,par1,par2,
     1        -1d0,1d0,40,eps,0.55d0,croots,nroots,errest,ier)

          if(nroots.ge.1) then
            errmin=1d20
            do i=1,nroots
              err=abs(croots(i)-ztrue)
              if(err.lt.errmin) then
                errmin=err
                ibest=i
              endif
            enddo
            write(6,'(a,e8.1,a,i3,a,e10.3,a,e10.3)')
     1          '    d=',d,' nroots=',nroots,' root err=',errmin,
     2          ' errest=',errest(ibest)
          else
            write(6,'(a,e8.1,a,i3,a,i3)')
     1          '    d=',d,' nroots=',nroots,' ier=',ier
          endif
        enddo

c
c test 10: sin(z) on [100,110], large roots
c
        print *,'=== sin(z) on [100,200], delta=0.1 ==='
        call zroots_cheb(1,zfun_sinz,par1,par2,
     1      100d0,200d0,400,eps,0.1d0,croots,nroots,errest,ier)
        call prinf('nroots *',nroots,1)
        call prinf('ier *',ier,1)
        call prin2('croots *',croots,2*nroots)
        call prin2('errest *',errest,nroots)
        print *,''

c
c test 11: delta effectiveness
c f(z) = (z-0.5)*(z-0.1i)*(z-0.3i)*(z-0.7i)*(z-1.5i)
c roots at 0.5 (real), 0.1i, 0.3i, 0.7i, 1.5i
c vary delta to see which roots are found
c
        print *,'=== delta effectiveness: roots at 0.5,0.1i,'
        print *,'    0.3i,0.7i,1.5i ==='
        do idelta=1,6
          dl=0.05d0*2**(idelta-1)
          call zroots_cheb(ifnewton,zfun_layers,par1,par2,
     1        -2d0,2d0,40,eps,dl,croots,nroots,errest,ier)
          write(6,'(a,f6.3,a,i3,a)',advance='no')
     1        '    delta=',dl,' nroots=',nroots,' roots='
          do i=1,nroots
            write(6,'(a,f7.4,a,f7.4,a)',advance='no')
     1          ' (',dreal(croots(i)),',',dimag(croots(i)),')'
          enddo
          write(6,*)
        enddo

        stop
        end


c
c test functions: fun(z, par1, par2, val, dval)
c z is complex*16, val=f(z), dval=f'(z)
c

        subroutine zfun_dist(z,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        complex *16 z,val,dval,ima
        data ima/(0d0,1d0)/
c       root at z = 0.5 + par1(1)*i
        val=z-0.5d0-par1(1)*ima
        dval=1
        return
        end

        subroutine zfun_dist2(z,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        complex *16 z,val,dval,ima
        data ima/(0d0,1d0)/
c       z^2 + (0.25 + par1(1)*i)
        val=z**2+0.25d0+par1(1)*ima
        dval=2*z
        return
        end

        subroutine zfun_x2p1(z,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        complex *16 z,val,dval
        val=z**2+1
        dval=2*z
        return
        end

        subroutine zfun_shift(z,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        complex *16 z,val,dval,ima
        data ima/(0d0,1d0)/
        val=(z-0.5d0)+0.1d0*ima
        dval=1
        return
        end

        subroutine zfun_sincos(z,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        complex *16 z,val,dval,ima
        data ima/(0d0,1d0)/
        val=sin(z)+ima*cos(z)
        dval=cos(z)-ima*sin(z)
        return
        end

        subroutine zfun_mixed(z,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        complex *16 z,val,dval,ima
        data ima/(0d0,1d0)/
c       roots at z=+/-0.1i and z=0.3-0.05i
        val=(z**2+0.01d0)*(z-0.3d0+0.05d0*ima)
        dval=2*z*(z-0.3d0+0.05d0*ima)+(z**2+0.01d0)
        return
        end

        subroutine zfun_expim1(z,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        complex *16 z,val,dval,ima
        data ima/(0d0,1d0)/
c       exp(iz)-1, roots at z=2*k*pi
        val=exp(ima*z)-1
        dval=ima*exp(ima*z)
        return
        end

        subroutine zfun_3near(z,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        complex *16 z,val,dval,ima
        data ima/(0d0,1d0)/
c       roots at 1-0.01i, -1+0.01i, 0.5i
        val=(z-1d0+0.01d0*ima)*(z+1d0-0.01d0*ima)
     1      *(z-0.5d0*ima)
        dval=(z+1d0-0.01d0*ima)*(z-0.5d0*ima)
     1      +(z-1d0+0.01d0*ima)*(z-0.5d0*ima)
     2      +(z-1d0+0.01d0*ima)*(z+1d0-0.01d0*ima)
        return
        end

        subroutine zfun_sinz(z,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        complex *16 z,val,dval
        val=sin(z)
        dval=cos(z)
        return
        end

        subroutine zfun_layers(z,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        complex *16 z,val,dval,ima
        complex *16 z1,z2,z3,z4,z5
        data ima/(0d0,1d0)/
c       roots at 0.5, 0.1i, 0.3i, 0.7i, 1.5i
        z1=z-0.5d0
        z2=z-0.1d0*ima
        z3=z-0.3d0*ima
        z4=z-0.7d0*ima
        z5=z-1.5d0*ima
        val=z1*z2*z3*z4*z5
        dval=z2*z3*z4*z5+z1*z3*z4*z5+z1*z2*z4*z5
     1      +z1*z2*z3*z5+z1*z2*z3*z4
        return
        end
