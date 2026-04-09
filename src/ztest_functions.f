c
c     Test functions for the complex root finder.
c     Each subroutine has the calling convention
c          test_funN(x, par1, par2, val, dval)
c     where val = f(x), dval = f'(x), x is complex*16.
c


        subroutine test_fun1(x,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        dimension par1(*),par2(*)
        complex *16 val,dval,x
        complex *16 ima
        save
        data ima/(0.0d0,1.0d0)/

        done=1

        val=(x-0.5)**5*(x-0.9)**4*(x+0.8)*(x-0.7*ima)*(x+0.1*ima)**3

        dval=5*(x-0.5)**4*(x-0.9)**4*(x+0.8)*(x-0.7*ima)*(x+0.1*ima)**3
        dval=dval+4*(x-0.5)**5*(x+0.8)*(x-0.9)
     1      **3*(x-0.7*ima)*(x+0.1*ima)**3
        dval=dval+ (x-0.5)**5*(x-0.9)**4*(x-0.7*ima)*(x+0.1*ima)**3
        dval=dval+ (x-0.5)**5*(x-0.9)**4*(x+0.8)*(x+0.1*ima)**3
        dval=dval+3*(x-0.5)**5*(x-0.9)**4*(x+0.8)
     1      *(x-0.7*ima)*(x+0.1*ima)**2


        return
        end

        subroutine test_fun2(x,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        dimension par1(*),par2(*)
        complex *16 val,dval,x
        complex *16 ima
        save
        data ima/(0.0d0,1.0d0)/

        done=1
        val=(x-0.5)*(x-0.9)*(x+0.8)*(x-0.7*ima)*(x+0.1*ima)


        dval=(x-0.9)*(x+0.8)*(x-0.7*ima)*(x+0.1*ima)
        dval=dval+(x-0.5)*(x+0.8)*(x-0.7*ima)*(x+0.1*ima)
        dval=dval+ (x-0.5)*(x-0.9)*(x-0.7*ima)*(x+0.1*ima)
        dval=dval+ (x-0.5)*(x-0.9)*(x+0.8)*(x+0.1*ima)
        dval=dval+ (x-0.5)*(x-0.9)*(x+0.8)*(x-0.7*ima)


        return
        end


        subroutine test_fun3(x,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        dimension par1(*),par2(*)
        complex *16 val,dval,x
        complex *16 ima
        save
        data ima/(0.0d0,1.0d0)/

        done=1
        const=1
        pi=atan(done)*4
        val=sin(x*const*3*pi/2)/(x-done*2)

        dval=const*3*pi/2*cos(x*const*3*pi/2)/(x-done*2)
        dval=dval-sin(x*const*3*pi/2)/(x-done*2)**2


        return
        end

        subroutine test_fun4(x,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        dimension par1(*),par2(*)
        complex *16 val,dval,x
        complex *16 ima
        save
        data ima/(0.0d0,1.0d0)/

        done=1
        pi=atan(done)*4
        val=sinh(x*3*pi/2)/(x-done*2)

        dval=3*pi/2*cosh(x*3*pi/2)/(x-done*2)
        dval=dval-sinh(x*3*pi/2)/(x-done*2)**2


        return
        end


        subroutine test_fun5(x,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        dimension par1(*),par2(*)
        complex *16 val,dval,x
        complex *16 ima
        save
        data ima/(0.0d0,1.0d0)/

        done=1
        pi=atan(done)*4
        val=sin(x*3*pi)/(x-done*2)

        dval=3*pi*cos(x*3*pi)/(x-done*2)
        dval=dval-sin(x*3*pi)/(x-done*2)**2


        return
        end

        subroutine test_fun6(x,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        dimension par1(*),par2(*)
        complex *16 val,dval,x,r0
        complex *16 ima
        save
        data ima/(0.0d0,1.0d0)/

        done=1
        d=100
        x0=2
        r0=1
        pi=atan(done)*4
        val=sin(d/(x*r0-x0*done))

        dval=-cos(d/(x*r0-x0*done))*(d*r0/(x*r0-x0*done)**2)


        return
        end


        subroutine test_fun7(x,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        dimension par1(*),par2(*)
        complex *16 val,dval,x,r0
        complex *16 ima
        save
        data ima/(0.0d0,1.0d0)/

        done=1
        d=100
        x0=2
        r0=exp(0.4*ima)
        pi=atan(done)*4
        val=sin(d/(x*r0-x0*done))

        dval=-cos(d/(x*r0-x0*done))*(d*r0/(x*r0-x0*done)**2)


        return
        end

        subroutine test_fun8(x,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        dimension par1(*),par2(*)
        complex *16 val,dval,x,r0
        complex *16 ima
        save
        data ima/(0.0d0,1.0d0)/

        done=1
        pi=atan(done)*4
        d=100
        x0=2
        r0=exp(ima*pi/4)
        val=sin(d/(x*r0-x0*done))

        dval=-cos(d/(x*r0-x0*done))*(d*r0/(x*r0-x0*done)**2)


        return
        end


        subroutine test_fun9(x,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        dimension par1(*),par2(*)
        complex *16 val,dval,x,r0
        complex *16 ima
        save
        data ima/(0.0d0,1.0d0)/

        done=1
        pi=atan(done)*4
        d=50
        phi=1.1
        x0=2
        r0=exp(phi*ima*x-ima*pi/4)
        val=sin(d/(x*r0-x0*done))

        dval=-cos(d/(x*r0-x0*done))*(d/(x*r0-x0*done)**2)*
     1      (r0*ima*phi*x+r0)


        return
        end
