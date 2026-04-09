c
c     Test droots_cheb: real roots on [a,b]
c
        implicit real *8 (a-h,o-z)

        call prini(6,13)

c
c test 1: funuser8 on [-1,1], real roots
c
        print *,'=== funuser8 on [-1,1], real ==='
        call test_real_ab(-1d0,1d0,100)

c
c test 2: sin(x) on [0, 10], real roots
c
        print *,'=== sin(x) on [0,10], real ==='
        call test_real_sin(0d0,10d0,100)

c
c test 3: sin(x) on [-20, 20], real roots
c
        print *,'=== sin(x) on [-20,20], real ==='
        call test_real_sin(-20d0,20d0,200)

c
c test 4: sin(100x) on [-1,1] with n=20 (should fail convergence)
c
        print *,'=== sin(100x) on [-1,1], n=20 (expect ier=10) ==='
        call test_noconv(-1d0,1d0,20)

c
c test 5: sin(100x) on [-1,1] with n=200 (should converge)
c
        print *,'=== sin(100x) on [-1,1], n=200 ==='
        call test_noconv(-1d0,1d0,200)

        stop
        end


        subroutine test_real_ab(a,b,n)
        implicit real *8 (a-h,o-z)
        real *8 roots(10 000),errest(10 000),par1(1),par2(1)
        external funuser8

        eps=epsilon(1d0)

        call droots_cheb(1,funuser8,par1,par2,a,b,n,eps,
     1      roots,nroots,errest,ier)

        call prinf('nroots *',nroots,1)
        call prinf('ier *',ier,1)
        call prin2('roots *',roots,nroots)
        call prin2('errest *',errest,nroots)
        print *,''
        return
        end


        subroutine test_real_sin(a,b,n)
        implicit real *8 (a-h,o-z)
        real *8 roots(10 000),errest(10 000),par1(1),par2(1)
        external dtest_sin

        eps=epsilon(1d0)

        call droots_cheb(1,dtest_sin,par1,par2,a,b,n,eps,
     1      roots,nroots,errest,ier)

        call prinf('nroots *',nroots,1)
        call prinf('ier *',ier,1)
        call prin2('roots *',roots,nroots)
        call prin2('errest *',errest,nroots)
        print *,''
        return
        end


        subroutine test_noconv(a,b,n)
        implicit real *8 (a-h,o-z)
        real *8 roots(10 000),errest(10 000),par1(1),par2(1)
        external dtest_sin100

        eps=epsilon(1d0)

        call droots_cheb(1,dtest_sin100,par1,par2,a,b,n,eps,
     1      roots,nroots,errest,ier)

        call prinf('nroots *',nroots,1)
        call prinf('ier *',ier,1)
        if(nroots.gt.0) call prin2('roots *',roots,nroots)
        print *,''
        return
        end


        subroutine dtest_sin100(x,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        val=sin(100*x)
        dval=100*cos(100*x)
        return
        end


        subroutine dtest_sin(x,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*)
        val=sin(x)
        dval=cos(x)
        return
        end

        subroutine funuser8(x,par1,par2,val,dval)
        implicit real *8 (a-h,o-z)
        real *8 par1(*),par2(*),roots(100)

        one=1
        roots(1)=-one/2+1.0d-6
        roots(2)=-one/3+2*1.0d-6
        roots(3)=-0.61d0
        roots(4)=0.121d0
        roots(5)=one-1.0d-3
        roots(6)=one-1.0d-3
        roots(7)=one-1.0d-3
        roots(8)=one-1.0d-3
        roots(9)=one-1.0d-3
        roots(10)=one-1.0d-3
        roots(11)=one-1.0d-3
        roots(12)=one-1.0d-3
        roots(13)=one-1.0d-3
        nrts=13
        val=1
        dval=0
        do i=1,nrts
          dval=dval*( x-roots(i))+val
          val=val*(x-roots(i))
        enddo
        return
        end

