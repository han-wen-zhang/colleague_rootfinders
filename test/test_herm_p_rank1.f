        implicit real *8 (a-h,o-z)
        dimension w(1 000 000)
        real *8 ca(10 000 000),cb(10 000 000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        n=100
        call prinf('n=*',n,1)
c 
c
        call testit(n,ca,cb)

        stop
        end
c
c
c
c 
c
        subroutine testit(n,a,b)
        implicit real *8 (a-h,o-z)
        complex *16 p(10 000),q(10 000),a(n,n),b(n,n),
     1      w(1 000 000),w2(1 000 000),clams(10 000),
     2      supdiag(10 000),diag(10 000),ima
        data ima/(0.0d0,1.0d0)/
c
c        construct the vectors p, q and the symmetric matrix
c
        call set_zero(p,n*2)
        call set_zero(q,n*2)
        call set_zero(a,n*n*2)
        call set_zero(b,n*n*2)
c
        call corrand_norm(n*10,w,w2)

        
        do 1200 i=1,n-1
c
        a(i,i+1)=w(i)
        a(i+1,i)=conjg(w(i))

        supdiag(i)=a(i,i+1)
 1200 continue

        do 1300 i=1,n
c
        a(i,i)=sqrt(i*1.1d0) +0.3d0*ima
        diag(i)=a(i,i)
c
        q(i)=i
        q(i)=i+ima/i
 1300 continue

        p(n)=1
        p(n-1)=1
        p(n-1)=0.321d0
        p(n-1)=0.321d0 +0.2*ima

c
c        form the matrix b = a + p \circ q*
c
        do 1600 i=1,n
        do 1400 j=1,n
c
        b(i,j)=a(i,j)+p(i)*conjg(q(j))
 1400 continue
 1600 continue

c
c        compute the eigenvalues naively
c
        call zgeigvals(n,b,clams,ier)

        call bubble(clams,n)

        call prin2('and clams=*',clams,n*2)
c
c        and using the specialized eigensolver
c
        eps=epsilon(1d0)
        call herm_p_rank1(ier,diag,supdiag,p,q,n,eps,
     1      w,niter_max,aver_niter)

        call prinf('after _rank1, ier=',ier,1)
        call prinf('and niter_max=',niter_max,1)
        call prin2('and aver_niter=',aver_niter,1)
c
        call bubble(diag,n)
        call prin2('and sorted, diag=*',diag,n*2)
c
c        ... compare
c
        call arr_subtr(clams,diag,w,n*2)

        call scapro(w,w,n*2,d)

        d=sqrt(d)

        call prin2('and norm of clams4-diag=*',d,1)

        return
        end
c       
c 
c 
c 
c 
        subroutine set_zero(x,n)
        implicit real *8 (a-h,o-z)
        save
        dimension x(n)
c 
        do 2000 i=1,n
        x(i)=0
 2000 continue
        return
        end
c
c
c
c
c
        subroutine arr_subtr(a,b,c,n)
        implicit real *8 (a-h,o-z)
        dimension b(*),a(*),c(*)

        do 1400 i=1,n
c
        c(i)=a(i)-b(i)
 1400 continue
c
        return
        end
c
c
c 
c 
c  
        subroutine scapro(x,y,n,prod)
        implicit real *8 (a-h,o-z)
        save
        real *8 x(*),y(*),prod
c
c       evaluate the polynomial at the point z
c
        prod=0
        do 1200 i=1,n
c
        prod=prod+x(i)*y(i)
 1200 continue
c
        return
        end
c
c
c
        subroutine zgeigvals(n,a,w,ier)
        implicit none
        integer n,ier,info,lwork,j,k
        complex*16 a(n,n),w(n)
        complex*16, allocatable :: acopy(:,:),work(:)
        complex*16 vl(1,1),vr(1,1)
        real*8, allocatable :: rwork(:)
        complex*16 wopt
c
        allocate(acopy(n,n),rwork(2*n))
        do k = 1,n
          do j = 1,n
            acopy(j,k) = a(j,k)
          enddo
        enddo
c
        call zgeev('N','N',n,acopy,n,w,vl,1,vr,1,
     1    wopt,-1,rwork,info)
        lwork = int(dble(wopt))
        allocate(work(lwork))
        call zgeev('N','N',n,acopy,n,w,vl,1,vr,1,
     1    work,lwork,rwork,info)
        ier = info
        deallocate(acopy,work,rwork)
        return
        end
c
c
c
        subroutine bubble(cs,n)
        implicit real *8 (a-h,o-z)
        complex *16 cs(n),cc,ima
        data ima/(0.0d0,1.0d0)/
c
        do 1400 i=1,n
c
        do 1200 j=1,n-1
c
        d=cs(j)
        dd=cs(j+1)

        ddd=abs(d-dd)
        if(ddd .lt. 1.0d-10) d=d+ima*cs(j)
        if(abs(ddd) .lt. 1.0d-10) dd=dd+ima*cs(j+1)
c
        if(dd .le. d) goto 1200
        cc=cs(j)
        cs(j)=cs(j+1)
        cs(j+1)=cc
 1200 continue
 1400 continue
c
        return
        end
c
