c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging code, and the beginning of the
c       diagonalization code proper. This file contains a single
c       user-callable subroutine herm_p_rank1.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
        subroutine herm_p_rank1(ier,diag,supdiag,p,q,n,eps,
     1      w,niter_max,aver_niter)
        implicit real *8 (a-h,o-z)
        complex *16 diag(*),supdiag(*),p(*),q(*),w(*)
c
c       This subroutine finds the spectrum of the lower Hessenberg
c       matrix B of the form 
c
c                    B = A + P \circ Q^*,                               (1)
c
c       where A is a Hermitian matrix, and P, Q are a pair of vectors.
c       The claim to fame of this subroutine is that it costs O(n^2)
c       operations and has a special form of backward stability, where
c       the computed eigenvalues are the exact eigenvalues of the matrix 
c
c                   (A+dA) + (P+dP) \circ (Q+dQ)^*                      (2)
c
c       with dA \approx eps||A||, dP \approx eps||P||, and 
c       dQ \approx eps||Q||.
c
c       PLEASE NOTE that this is a very specialized subroutine, in 
c       in that the matrix B is BOTH of the form (1) and lower 
c       Hessenberg, so that 
c
c                    B(i,j)=0                                           (3)
c
c       for all 
c
c                    j > i+1.                                           (4)
c
c       Note also that, for matrices of this form, A is determined
c       entirely by its diagonal and superdiagonal together with the
c       vectors P and Q.
c
c                    Input parameters:
c
c  diag - the diagonal of the matrix A in (1);
c       destroyed by this subroutine utterly
c  supdiag - the superdiagonal of the matrix A in (1);
c       destroyed by this subroutine utterly
c  p - the vector P in (1); destroyed by this subroutine utterly
c  q - the vector Q in (1); destroyed by this subroutine utterly
c  n - the dimensionality of all matrices and vectors involved
c  eps - the ABSOLUTE precision to which the eigenvalues will be
c       determined.
c
c                    Output parameters:
c
c  ier - error return code:
c    ier=0 means successful execution
c    ier=128 means that at some point the QR process failed to converge;
c      this should not happen, and most likely indicates a user bug.
c  diag - the eigenvalues of the matrix B in (1), in no particular
c      order
c  niter_max - the maximum number of iterations taken by the QR process
c      on any iteration; provided solely for a curious user's edification
c  aver_niter - the average number of iterations taken by the QR process;
c      provided solely for a curious user's edification
c       
c                    Work array:
c
c  w - must be at least 12*n+100 real *8 locations long
c

c
c       . . . find the first n-1 eigenvalues
c
        ier=0
        niter_max=0
c
        aver_niter=0
c
        do 2000 i=1,n-1
c
        nn=n-i+1
        call herm_p_rank1_one_dimen(jer,diag(i),supdiag(i),
     1      p(i),q(i),nn,eps,w,niter)
c
ccc        call prinf('deflating!!! i=',i,1)
c
        diag(i)=diag(i)+p(i)*conjg(q(i))
c
        aver_niter=aver_niter+niter
c
        if(jer .ne. 0) then
            ier=128
            return
        endif
c
        if(niter_max .lt. niter) niter_max=niter
 2000 continue
c
c        recover the last eigenvalue
c
        diag(n)=diag(n)+p(n)*conjg(q(n))
c
        aver_niter=aver_niter/n
c
        return
        end
c
c
c
c 
c
        subroutine herm_p_rank1_one_dimen(ier,diag,supdiag,
     1      p,q,n,eps,w,niter)
        implicit real *8 (a-h,o-z)
        complex *16 diag(*),supdiag(*),p(*),q(*),
     1      w(*),aa,bb,cc,clam,clam1,clam2,clamsum,d1,d2,s1,s2
c
c       conduct QR iterations
c
        ier=16
        ifout=0
        clamsum=0
        do 2000 ijk=1,1000
c
c       . . . find the shift
c
        d1=diag(1)+p(1)*conjg(q(1)) 
        d2=diag(2)+p(2)*conjg(q(2)) 
        s1=supdiag(1)+p(1)*conjg(q(2)) 
        s2=conjg(supdiag(1))+p(2)*conjg(q(1)) 
c
        aa=1
        bb=-(d1+d2)
        cc=d1*d2 - s1*s2
c
        clam1=(-bb+sqrt(bb**2-4*aa*cc)) / (2*aa)
        clam2=(-bb-sqrt(bb**2-4*aa*cc)) / (2*aa)
ccc        call prin2('clam1=',clam1,2)
ccc        call prin2('clam2=',clam2,2)
c
        clam=clam1
        if( abs(clam2-d1) .lt. abs(clam1-d1) ) clam=clam2
c
        clamsum=clamsum+clam
c
        do 1400 i=1,n
c
        diag(i)=diag(i)-clam
 1400 continue
c
c       conduct a QR ireration
c
        isubdiag=1
        lsubdiag=n+4
c
        iq2=isubdiag+lsubdiag
        lq2=n+4
c
        iuus=iq2+lq2
        luus=n*4+16
c
        ltot=iuus+luus
c
        call herm_p_rank1_iter(n,p,q,diag,supdiag,w(isubdiag),
     1      w(iq2),w(iuus))
c
c       check if the time has come to terminate the iterations
c
        d=abs(supdiag(1) + p(1)*conjg(q(2)))
ccc        call prin2('and the off-diagonal d=',d,1)
c
        if(d .lt. eps) ifout=ifout+1
        if(ifout .eq. 2) goto 2400
c
 2000 continue
c
        return
c
 2400 continue
c
        ier=0
        niter=ijk
c
        do 2600 i=1,n
        diag(i)=diag(i)+clamsum
 2600 continue
c
        return
        end
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     herm_p_rank1_iter_readable: unoptimized but readable
c     version of herm_p_rank1_iter. Kept for reference.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c$$$        subroutine herm_p_rank1_iter_readable(n,p,q,diag,
c$$$     1      supdiag,subdiag,q2,uus)
c$$$        implicit real *8 (a-h,o-z)
c$$$        complex *16 xy(2),uu(2,2),uus(2,2,*),p(*),q(*),
c$$$     1      diag(*),subdiag(*),supdiag(*),q2(*),subsub
c$$$c
c$$$        do 1200 i=1,n
c$$$        q2(i)=q(i)
c$$$ 1200 continue
c$$$        do 1210 i=1,n-1
c$$$        subdiag(i)=conjg(supdiag(i))
c$$$ 1210 continue
c$$$c
c$$$        ii=0
c$$$        do 2000 k=n,2,-1
c$$$        ii=ii+1
c$$$        xy(1)=supdiag(k-1)+p(k-1)*conjg(q(k))
c$$$        xy(2)=diag(k)+p(k)*conjg(q(k))
c$$$        call herm_p_rank1_rotfnd4(xy,uus(1,1,ii))
c$$$        if(k .eq. 2) goto 1300
c$$$        subsub=-q2(k)*conjg(p(k-2))
c$$$        call herm_p_rank1_two_elems_rotate(uus(1,1,ii),
c$$$     1      subdiag(k-2),subsub)
c$$$ 1300 continue
c$$$        call herm_p_rank1_two_elems_rotate(uus(1,1,ii),
c$$$     1      supdiag(k-1),diag(k))
c$$$        call herm_p_rank1_two_elems_rotate(uus(1,1,ii),
c$$$     1      diag(k-1),subdiag(k-1))
c$$$        call herm_p_rank1_two_elems_rotate(uus(1,1,ii),q2(k-1),q2(k))
c$$$        call herm_p_rank1_two_elems_rotate(uus(1,1,ii),p(k-1),p(k))
c$$$        dd1=abs(supdiag(k-1))**2 + abs(diag(k))**2
c$$$        dd2=abs(p(k-1)*conjg(q(k)))**2 + abs(p(k)*conjg(q(k)))**2
c$$$        if (dd2 .gt. dd1)
c$$$     1      p(k-1)=-supdiag(k-1)/conjg(q(k))
c$$$ 2000 continue
c$$$c
c$$$        iii=0
c$$$        do 3000 k=n,2,-1
c$$$        iii=iii+1
c$$$        uu(2,1)=conjg(uus(2,1,iii))
c$$$        uu(2,2)=conjg(uus(2,2,iii))
c$$$        uu(1,1)=conjg(uus(1,1,iii))
c$$$        uu(1,2)=conjg(uus(1,2,iii))
c$$$        supdiag(k-1)=-p(k-1)*conjg(q(k))
c$$$        call herm_p_rank1_two_elems_rotate(uus(1,1,iii),q(k-1),q(k))
c$$$        call herm_p_rank1_two_elems_rotate(uu,diag(k-1),supdiag(k-1))
c$$$        call herm_p_rank1_two_elems_rotate(uu,subdiag(k-1),diag(k))
c$$$ 3000 continue
c$$$        return
c$$$        end
c
c
c 
        subroutine herm_p_rank1_iter(n,p,q,diag,
     1      supdiag,subdiag,q2,uus)
        implicit real *8 (a-h,o-z)
        complex *16 xy(2),uu(2,2),uus(2,2,*),p(*),q(*),
     1      diag(*),subdiag(*),supdiag(*),q2(*),subsub
c
c        construct the vectors we will use
c
        do 1200 i=1,n
c
        q2(i)=q(i)
 1200 continue

        do 1210 i=1,n-1
c
        subdiag(i)=conjg(supdiag(i))
 1210 continue
c       
c        eliminate the superdiagonal
c
        ii=0
        do 2000 k=n,2,-1
c
        ii=ii+1
        xy(1)=supdiag(k-1)+p(k-1)*conjg(q(k))
        xy(2)=diag(k)+p(k)*conjg(q(k))
c
        call herm_p_rank1_rotfnd4(xy,uus(1,1,ii))
c
c
        if(k .eq. 2) goto 1300
c
        subsub=-q2(k)*conjg(p(k-2))
c
        call herm_p_rank1_two_elems_rotate10(uus(1,1,ii),
     1      subdiag(k-2),subsub)
 1300 continue
c
        call herm_p_rank1_two_elems_rotate(uus(1,1,ii),
     1      diag(k-1),subdiag(k-1))
c
c
        dd1=supdiag(k-1)*conjg(supdiag(k-1))
     1      + diag(k)*conjg(diag(k))
        dd2=p(k-1)*conjg(q(k))*conjg(p(k-1))*q(k)
     1      + p(k)*conjg(q(k))*conjg(p(k))*q(k)
c
c
        if (dd2 .le. dd1) then
c
        call herm_p_rank1_two_elems_rotate01(uus(1,1,ii),
     1      supdiag(k-1),diag(k))
        call herm_p_rank1_two_elems_rotate(uus(1,1,ii),p(k-1),p(k))
c
        else
c
        call herm_p_rank1_two_elems_rotate(uus(1,1,ii),
     1      supdiag(k-1),diag(k))
        call herm_p_rank1_two_elems_rotate01(uus(1,1,ii),p(k-1),p(k))
        p(k-1)=-supdiag(k-1)/conjg(q(k))
c
        endif
c
c
        call herm_p_rank1_two_elems_rotate10(uus(1,1,ii),q2(k-1),q2(k))
c
 2000 continue

c
c        restore the matrix to its tridiagonal form
c
        iii=0
        do 3000 k=n,2,-1
c
        iii=iii+1
c
        uu(2,1)=conjg(uus(2,1,iii))
        uu(2,2)=conjg(uus(2,2,iii))
        uu(1,1)=conjg(uus(1,1,iii))
        uu(1,2)=conjg(uus(1,2,iii))
c
c
        supdiag(k-1)=-p(k-1)*conjg(q(k))
c
        call herm_p_rank1_two_elems_rotate(uus(1,1,iii),q(k-1),q(k))

        call herm_p_rank1_two_elems_rotate(uu,diag(k-1),supdiag(k-1))
        call herm_p_rank1_two_elems_rotate01(uu,subdiag(k-1),diag(k))
c
 3000 continue
c
        return
        end
c 
c 
c 
c 
c 
        subroutine herm_p_rank1_two_elems_rotate(u,x,y)
        implicit complex *16 (a-h,o-z)
        dimension u(2,2)
c
        d1=u(1,1)*x+u(1,2)*y
        y=u(2,1)*x+u(2,2)*y
c 
        x=d1
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine herm_p_rank1_two_elems_rotate01(u,x,y)
        implicit complex *16 (a-h,o-z)
        dimension u(2,2)
c
        y=u(2,1)*x+u(2,2)*y
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine herm_p_rank1_two_elems_rotate10(u,x,y)
        implicit complex *16 (a-h,o-z)
        dimension u(2,2)
c
        x=u(1,1)*x+u(1,2)*y
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine herm_p_rank1_rotfnd4(a,u)
        implicit complex *16 (a-h,o-z)
        dimension a(2),u(2,2)
        real *8 d,d2,done
        data done/1.0d0/
c 
        alpha=a(1)
        beta=a(2)
c 
        d=(a(1)*conjg(a(1))+a(2)*conjg(a(2)))
        d=sqrt(d)
c 
        if(d .eq. 0) then
c 
            u(2,2)=1
            u(1,2)=0
            u(1,1)=1
            u(2,1)=0
            return
        endif
c
        d2=done/d
        u(1,1)=beta*d2
        u(2,2)=conjg(u(1,1))
        u(1,2)=-alpha*d2
        u(2,1)=-conjg(u(1,2))
c
        return
        end
c 
c 
c 
c 
c 
        subroutine herm_p_rank1_rarrcopy(x,y,n)
        implicit real *8 (a-h,o-z)
        real *8 x(*),y(*)
c
        do 2400 i=1,n
        y(i)=x(i)
 2400 continue
        return
        end  
c 
c 
c 
c 
c  
        subroutine herm_p_rank1_mach_zero(zero_mach)
        implicit real *8 (a-h,o-z)
c
        zero_mach=100       
c
        d1=1.1d0
        d=1.11d0
        do 1200 i=1,10000
c
        d=d/2
        d2=d1+d
        d4=sin(sin(d1)-sin(d2))
c
        if(d4 .eq. 0) goto 1400
c
 1200 continue
 1400 continue
c
        zero_mach=d
        return
        end
c
c
c
c
c
        subroutine herm_p_rank1_set_zero(x,n)
        implicit real *8 (a-h,o-z)
        real *8 x(n)
c 
        do 2000 i=1,n
        x(i)=0
 2000 continue
        return
        end
