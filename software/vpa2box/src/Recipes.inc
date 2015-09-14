!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! All of these routines are (C) Copr. 1986-92 Numerical Recipes Software #Y9%012#(.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!--------------------------------------------------------------------------------      
      SUBROUTINE AMOEBA(P,Y,MP,NP,NDIM,FTOL,FUNK,ITER)
!--------------------------------------------------------------------------------      
!     Modified by Clay Porch (4/6/94) to allow control of parameter bounds from
! subroutine FUNK
!     Multidimensional minimization of the function FUNK(X,NP) where X is      
! NDIM-dimensional vector by the downhill simplex method of Nelder and Mead
! (1965) as implemented in Numerical Recipes.
!     P(NDIM+1,NDIM) is the matrix of NDIM+1 vertices for the NDIM number of      
! of parameters.  Y(NDIM+1) is the corresponding value of FUNK for each set of
! vertices.  PBAR is the vector average of all points in the simplex other than
! the highest (worst)-- P(IHI,J).  
      INTEGER ITER, NDIM, NMAX, NP, MPTS, ILO, IHI, INHI, I,ITMAX,J,mp
      REAL (KIND=8) :: ALPHA,BETA,GAMMA
      PARAMETER (NMAX=800,ALPHA=1.D0,BETA=0.5d0,GAMMA=2.D0,ITMAX=500)
      REAL (KIND=8) :: P(MP,NP),Y(MP),PR(NMAX),PRR(NMAX),PBAR(NMAX),RTOL,YPR,YPRR,FTOL,DENOM
      EXTERNAL FUNK
      MPTS=NDIM+1
      ITER=0
1     ILO=1
      IF(Y(1).GT.Y(2))THEN
        IHI=1
        INHI=2
      ELSE
        IHI=2
        INHI=1
      ENDIF
      DO 11 I=1,MPTS
        IF(Y(I).LT.Y(ILO)) ILO=I
        IF(Y(I).GT.Y(IHI))THEN
          INHI=IHI
          IHI=I
        ELSE IF(Y(I).GT.Y(INHI))THEN
          IF(I.NE.IHI) INHI=I
        ENDIF
11    CONTINUE
      DENOM=DABS(Y(IHI))+DABS(Y(ILO))
      IF(DENOM.EQ.0.) RETURN
      RTOL=2.d0*dABS(Y(IHI)-Y(ILO))/DENOM
      IF(RTOL.LT.FTOL) RETURN
      IF(ITER.EQ.ITMAX) return
      ITER=ITER+1
      DO 12 J=1,NDIM
        PBAR(J)=0.d0
12    CONTINUE
      DO 14 I=1,MPTS
        IF(I.NE.IHI)THEN
          DO 13 J=1,NDIM
            PBAR(J)=PBAR(J)+P(I,J)
13        CONTINUE
        ENDIF
14    CONTINUE
      DO 15 J=1,NDIM
        PBAR(J)=PBAR(J)/NDIM
        PR(J)=(1.D0+ALPHA)*PBAR(J)-ALPHA*P(IHI,J)
15    CONTINUE
      CALL FUNK(PR,YPR,NP)
      IF(YPR.LE.Y(ILO))THEN
        DO 16 J=1,NDIM
          PRR(J)=GAMMA*PR(J)+(1.D0-GAMMA)*PBAR(J)
16      CONTINUE
        CALL FUNK(PRR,YPRR,NP)
        IF(YPRR.LT.Y(ILO))THEN
          DO 17 J=1,NDIM
            P(IHI,J)=PRR(J)
17        CONTINUE
          Y(IHI)=YPRR
        ELSE
          DO 18 J=1,NDIM
            P(IHI,J)=PR(J)
18        CONTINUE
          Y(IHI)=YPR
        ENDIF
      ELSE IF(YPR.GE.Y(INHI))THEN
        IF(YPR.LT.Y(IHI))THEN
          DO 19 J=1,NDIM
            P(IHI,J)=PR(J)
19        CONTINUE
          Y(IHI)=YPR
        ENDIF
        DO 21 J=1,NDIM
          PRR(J)=BETA*P(IHI,J)+(1.D0-BETA)*PBAR(J)
21      CONTINUE
        CALL FUNK(PRR,YPRR,NP)
        IF(YPRR.LT.Y(IHI))THEN
          DO 22 J=1,NDIM
            P(IHI,J)=PRR(J)
22        CONTINUE
          Y(IHI)=YPRR
        ELSE
          DO 24 I=1,MPTS
            IF(I.NE.ILO)THEN
!             {this set of loops differs from Numerical Recipes in that the
!              assignment of P(I,J) to PR(J) is made after the evaluation of
!              FUNK, which can have parameter limits}
              DO 23 J=1,NDIM
23              PR(J)=0.5d0*(P(I,J)+P(ILO,J))
              CALL FUNK(PR,Y(i),NP)
              DO 26 J=1,NDIM  
26              P(I,J)=PR(J)
            ENDIF
24        CONTINUE
        ENDIF
      ELSE
        DO 25 J=1,NDIM
          P(IHI,J)=PR(J)
25      CONTINUE
        Y(IHI)=YPR
      ENDIF
      GO TO 1
      END

!--------------------------------------------------------------------------------      
      SUBROUTINE newt(x,n,check,funcv)
! This is a quasi-newton method which provides a more global search than the 
! regular Newton iteration.  It is very likely to converge to the true answer
! provided the variates x are scaled to the order of unity.  However, it is
! also much slower than the Newton method.  I (Clay Porch) recommend trying
! the Newton method first if a reasonably good first quess can be provided 
! and then resorting to this routine if Newton's method does not converge
! properly.  Note that I have made several modifications to accomodate
! Fortran 90 and have also allowed different functions to be called by including
! funcv in the list of call variables (rather than forcing the user to name
! a single function 'funcv' as originally required). 
!--------------------------------------------------------------------------------      
      INTEGER n,nn,NP,MAXITS
      LOGICAL check
      DOUBLE PRECISION x(n),fvec(n),TOLF,TOLMIN,TOLX,STPMX
      PARAMETER (NP=40,MAXITS=200,TOLF=1.d-6,TOLMIN=1.d-6,TOLX=3.d-6,STPMX=100.d0)
!U    USES funcv,fmin,lnsrch,lubksb,ludcmp
      INTEGER i,its,j,indx(NP)
      DOUBLE PRECISION d,den,f,fold,stpmax,sum,temp,test,fjac(NP,NP),g(NP),p(NP),xold(NP),fmin
      EXTERNAL fmin,funcv
      nn=n ; f=fmin(x,fvec,nn,funcv) ; test=0.d0
      do 11 i=1,n
        if(abs(fvec(i)).gt.test)test=abs(fvec(i))
11    continue
      if(test.lt..01d0*TOLF) return
      sum=0.d0
      do 12 i=1,n
        sum=sum+x(i)**2
12    continue
      stpmax=STPMX*max(sqrt(sum), dble(n))
      do 21 its=1,MAXITS
        call funcv(n,x,fvec,fjac,its)
        do 14 i=1,n
          sum=0.d0
          do 13 j=1,n
            sum=sum+fjac(j,i)*fvec(j)
13        continue
          g(i)=sum
14      continue
        do 15 i=1,n
          xold(i)=x(i)
15      continue
        fold=f
        do 16 i=1,n
          p(i)=-fvec(i)
16      continue
        call ludcmp(fjac,n,NP,indx,d) ; call lubksb(fjac,n,NP,indx,p)
        call lnsrch(n,xold,fold,g,p,x,f,stpmax,check,fmin,funcv)
        test=0.d0
        do 17 i=1,n
          if(abs(fvec(i)).gt.test)test=abs(fvec(i))
17      continue
        if(test.lt.TOLF)then ; check=.false. ; return ; endif
        if(check)then
          test=0.d0 ; den=max(f,.5d0*n)
          do 18 i=1,n
            temp=abs(g(i))*max(abs(x(i)),1.d0)/den
            if(temp.gt.test)test=temp
18        continue
          if(test.lt.TOLMIN)then ; check=.true. ; else ; check=.false. ; endif
          return
        endif
        test=0.d0
        do 19 i=1,n
          temp=(abs(x(i)-xold(i)))/max(abs(x(i)),1.d0)
          if(temp.gt.test)test=temp
19      continue
        if(test.lt.TOLX)return
21    continue
      END
      
      SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,stpmax,check,func,funcv)
      INTEGER n,i
      LOGICAL check
      DOUBLE PRECISION f,fold,stpmax,g(n),p(n),x(n),xold(n),func,ALF,TOLX,FVEC(n)
      PARAMETER (ALF=1.d-4,TOLX=3.d-16)
      EXTERNAL func, funcv
      DOUBLE PRECISION a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam
      check=.false. ; sum=0.d0
      do 11 i=1,n
        sum=sum+p(i)*p(i)
11    continue
      sum=sqrt(sum)
      if(sum.gt.stpmax)then
        do 12 i=1,n
          p(i)=p(i)*stpmax/sum
12      continue
      endif
      slope=0.d0
      do 13 i=1,n
        slope=slope+g(i)*p(i)
13    continue
      test=0.d0
      do 14 i=1,n
        temp=abs(p(i))/max(abs(xold(i)),1.d0)
        if(temp.gt.test)test=temp
14    continue
      alamin=TOLX/test ; alam=1.d0
1     continue
        do 15 i=1,n
          x(i)=xold(i)+alam*p(i)
15      continue
      f=func(x,fvec,n,funcv)
        if(alam.lt.alamin)then
          do 16 i=1,n
            x(i)=xold(i)
16        continue
          check=.true.
          return
        else if(f.le.fold+ALF*alam*slope)then
          return
        else
          if(alam.eq.1.d0)then
            tmplam=-slope/(2.d0*(f-fold-slope))
          else
            rhs1=f-fold-alam*slope ; rhs2=f2-fold2-alam2*slope
            a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2) ; b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
            if(a.eq.0.d0)then
              tmplam=-slope/(2.d0*b)
            else
              disc=b*b-3.d0*a*slope
              tmplam=(-b+sqrt(disc))/(3.d0*a)
            endif
            if(tmplam.gt..5d0*alam)tmplam=.5d0*alam
          endif
        endif
      alam2=alam ; f2=f ; fold2=fold ; alam=max(tmplam,.1d0*alam)
      goto 1
      END

!--------------------------------------------------------------------------------      
      SUBROUTINE lubksb(a,n,np,indx,b)
!--------------------------------------------------------------------------------      
      INTEGER n,np,indx(n)
      DOUBLE PRECISION a(np,np),b(n)
      INTEGER i,ii,j,ll
      DOUBLE PRECISION sum
      ii=0
      do 12 i=1,n
        ll=indx(i) ; sum=b(ll) ; b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.d0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return ; END
      
!--------------------------------------------------------------------------------      
      SUBROUTINE ludcmp(a,n,np,indx,d)
!--------------------------------------------------------------------------------      
      INTEGER n,np,indx(n),NMAX
      DOUBLE PRECISION d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0d-20)
      INTEGER i,imax,j,k
      DOUBLE PRECISION aamax,dum,sum,vv(NMAX)
      d=1.d0
      do 12 i=1,n
        aamax=0.d0
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
!        if (aamax.eq.0.d0) pause 'singular matrix in ludcmp'
        vv(i)=1.d0/aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.d0
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k) ; a(imax,k)=a(j,k) ; a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.d0)a(j,j)=TINY
        if(j.ne.n)then
          dum=1.d0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return ; END
!--------------------------------------------------------------------------------      
      FUNCTION fmin(x,fvec,n,funcv)
!--------------------------------------------------------------------------------      
      INTEGER n
      DOUBLE PRECISION fmin,x(n),fvec(n)
!U    USES funcv
      INTEGER i
      DOUBLE PRECISION sum,DF
      EXTERNAL FUNCV
      call funcv(n,x,fvec,df,1)
      sum=0.d0
      do 11 i=1,n
        sum=sum+fvec(i)**2
11    continue
      fmin=0.5d0*sum
      return ; END


      SUBROUTINE SORT2(N,RA,RB)
      INTEGER IR,L,N,J,I
      REAL (KIND=8) :: RA(N),RB(N),RRA,RRB
      L=N/2+1 ; IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1 ; RRA=RA(L) ; RRB=RB(L)
        ELSE
          RRA=RA(IR) ; RA(IR)=RA(1) ; RRB=RB(IR) ; RB(IR)=RB(1) ;  IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA ; RB(1)=RRB ; RETURN
          ENDIF
        ENDIF
        I=L ; J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1)) J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J) ; RB(I)=RB(J) ; I=J ; J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA ; RB(I)=RRB
      GO TO 10
      END

!---------------------------------------------------------------------
    SUBROUTINE gaussj(a,n,np,b,m,mp)
! Linear equation solution by Gauss-Jordan elimination. a(1:n,1:n) is
! an input matrix representing the left side A of AX=B. b(1:n,1:m) is
! an input matrix representing B of the same equation.  On output A is
! replaced by its matrix inverse and b is replaced by the corresponding
! solution.  Note that if B is the identity matrix on input then on output
! both a and b will contain the matrix inverse.  NMAX is the largest
! anticipated value of {n,np,mp,m}
!---------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: m,mp,n,np,NMAX
      REAL(KIND=8) :: a(np,np),b(np,mp)
      PARAMETER (NMAX=5000)
      INTEGER :: i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      REAL(KIND=8) :: big,dum,pivinv
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                WRITE(*,*) 'singular matrix in gaussj'
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.) WRITE(*,*) 'singular matrix in gaussj'
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software #Y9%012#(.

