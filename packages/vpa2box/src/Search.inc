!---------------------------------------------------------------
      SUBROUTINE ESTIMATE2 (p,step,nop,func,max,iprint,stopcr,nloop,iquad,simp,var,functn,ifault)
!---------------------------------------------------------------
!     A PROGRAM FOR FUNCTION MINIMIZATION USING THE SIMPLEX METHOD.
!     The minimum found will often be a local, not a global, minimum.
!     Note that the present routine has been modified at line
!     to allow runs reaching the maximum iteration to be reported
! 
!     FOR DETAILS, SEE NELDER & MEAD, THE COMPUTER JOURNAL, JANUARY 1965
! 
!     PROGRAMMED BY D.E.SHAW,
!     CSIRO, DIVISION OF MATHEMATICS & STATISTICS
!     P.O. BOX 218, LINDFIELD, N.S.W. 2070
! 
!     WITH AMENDMENTS BY R.W.M.WEDDERBURN
!     ROTHAMSTED EXPERIMENTAL STATION
!     HARPENDEN, HERTFORDSHIRE, ENGLAND
! 
!     Further amended by Alan Miller,
!     CSIRO, Division of Mathematics & Statistics
!     Private Bag 10, CLAYTON, VIC. 3168
! 
!     ARGUMENTS:-
!     P()     = INPUT, STARTING VALUES OF PARAMETERS
!               OUTPUT, FINAL VALUES OF PARAMETERS
!     STEP()  = INPUT, INITIAL STEP SIZES
!     NOP     = INPUT, NO. OF PARAMETERS, INCL. ANY TO BE HELD FIXED
!     FUN!    = OUTPUT, THE FUNCTION VALUE CORRESPONDING TO THE FINAL
!               PARAMETER VALUES
!     MAX     = INPUT, THE MAXIMUM NO. OF FUNCTION EVALUATIONS ALLOWED
!     IPRINT  = INPUT, PRINT CONTROL PARAMETER
!                     < 0 NO PRINTING
!                     = 0 PRINTING OF PARAMETER VALUES AND THE FUNCTION
!                         VALUE AFTER INITIAL EVIDENCE OF CONVERGENCE.
!                     > 0 AS FOR IPRINT = 0 PLUS PROGRESS REPORTS AFTER
!                         EVERY IPRINT EVALUATIONS, PLUS PRINTING FOR TH
!                         INITIAL SIMPLEX.
!     STOPCR  = INPUT, STOPPING CRITERION
!     NLOOP   = INPUT, THE STOPPING RULE IS APPLIED AFTER EVERY NLOOP
!               FUNCTION EVALUATIONS.
!     IQUAD   = INPUT, = 1 IF THE FITTING OF A QUADRATIC SURFACE IS REQU
!                      = 0 IF NOT
!     SIMP    = INPUT, CRITERION FOR EXPANDING THE SIMPLEX TO OVERCOME
!               ROUNDING ERRORS BEFORE FITTING THE QUADRATIC SURFACE.
!     VAR()   = OUTPUT, CONTAINS THE DIAGONAL ELEMENTS OF THE INVERSE OF
!               THE INFORMATION MATRIX.
!     FUNCTN  = INPUT, NAME OF THE USER'S SUBROUTINE - ARGUMENTS (P,FUNC
!               WHICH RETURNS THE FUNCTION VALUE FOR A GIVEN SET OF
!               PARAMETER VALUES IN ARRAY P.
!****   FUNCTN MUST BE DECLARED EXTERNAL IN THE CALLING PROGRAM.
!       IFAULT  = OUTPUT, = 0 FOR SUCCESSFUL TERMINATION
!                         = 1 IF MAXIMUM NO. OF FUNCTION EVALUATIONS EXC
!                         = 2 IF INFORMATION MATRIX IS NOT +VE SEMI-DEFI
!                         = 3 IF NOP < 1
!                         = 4 IF NLOOP < 1
! 
!       Advice on usage:
!       If the function minimized can be expected to be smooth in the vi
!       of the minimum, users are strongly urged to use the quadratic-su
!       fitting option.   This is the only satisfactory way of testing t
!       minimum has been found.   The value of SIMP should be set to at
!       1000 times the rounding error in calculating the fitted function
!       e.g. in double precision on a micro- or mini-computer with about
!       decimal digit representation of floating-point numbers, the roun
!       errors in calculating the objective function may be of the order
!       1.E-12 say in a particular case.   A suitable value for SIMP wou
!       be 1.E-08.   However, if numerical integration is required in th
!       calculation of the objective function, it may only be accurate t
!       1.E-05 and an appropriate value for SIMP would be about 0.1.
!       If the fitted quadratic surface is not +ve definite (and the fun
!       should be smooth in the vicinity of the minimum), it probably me
!       that the search terminated prematurely and you have not found th
!       minimum.
! 
!       N.B. P, STEP AND VAR (IF IQUAD = 1) MUST HAVE DIMENSION AT LEAST
!            IN THE CALLING PROGRAM.
!       THE DIMENSIONS BELOW ARE FOR A MAXIMUM OF 100 PARAMETERS.
!      The dimension of BMAT should be at least NOP*(NOP+1)/2.
! 
!****      N.B. This version is in DOUBLE PRECISION throughout
! 
!       LATEST REVISION - 11 August 1991
! 
!***********************************************************************
! 
      implicit double precision (a-h,o-z)
      external functn
      dimension p(nop), step(nop), var(nop)
      dimension g(401,400), h(401), pbar(400), pstar(400), pstst(400),aval(400), bmat(600000), pmin(400), vc(600000), temp(400)
      data zero/0.d0/,one/1.d0/,two/2.d0/,three/3.d0/,half/0.5d0/
! 
!     A = REFLECTION COEFFICIENT, B = CONTRACTION COEFFICIENT, AND
!     C = EXPANSION COEFFICIENT.
! 
      data a,b,c/1.d0,0.5d0,2.d0/
! 
!     SET LOUT = LOGICAL UNIT NO. FOR OUTPUT
! 
      data lout/12/
! 
!     IF PROGRESS REPORTS HAVE BEEN REQUESTED, PRINT HEADING
! 
      if (iprint.gt.0) write (lout,870) iprint
      if (iprint.gt.0) write (*,880)
! 
!     CHECK INPUT ARGUMENTS
! 
      ifault=0
      if (nop.le.0) ifault=3
      if (nloop.le.0) ifault=4
      if (ifault.ne.0) return
! 
!     SET NAP = NO. OF PARAMETERS TO BE VARIED, I.E. WITH STEP.NE.0
! 
      nap=0
      loop=0
      iflag=0
      do 10 i=1,nop
	 if (step(i).ne.zero) nap=nap+1
10    continue
! 
!     IF NAP = 0 EVALUATE FUNCTION AT THE STARTING POINT AND RETURN
! 
      if (nap.gt.0) go to 20
      call functn (p,func,nop)
      return
! 
!     SET UP THE INITIAL SIMPLEX
! 
20    do 30 i=1,nop
30       g(1,i)=p(i)
      irow=2
      do 50 i=1,nop
	 if (step(i).eq.zero) go to 50
	 do 40 j=1,nop
40          g(irow,j)=p(j)
	 g(irow,i)=p(i)+step(i)
	 irow=irow+1
50    continue
      np1=nap+1
      neval=0
      do 70 i=1,np1
	 do 60 j=1,nop
60          p(j)=g(i,j)
	 call functn (p,h(i),nop)
	 neval=neval+1
	 if (iprint.le.0) go to 70
	 write (lout,890) neval,h(i),(p(j),j=1,nop)
	 write (*,900) neval,h(i),(p(j),j=1,nop)
70    continue
! 
!     START OF MAIN CYCLE.
! 
!     FIND MAX. & MIN. VALUES FOR CURRENT SIMPLEX (HMAX & HMIN).
! 
80    loop=loop+1
      imax=1
      imin=1
      hmax=h(1)
      hmin=h(1)
      do 100 i=2,np1
	 if (h(i).le.hmax) go to 90
	 imax=i
	 hmax=h(i)
	 go to 100
90       if (h(i).ge.hmin) go to 100
	 imin=i
	 hmin=h(i)
100   continue
! 
!     FIND THE CENTROID OF THE VERTICES OTHER THAN P(IMAX)
! 
      do 110 i=1,nop
110      pbar(i)=zero
      do 130 i=1,np1
	 if (i.eq.imax) go to 130
	 do 120 j=1,nop
120         pbar(j)=pbar(j)+g(i,j)
130   continue
      do 140 j=1,nop
	 fnap=nap
140      pbar(j)=pbar(j)/fnap
! 
!     REFLECT MAXIMUM THROUGH PBAR TO PSTAR,
!     HSTAR = FUNCTION VALUE AT PSTAR.
! 
      do 150 i=1,nop
150      pstar(i)=a*(pbar(i)-g(imax,i))+pbar(i)
      call functn (pstar,hstar,nop)
      neval=neval+1
      if (iprint.le.0) go to 160
      if (mod(neval,iprint).eq.0) write (lout,890) neval,hstar,(pstar(j),j=1,nop)
      if (mod(neval,iprint).eq.0) write (*,900) neval,hstar,(pstar(j),j=1,nop)
! 
!     IF HSTAR < HMIN, REFLECT PBAR THROUGH PSTAR,
!     HSTST = FUNCTION VALUE AT PSTST.
! 
160   if (hstar.ge.hmin) go to 200
      do 170 i=1,nop
170      pstst(i)=c*(pstar(i)-pbar(i))+pbar(i)
      call functn (pstst,hstst,nop)
      neval=neval+1
      if (iprint.le.0) go to 180
      if (mod(neval,iprint).eq.0) write (lout,890) neval,hstst,(pstst(j),j=1,nop)
      if (mod(neval,iprint).eq.0) write (*,900) neval,hstst,(pstst(j),j=1,nop)
! 
!     IF HSTST < HMIN REPLACE CURRENT MAXIMUM POINT BY PSTST AND
!     HMAX BY HSTST, THEN TEST FOR CONVERGENCE.
! 
180   if (hstst.ge.hmin) go to 300
      do 190 i=1,nop
	 if (step(i).ne.zero) g(imax,i)=pstst(i)
190   continue
      h(imax)=hstst
      go to 320
! 
!     HSTAR IS NOT < HMIN.
!     TEST WHETHER IT IS < FUNCTION VALUE AT SOME POINT OTHER THAN
!     P(IMAX).   IF IT IS REPLACE P(IMAX) BY PSTAR & HMAX BY HSTAR.
! 
200   do 210 i=1,np1
	 if (i.eq.imax) go to 210
	 if (hstar.lt.h(i)) go to 300
210   continue
! 
!     HSTAR > ALL FUNCTION VALUES EXCEPT POSSIBLY HMAX.
!     IF HSTAR <= HMAX, REPLACE P(IMAX) BY PSTAR & HMAX BY HSTAR.
! 
      if (hstar.gt.hmax) go to 230
      do 220 i=1,nop
	 if (step(i).ne.zero) g(imax,i)=pstar(i)
220   continue
      hmax=hstar
      h(imax)=hstar
! 
!     CONTRACTED STEP TO THE POINT PSTST,
!     HSTST = FUNCTION VALUE AT PSTST.
! 
230   do 240 i=1,nop
240      pstst(i)=b*g(imax,i)+(1.d0-b)*pbar(i)
      call functn (pstst,hstst,nop)
      neval=neval+1
      if (iprint.le.0) go to 250
      if (mod(neval,iprint).eq.0) write (lout,890) neval,hstst,(pstst(j),j=1,nop)
      if (mod(neval,iprint).eq.0) write (*,900) neval,hstst,(pstst(j),j=1,nop)
! 
!     IF HSTST < HMAX REPLACE P(IMAX) BY PSTST & HMAX BY HSTST.
! 
250   if (hstst.gt.hmax) go to 270
      do 260 i=1,nop
	 if (step(i).ne.zero) g(imax,i)=pstst(i)
260   continue
      h(imax)=hstst
      go to 320
! 
!     HSTST > HMAX.
!     SHRINK THE SIMPLEX BY REPLACING EACH POINT, OTHER THAN THE CURRENT
!     MINIMUM, BY A POINT MID-WAY BETWEEN ITS CURRENT POSITION AND THE
!     MINIMUM.
! 
270   do 290 i=1,np1
	 if (i.eq.imin) go to 290
	 do 280 j=1,nop
	    if (step(j).ne.zero) g(i,j)=(g(i,j)+g(imin,j))*half
	    p(j)=g(i,j)
280      continue
	 call functn (p,h(i),nop)
	 neval=neval+1
	 if (iprint.le.0) go to 290
	 if (mod(neval,iprint).eq.0) write (lout,890) neval,h(i),(p(j),j=1,nop)
	 if (mod(neval,iprint).eq.0) write (*,900) neval,h(i),(p(j),j=1,nop)
290   continue
      go to 320
! 
!     REPLACE MAXIMUM POINT BY PSTAR & H(IMAX) BY HSTAR.
! 
300   do 310 i=1,nop
	 if (step(i).ne.zero) g(imax,i)=pstar(i)
310   continue
      h(imax)=hstar
! 
!     IF LOOP = NLOOP TEST FOR CONVERGENCE, OTHERWISE REPEAT MAIN CYCLE.
! 
320   if (loop.lt.nloop) go to 80
! 
!     CALCULATE MEAN & STANDARD DEVIATION OF FUNCTION VALUES FOR THE
!     CURRENT SIMPLEX.
! 
      hstd=zero
      hmean=zero
      do 330 i=1,np1
330      hmean=hmean+h(i)
      fnp1=np1
      hmean=hmean/fnp1
      do 340 i=1,np1
340      hstd=hstd+(h(i)-hmean)**2
      hstd=dsqrt(hstd/dble(np1))
! 
!     IF THE RMS > STOPCR, SET IFLAG & LOOP TO ZERO AND GO TO THE
!     START OF THE MAIN CYCLE AGAIN.
! 
      if (hstd.le.stopcr.or.neval.gt.max) go to 350
      iflag=0
      loop=0
      go to 80
!  
!     FIND THE CENTROID OF THE CURRENT SIMPLEX AND THE FUNCTION VALUE TH
!  
350   do 370 i=1,nop
	 if (step(i).eq.zero) go to 370
	 p(i)=zero
	 do 360 j=1,np1
360         p(i)=p(i)+g(j,i)
	 fnp1=np1
	 p(i)=p(i)/fnp1
370   continue
      call functn (p,func,nop)
      neval=neval+1
      if (iprint.le.0) go to 380
      if (mod(neval,iprint).eq.0) write (lout,890) neval,func,(p(j),j=1,nop)
      if (mod(neval,iprint).eq.0) write (*,900) neval,func,(p(j),j=1,nop)
!  
!     TEST WHETHER THE NO. OF FUNCTION VALUES ALLOWED, MAX, HAS BEEN
!     OVERRUN; IF SO, EXIT WITH IFAULT = 1.
!  
380   if (neval.le.max) go to 390
      ifault=1
      if (iprint.lt.0) return
      write (lout,910) max
      write (lout,920) hstd
      write (lout,930) (p(i),i=1,nop)
      write (lout,940) func
!     folling line added by Clay Porch, 7/15/96 
       goto 430
!      return
!  
!     CONVERGENCE CRITERION SATISFIED.
!     IF IFLAG = 0, SET IFLAG & SAVE HMEAN.
!     IF IFLAG = 1 & CHANGE IN HMEAN <= STOPCR THEN SEARCH IS COMPLETE.
!  
390   if (iprint.lt.0) go to 400
      write (lout,950)
      write (lout,930) (p(i),i=1,nop)
      write (lout,940) func
400   if (iflag.gt.0) go to 420
      iflag=1
410   savemn=hmean
      loop=0
      go to 80
420   if (abs(savemn-hmean).ge.stopcr) go to 410
      if (iprint.lt.0) go to 430
      write (lout,960) neval
      write (lout,970) (p(i),i=1,nop)
      write (lout,980) func
430   if (iquad.le.0) return
!-----------------------------------------------------------------
!  
!     QUADRATIC SURFACE FITTING
!  
      if (iprint.ge.0) write (lout,990)
!  
!     EXPAND THE FINAL SIMPLEX, IF NECESSARY, TO OVERCOME ROUNDING
!     ERRORS.
!  
      neval=0
      do 460 i=1,np1
440      test=abs(h(i)-func)
	 if (test.ge.simp) go to 460
	 do 450 j=1,nop
	    if (step(j).ne.zero) g(i,j)=(g(i,j)-p(j))+g(i,j)
	    pstst(j)=g(i,j)
450      continue
	 call functn (pstst,h(i),nop)
	 neval=neval+1
	 go to 440
460   continue
!  
!     FUNCTION VALUES ARE CALCULATED AT AN ADDITIONAL NAP POINTS.
!  
      do 480 i=1,nap
	 i1=i+1
	 do 470 j=1,nop
470         pstar(j)=(g(1,j)+g(i1,j))*half
	 call functn (pstar,aval(i),nop)
	 neval=neval+1
480   continue
!  
!     THE MATRIX OF ESTIMATED SECOND DERIVATIVES IS CALCULATED AND ITS
!     LOWER TRIANGLE STORED IN BMAT.
!  
      a0=h(1)
      do 510 i=1,nap
	 i1=i-1
	 i2=i+1
	 if (i1.lt.1) go to 510
	 do 500 j=1,i1
	    j1=j+1
	    do 490 k=1,nop
490            pstst(k)=(g(i2,k)+g(j1,k))*half
	    call functn (pstst,hstst,nop)
	    neval=neval+1
	    l=i*(i-1)/2+j
	    bmat(l)=two*(hstst+a0-aval(i)-aval(j))
500      continue
510   continue
      l=0
      do 520 i=1,nap
	 i1=i+1
	 l=l+i
	 bmat(l)=two*(h(i1)+a0-two*aval(i))
520   continue
!  
!     THE VECTOR OF ESTIMATED FIRST DERIVATIVES IS CALCULATED AND
!     STORED IN AVAL.
!  
      do 530 i=1,nap
	 i1=i+1
	 aval(i)=two*aval(i)-(h(i1)+three*a0)*half
530   continue
!  
!     THE MATRIX Q OF NELDER & MEAD IS CALCULATED AND STORED IN G.
!  
      do 540 i=1,nop
540      pmin(i)=g(1,i)
      do 550 i=1,nap
	 i1=i+1
	 do 550 j=1,nop
	 g(i1,j)=g(i1,j)-g(1,j)
550   continue
      do 560 i=1,nap
	 i1=i+1
	 do 560 j=1,nop
	 g(i,j)=g(i1,j)
560   continue
!  
!     INVERT BMAT
!  
      call syminv (bmat,nap,bmat,temp,nullty,ifault,rmax)
      if (ifault.ne.0) go to 570
      irank=nap-nullty
      go to 580
570   if (iprint.ge.0) write (lout,1000)
      ifault=2
      return
!  
!     BMAT*A/2 IS CALCULATED AND STORED IN H.
!  
580   do 620 i=1,nap
	 h(i)=zero
	 do 610 j=1,nap
	    if (j.gt.i) go to 590
	    l=i*(i-1)/2+j
	    go to 600
590         l=j*(j-1)/2+i
600         h(i)=h(i)+bmat(l)*aval(j)
610      continue
620   continue
!  
!     FIND THE POSITION, PMIN, & VALUE, YMIN, OF THE MINIMUM OF THE
!     QUADRATIC.
!  
      ymin=zero
      do 630 i=1,nap
630      ymin=ymin+h(i)*aval(i)
      ymin=a0-ymin
      do 640 i=1,nop
	 pstst(i)=zero
	 do 640 j=1,nap
640      pstst(i)=pstst(i)+h(j)*g(j,i)
      do 650 i=1,nop
650      pmin(i)=pmin(i)-pstst(i)
      if (iprint.lt.0) go to 660
      write (lout,1010) ymin,(pmin(i),i=1,nop)
      write (lout,1020)
!  
!     Calculate true function value at the minimum of the quadratic.
!  
660   neval=neval+1
      call functn (pmin,hstar,nop)
!  
!     If HSTAR < FUNC, replace search minimum with quadratic minimum.
!  
      if (hstar.ge.func) go to 680
      func=hstar
      do 670 i=1,nop
670      p(i)=pmin(i)
      write (lout,1030) func
!  
!     Q*BMAT*Q'/2 IS CALCULATED & ITS LOWER TRIANGLE STORED IN VC
!  
680   do 750 i=1,nop
	 do 720 j=1,nap
	    h(j)=zero
	    do 710 k=1,nap
	       if (k.gt.j) go to 690
	       l=j*(j-1)/2+k
	       go to 700
690            l=k*(k-1)/2+j
700            h(j)=h(j)+bmat(l)*g(k,i)*half
710         continue
720      continue
	 do 740 j=i,nop
	    l=j*(j-1)/2+i
	    vc(l)=zero
	    do 730 k=1,nap
730            vc(l)=vc(l)+h(k)*g(k,j)
740      continue
750   continue
!  
!     THE DIAGONAL ELEMENTS OF VC  ARE COPIED INTO VAR.
!  
      j=0
      do 760 i=1,nop
	 j=j+i
	 var(i)=vc(j)
760   continue
      if (iprint.lt.0) return
      write (lout,1040) irank
      ijk=1
      go to 830
770   continue
      write (lout,1050)
      call syminv (vc,nap,bmat,temp,nullty,ifault,rmax)
!  
!     BMAT NOW CONTAINS THE INFORMATION MATRIX
!  
      write (lout,1060)
      ijk=3
      go to 830
!  
!     Calculate correlations of parameter estimates, put into VC.
!  
780   ijk=2
      ii=0
      ij=0
      do 800 i=1,nop
	 ii=ii+i
	 if (vc(ii).gt.zero) then
	    vc(ii)=one/sqrt(vc(ii))
	 else
	    vc(ii)=zero
	 end if
	 jj=0
	 do 790 j=1,i-1
	    jj=jj+j
	    ij=ij+1
	    vc(ij)=vc(ij)*vc(ii)*vc(jj)
790      continue
	 ij=ij+1
800   continue
      write (lout,1070)
      write (lout,1071) (ii,ii=1,nop)
      ii=0
      do 810 i=1,nop
	 ii=ii+i
	 if (vc(ii).ne.zero) vc(ii)=one
810   continue
      go to 830
820   write (lout,1080) neval
      return
!  
!     Pseudo-subroutine to print VC if IJK = 1 or 2, or
!     BMAT if IJK = 3.
!  
830   l=1
840   if (l.gt.nop) go to (770,820,780),ijk
      ii=l*(l-1)/2
      do 860 i=l,nop
	 i1=ii+l
	 ii=ii+i
	 i2=min(ii,i1+19)
	 if (ijk.eq.3) go to 850
	 write (lout,1090) i,(vc(j),j=i1,i2)
	 go to 860
850      write (lout,1091) i,(bmat(j),j=i1,i2)
860   continue
      write (lout,1100)
      l=l+20
      go to 840
 

870   format (' PROGRESS REPORT EVERY',i4,' FUNCTION EVALUATIONS',/,' EVAL.  FUNC.',15x,'PARAMETER VALUES')
880   format (' EVAL.  FUNC.',15x,'PARAMETER VALUES')
890   format (/1x,i6,2x,g12.5,2x,5g12.5,3(/20x,5g12.5))
900   format (1x,i6,2x,g12.5,2x,5g12.5,3(/20x,5g12.5))
910   format (' NO. OF FUNCTION EVALUATIONS EXCEEDS',i5)
920   format (' RMS OF FUNCTION VALUES OF LAST SIMPLEX =',g14.6)
930   format (' CENTROID OF LAST SIMPLEX =',4(/1x,6g13.5))
940   format (' FUNCTION VALUE AT CENTROID =',g14.6)
950   format (/' EVIDENCE OF CONVERGENCE')
960   format (//' MINIMUM FOUND AFTER',i5,' FUNCTION EVALUATIONS')
970   format (' MINIMUM AT',4(/1x,6g13.6))
980   format (' FUNCTION VALUE AT MINIMUM =',g14.6)
990   format (/' QUADRATIC SURFACE FITTING ABOUT SUPPOSED MINIMUM'/)
1000  format (/' MATRIX OF ESTIMATED SECOND DERIVATIVES NOT +VE DEFN.'/'MINIMUM PROBABLY NOT FOUND'/)
1010  format (' MINIMUM OF QUADRATIC SURFACE =',g14.6,' AT',4(/1x,6g13.5))
1020  format (' IF THIS DIFFERS BY MUCH FROM THE MINIMUM ESTIMATED',1x, &
	 'FROM THE MINIMIZATION',/,' THE MINIMUM MAY BE FALSE OR THE INFORMATION MATRIX MAY BE',1x,'INACCURATE')
1030  format (' True func. value at minimum of quadratic = ',g14.6/)
1040  format (' RANK OF INFORMATION MATRIX =',i3/,' GENERALIZED INVERSE OF INFORMATION MATRIX:-')
1050  format (/,' IF THE FUNCTION MINIMIZED WAS -LOG(LIKELIHOOD)',/,' THIS IS THE COVARIANCE MATRIX OF THE PARAMETERS',/, &
	 'IF THE FUNCTION WAS A SUM OF SQUARES OF RESIDUALS',/,' THIS MATRIX MUST BE MULTIPLIED BY TWICE THE ESTIMATED', & 
	 1x,'RESIDUAL VARIANCE',/,' TO OBTAIN THE COVARIANCE MATRIX.')
1060  format (' INFORMATION MATRIX:-'/)
1070  format (/' CORRELATION MATRIX:-')
1071  format (7x,20i6)
1080  format (/' A FURTHER',i4,' FUNCTION EVALUATIONS HAVE BEEN USED'/)
1090  format (1x,i4,2x,20f13.4)
1091  format (1x,i4,2x,20f6.4)
1100  format (/)
      end


!---------------------------------------------------------------
      subroutine syminv (a,n,c,w,nullty,ifault,rmax)
!---------------------------------------------------------------
!  
!     ALGORITHM AS7, APPLIED STATISTICS, VOL.17, 1968.
!  
!     ARGUMENTS:-
!     A()     = INPUT, THE SYMMETRIC MATRIX TO BE INVERTED, STORED IN
!               LOWER TRIANGULAR FORM
!     N       = INPUT, ORDER OF THE MATRIX
!     C()     = OUTPUT, THE INVERSE OF A (A GENERALIZED INVERSE IF C IS
!               SINGULAR), ALSO STORED IN LOWER TRIANGULAR.
!               C AND A MAY OCCUPY THE SAME LOCATIONS.
!     W()     = WORKSPACE, DIMENSION AT LEAST N.
!     NULLTY  = OUTPUT, THE RANK DEFICIENCY OF A.
!     IFAULT  = OUTPUT, ERROR INDICATOR
!                     = 1 IF N < 1
!                     = 2 IF A IS NOT +VE SEMI-DEFINITE
!                     = 0 OTHERWISE
!     RMAX    = OUTPUT, APPROXIMATE BOUND ON THE ACCURACY OF THE DIAGONA
!               ELEMENTS OF C.  E.G. IF RMAX = 1.E-04 THEN THE DIAGONAL
!               ELEMENTS OF C WILL BE ACCURATE TO ABOUT 4 DEC. DIGITS.
!  
!     LATEST REVISION - 18 October 1985
!  
!***********************************************************************
!  
      implicit double precision (a-h,o-z)
      dimension a(*), c(*), w(n)
      data zero/0.d0/,one/1.d0/
!  
      nrow=n
      ifault=1
      if (nrow.le.0) go to 90
      ifault=0
!  
!     CHOLESKY FACTORIZATION OF A, RESULT IN C
!  
      call chola (a,nrow,c,nullty,ifault,rmax,w)
      if (ifault.ne.0) go to 90
!  
!     INVERT C & FORM THE PRODUCT (CINV)'*CINV, WHERE CINV IS THE INVERS
!     OF C, ROW BY ROW STARTING WITH THE LAST ROW.
!     IROW = THE ROW NUMBER, NDIAG = LOCATION OF LAST ELEMENT IN THE ROW
!  
      nn=nrow*(nrow+1)/2
      irow=nrow
      ndiag=nn
10    if (c(ndiag).eq.zero) go to 60
      l=ndiag
      do 20 i=irow,nrow
	 w(i)=c(l)
	 l=l+i
20    continue
      icol=nrow
      jcol=nn
      mdiag=nn
30    l=jcol
      x=zero
      if (icol.eq.irow) x=one/w(irow)
      k=nrow
40    if (k.eq.irow) go to 50
      x=x-w(k)*c(l)
      k=k-1
      l=l-1
      if (l.gt.mdiag) l=l-k+1
      go to 40
50    c(l)=x/w(irow)
      if (icol.eq.irow) go to 80
      mdiag=mdiag-icol
      icol=icol-1
      jcol=jcol-1
      go to 30
!  
!     Special case, zero diagonal element.
!  
60    l=ndiag
      do 70 j=irow,nrow
	 c(l)=zero
	 l=l+j
70    continue
!  
!      End of row.
!  
80    ndiag=ndiag-irow
      irow=irow-1
      if (irow.ne.0) go to 10
90    return
      end

!-----------------------------------------------------------------------------      
      subroutine chola (a,n,u,nullty,ifault,rmax,r)
!-----------------------------------------------------------------------------      
!  
!     ALGORITHM AS6, APPLIED STATISTICS, VOL.17, 1968, WITH
!     MODIFICATIONS BY A.J.MILLER
!  
!     ARGUMENTS:-
!     A()     = INPUT, A +VE DEFINITE MATRIX STORED IN LOWER-TRIANGULAR
!               FORM.
!     N       = INPUT, THE ORDER OF A
!     U()     = OUTPUT, A LOWER TRIANGULAR MATRIX SUCH THAT U*U' = A.
!               A & U MAY OCCUPY THE SAME LOCATIONS.
!     NULLTY  = OUTPUT, THE RANK DEFICIENCY OF A.
!     IFAULT  = OUTPUT, ERROR INDICATOR
!                     = 1 IF N < 1
!                     = 2 IF A IS NOT +VE SEMI-DEFINITE
!                     = 0 OTHERWISE
!     RMAX    = OUTPUT, AN ESTIMATE OF THE RELATIVE ACCURACY OF THE
!               DIAGONAL ELEMENTS OF U.
!     R()     = OUTPUT, ARRAY CONTAINING BOUNDS ON THE RELATIVE ACCURACY
!               OF EACH DIAGONAL ELEMENT OF U.
!  
!     LATEST REVISION - 18 October 1985
!  
!***********************************************************************
!  
      implicit double precision (a-h,o-z)
      dimension a(*), u(*), r(n)
!  
!     ETA SHOULD BE SET EQUAL TO THE SMALLEST +VE VALUE SUCH THAT
!     1.0 + ETA IS CALCULATED AS BEING GREATER THAN 1.0 IN THE ACCURACY
!     BEING USED.
!  
      data eta/1.d-16/,zero/0.d0/,five/5.d0/,tiny/1.d-20/
!  
      ifault=1
      if (n.le.0) go to 90
      ifault=2
      nullty=0
      rmax=eta
      r(1)=eta
      j=1
      k=0
!  
!     FACTORIZE COLUMN BY COLUMN, ICOL = COLUMN NO.
!  
      do 80 icol=1,n
	 l=0
!  
!     IROW = ROW NUMBER WITHIN COLUMN ICOL
!  
	 do 40 irow=1,icol
	    k=k+1
	    w=a(k)
	    if (irow.eq.icol) rsq=(w*eta)**2
	    m=j
	    do 10 i=1,irow
	       l=l+1
	       if (i.eq.irow) go to 20
	       w=w-u(l)*u(m)
	       if (irow.eq.icol) rsq=rsq+(u(l)**2*r(i))**2
	       m=m+1
10          continue
20          if (irow.eq.icol) go to 50
	    if (u(l).eq.zero) go to 30
	    u(k)=w/u(l)
	    go to 40
30          u(k)=zero
	    if (abs(w).gt.abs(rmax*a(k))) go to 90
40       continue
!  
!     END OF ROW, ESTIMATE RELATIVE ACCURACY OF DIAGONAL ELEMENT.
!  
50       rsq=sqrt(rsq)
	 if (abs(w).le.five*rsq) go to 60
	 if (w.lt.zero) w=tiny
	 u(k)=sqrt(w)
	 r(i)=rsq/w
	 if (r(i).gt.rmax) rmax=r(i)
	 go to 70
60       u(k)=zero
	 nullty=nullty+1
70       j=j+icol
80    continue
      ifault=0
!  
90    return
      end
