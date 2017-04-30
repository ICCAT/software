!
! Last change:  CEP  13 Apr 2017   12:29 pm
! program VPA-2BOX.F90 , VERSION 4.01
! A sequential population analysis tool
! programmed by Clay E. Porch in FORTRAN 90
! 11/1/2000
!
! changes for 3.01: 1) uses bisection routine for 1 dimensional problems rather than Newton's method
!                   2) corrects formula for chi-square discrepancy of lognormal deviates
!                   3) changes AIC and BIC formulae from deviance-based to full-likelihood expressions to
!                      allow comparison of models where the variances are estimated
!                   4) puts in error message to prevent use of overlap model when estimating terminal parameters on N
!                      scale but using different age groups for each stock (in that case cannot solve backwards recursion)
!                   5) some format changes (allowing big numbers to be represented in exponential format)
! changes for 4.0:  1) generalizes tag-attrition model to allow for pop-up satellite tags etc...
!                   2) uses new formats for the tagging data and tagging parameter definitions
!                      to make them easier to work with
!                   3) use of MFEXP(arg) function to handle large arguments that cause overflows in EXP(arg)
!                   4) added switch read from control file to make computing the covariance matrix optional
! changes for 4.01  1) option to prorate plus group out to older ages (makes equilibrium assumption in first year and
!                      then tracks each cohort)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! VARIABLE MODULE DEFINITIONS !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE DIMS
        INTEGER :: As, Ys, Gs, Bs, Ts  ! maximum age (or max. # length bins), maximum fisheries or surveys, maximum year blocks, maximum tag cohorts
        PARAMETER (As=101 , Ys=201, Gs=50 , Bs=2, Ts=10000)
      END MODULE

      MODULE PARAMETERS
        INTEGER NPAR
        REAL (KIND=8) PI
        PARAMETER (NPAR=2000, PI=3.1415926536)
        REAL (KIND=8) :: PARM_SPECS(NPAR,5),PARM_EST(NPAR),PARM_VAR(NPAR),PDEV,SIGMA_FT,SIGMA_REC,SIGMA_STOCK,RATIO_STOCK, &
                         CV_OVERIDE,SIGMA_TAG_MOD,FISHING_PRESSURE(12)
        INTEGER :: PARM_KEY(0:NPAR+1), OPTION(10), SEED, MAXITER, CHECKFLAG, MODEL_TYPE, SCALEs,ICALL=0,NBOX,LINK_FT,FIRSTCALL=0, &
                   LINK_REC,LINK_STOCK,PDF_TAG,ADD_VAR,N_RETRO,N_RETRO_LOOP,LINK_YOUNGEST,LINK_OLDEST,STINE_CORR,DATE_VALUES(8), &
                   COMPUTE_COVAR,BOOT_FILE_TYPE
        CHARACTER (LEN=120) :: INFILE(0:10),CONFILE
        CHARACTER (LEN=9) :: MONTH_NAME(12)
        CHARACTER (LEN=2) :: SPACE(0:9), ENDL='EN'
        DATA MONTH_NAME/' January',' February','    March','    April','      May','     June','     July','   August','September',&
                        ' October',' November',' December'/
        DATA SPACE/'0','1X','2X','3X','4X','5X','6X','7X','8X','9X'/
      END MODULE PARAMETERS

      MODULE DATUM
        USE DIMS
        REAL (KIND=8) :: CATCH_DATA(Bs,0:As,0:Ys),EFFORT_DATA(Gs,Bs,0:Ys),PSEL(Gs,Bs,0:As,0:Ys),  &
                         SIGMA_CATCH(Bs),SIGMA_EFFORT(Gs,Bs,0:Ys),SIGMA_M(Bs,0:As,0:Ys),SIGMA_T(Bs,0:As,0:Ys), &
                         SIGMA_F(Bs,0:Ys),SIGMA_FRATIO(Bs,0:Ys),SIGMA_TERMINAL(Bs,0:As), SIGMA_SCALE_VARIANCE(Gs,Bs), &
                         SIGMA_TAG_REPORT(1000),SIGMA_TAG_SURV(1000),SIGMA_TAG_LOSS(1000), &
                         SIGMA_Q(Gs,Bs,0:Ys),INDEX_MEAN(Gs,Bs),WEIGHT_CATCH(Gs,Bs,0:As,0:Ys),FECUNDITY(Bs,0:As), &
                         WEIGHT_SSB(Bs,0:As,0:Ys),RECOVERY_DATA(Ts,BS),RECAPTURE_DATA(Ts,BS,YS), &
                         EFFORT_DATA_STORE(Gs,Bs,0:Ys), SUMRDATA(Ts),SIGMA_TAG_NOMIX(1000), &
                         SIGMA_TAG_NOMIX2(1000),T_PROPORTION(12),PSEL_MEAN(Gs,Bs,0:Ys),SAMPLE_SIZE(Gs,Bs,0:Ys), &
                         PRED_WEIGHT(Gs,Bs,0:As,0:Ys), PRED_LENGTH(Gs,Bs,0:As,0:Ys)
        INTEGER :: PDF_CATCH(Bs),PDF_EFFORT(Gs,Bs),BIO_EFFORT(Gs,Bs),BIO_CATCH(Bs), &
                   PDF_STOCKRECRUIT,SEL_TYPE(Gs,Bs),PDF_FRATIO(Bs),PDF_TERMINAL(Bs),PDF_M(Bs),PDF_T(Bs), &
                   BOTH_STOCKS(Gs,Bs),IGNORE_RECRUIT(2),N_DATA,N_Qs,NPI(GS,BS)
        CHARACTER (LEN=50) :: TITLE_FISHERY(Bs), TITLE_EFFORT(Gs,Bs)
      END MODULE DATUM

      MODULE STATISTICS
        USE DIMS
        TYPE GRWTH
          INTEGER :: CURVE
          REAL (KIND=8) :: LINF,K,T0,M,WA,WB,DATE
        END TYPE GRWTH
        TYPE (GRWTH), SAVE :: GROWTH(2,0:100)
        TYPE TAGCOHORT
          INTEGER :: Br,Sr,Yr
          REAL (KIND=8) :: Mr
          INTEGER :: Ye
          REAL (KIND=8) :: Me
          INTEGER :: Ar
          REAL (KIND=8) :: N,Sigma
          INTEGER :: Surv,Shed,Rpt0,Rpt1,Rpt2,Mix1,Mix2
        END TYPE TAGCOHORT
        TYPE (TAGCOHORT), SAVE :: TC(Ts)
        REAL (KIND=8) :: N_STOCK(Bs,0:As,0:Ys),N_AREA(Bs,0:As,0:Ys),N_temp(Bs,Bs,0:As,0:Ys), &
                         F(Bs,0:As,0:Ys),M(Bs,0:As,0:Ys),T(Bs,0:As,0:Ys), &
                         EFFORT(Gs,Bs,0:Ys),SEL_EFFORT(Gs,Bs,0:As,0:Ys),Q_EFFORT(Gs,Bs,Ys), & 
                         PLUSAGE(Bs,0:Ys),STOCK_RECRUIT(Bs,5),SSB(Bs,Ys),RECRUITS(Bs,Ys),SCALE_VARIANCE(Gs,Bs), &
                         SEL_TERMINAL(Bs,0:As,2),F_RATIO(Bs,Ys),TAG_LOSS(1000),TAG_REPORT(1000,0:Bs),TAG_SURV(1000), &
                         N_TAGS(Ts,Bs,Ys),RECAPTURES(Ts,Bs,Ys),RECOVERIES(Ts,BS),T_mod_month(Ts), &
                         ETA_Q(Gs,Bs,Ys),ETA_F(Bs,Ys),ETA_SR(Bs,Ys)=0,ETA_T(Bs,0:As,Ys),ETA_M(Bs,0:As,Ys),ETA_TERMINAL(Bs,0:As), &
                         ETA_TAG_SURV(1000),ETA_TAG_LOSS(1000),ETA_TAG_REPORT(1000),TAG_NOMIX(1000),TAG_NOMIX2(1000), &
                         ETA_TAG_NOMIX(1000),ETA_SCALE_VARIANCE(Gs,Bs)
        INTEGER :: FIRSTYEAR, LASTYEAR, DISPLAYYEAR, NYEARS, FIRSTAGE, LASTAGE, YEAR, AGE, GEAR, BOX, &
                   PLUSGROUP,EXPANDED_PLUSGROUP,NGEARS(Bs),SEASON_EFFORT(Gs,Bs),SEASON_SSB(Bs),AGE_EFFORT(2,Gs,Bs), &
                   STARTYEAR(Bs),ENDYEAR(Bs),N_TAG_COHORTS,N_SURV_PARS,N_SHED_PARS,N_REPT_PARS,N_MIX_PARS

      END MODULE STATISTICS                     

      MODULE LOGLIKELIHOODS
        USE DIMS
        REAL (KIND=8) :: LIKE_RECRUITMENT(Bs),LIKE_EFFORT(Gs,Bs),LIKE_TAG(Bs),LIKE_M(Bs),LIKE_T(Bs), &
                         PENALTY,LIKELIHOODS,POSTERIORS,CONSTRAINTS,LIKE_Q(Gs,Bs),LIKE_F(Bs),LIKE_TERMINAL(Bs), &
                         LIKE_TAG_REPORT,LIKE_TAG_SURV,LIKE_TAG_LOSS,LIKE_TAG_NOMIX,LIKE_V(Bs), &
                         SUM_EFFORT_DISCREPANCY(Gs,Bs),EFFORT_DISCREPANCY(Gs,Bs,Ys),DISCREPANCY
        CHARACTER (LEN=20) :: PDFNAME(0:12)
        CHARACTER (LEN=1) :: CONSTANTS
        DATA PDFNAME/'Not used','Lognormal dist.','Normal dist.','Multinomial dist.','Poisson dist.','Chi-square stat.', &
                     'Robust stat.','gamma','unavailable','unavailable','unavailable','unavailable','std. nrml dev.'/
        DATA CONSTANTS/'N'/
      END MODULE LOGLIKELIHOODS                     
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!       MAIN PROGRAM        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      PROGRAM VPA
      USE PARAMETERS ; USE STATISTICS, ONLY: LASTYEAR
      IMPLICIT NONE
      CHARACTER (LEN=12) :: RETRO_RESULTS(20),RETRO_ESTIMATES(20),RETRO_SPREAD(20)
      DATA RETRO_RESULTS/'MINUS1.R','MINUS2.R','MINUS3.R','MINUS4.R','MINUS5.R','MINUS6.R','MINUS7.R','MINUS8.R', &
                         'MINUS9.R','MINUS10.R','MINUS11.R','MINUS12.R','MINUS13.R','MINUS14.R','MINUS15.R', &
                         'MINUS16.R','MINUS17.R','MINUS18.R','MINUS19.R','MINUS20.R'/
      DATA RETRO_ESTIMATES/'MINUS1.EST','MINUS2.EST','MINUS3.EST','MINUS4.EST','MINUS5.EST','MINUS6.EST','MINUS7.EST','MINUS8.EST',&
                         'MINUS9.EST','MINUS10.EST','MINUS11.EST','MINUS12.EST','MINUS13.EST','MINUS14.EST','MINUS15.EST', &
                         'MINUS16.EST','MINUS17.EST','MINUS18.EST','MINUS19.EST','MINUS20.EST'/
      DATA RETRO_SPREAD/'MINUS1.SPD','MINUS2.SPD','MINUS3.SPD','MINUS4.SPD','MINUS5.SPD','MINUS6.SPD','MINUS7.SPD','MINUS8.SPD', &
                         'MINUS9.SPD','MINUS10.SPD','MINUS11.SPD','MINUS12.SPD','MINUS13.SPD','MINUS14.SPD','MINUS15.SPD', &
                         'MINUS16.SPD','MINUS17.SPD','MINUS18.SPD','MINUS19.SPD','MINUS20.SPD'/
      WRITE(*,'(/,/,/,/,/,/,/,/,/,/,/,/,/,/,/,/,/,/,/,/,/,/)')
      WRITE(*,*) '+-------------------------------------------------------------------------+'
      WRITE(*,*) '| Program VPA-2BOX.F90, Version 4.01 (April 15, 2017)                        |'
      WRITE(*,*) '|                                                                         |'
      WRITE(*,*) '| A virtual population analysis tool that uses catch-at-age, indices of   |'
      WRITE(*,*) '| abundance, indices of mortality rates, and tag-recoveries to estimate   |'
      WRITE(*,*) '| the abundance and mortality of one or two (intermixing) populations.    |'
      WRITE(*,*) '|                                                                         |'
      WRITE(*,*) '| based on the methods of                                                 |'
      WRITE(*,*) '|   Porch, C. E., Turner, S. C., and Powers, J. E. 2001                   |'
      WRITE(*,*) '|   Virtual population analyses of Atlantic bluefin tuna with alternative |'
      WRITE(*,*) '|   models of transatlantic migration: 1970-1997. Int. Comm. Conserv. Atl.|'
      WRITE(*,*) '|   Tunas, Coll. Vol. Sci. Pap. 52: 1022-1045                             |'
      WRITE(*,*) '|                                                                         |'
      WRITE(*,*) '| programmed by                                                           |'
      WRITE(*,*) '|   Clay E. Porch                                                         |'
      WRITE(*,*) '|   NOAA Fisheries                                                        |'
      WRITE(*,*) '|   75 Virginia Beach Drive                                               |'
      WRITE(*,*) '|   Miami, Fl 33149 (USA)                                                 |'
      WRITE(*,*) '+-------------------------------------------------------------------------+'
      WRITE(*,'(/,/)')

      OPEN(12,FILE='VPA-2BOX.LOG',STATUS='UNKNOWN')
      CALL CONTROL ! options controlling simulation
      CALL READ_DATA  ! user defined
      CALL READ_PARAMETERS
      CALL ESTIMATE(1+COMPUTE_COVAR)
      CALL RESULTS(1)  ! user defined
      WRITE(12,*) 'DETERMINISTIC RUN COMPLETED'
      IF(N_RETRO>0) THEN
          WRITE(*,'(/,/,/,/,/)') ; WRITE(*,*) '*******************************************'
          WRITE(*,*) 'RETROSPECTIVE PATTERN ANALYSES' ; WRITE(*,*) '*******************************************' ; WRITE(*,'(/)')
          DO N_RETRO_LOOP=1,N_RETRO
            CALL REDEFINE_PARAMETERS
            LASTYEAR=LASTYEAR-1 ; ICALL=0
            INFILE(5)=RETRO_SPREAD(N_RETRO_LOOP) ; INFILE(4)=RETRO_ESTIMATES(N_RETRO_LOOP) ; INFILE(3)=RETRO_RESULTS(N_RETRO_LOOP)
            CALL ESTIMATE(0)
            CALL RESULTS(0)  ! user defined
          END DO
          WRITE(12,*) 'RETROSPECTIVE RUNS COMPLETED'
        ELSE
          CALL BOOTSTRAP
          WRITE(12,*) 'BOOTSTRAP RUNS COMPLETED '
      ENDIF
!     CALL PROFILES  ! user defined
      CLOSE(12)
      STOP
      END
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SUBROUTINES AND FUNCTIONS !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!_________________________________________________________________
!     User prescribed subroutines    
!_________________________________________________________________      
      INCLUDE 'vpa-2box.inc'
      INCLUDE 'RECIPES.INC'
!_________________________________________________________________
      SUBROUTINE CONTROL
! reads the control information
!_________________________________________________________________      
      USE PARAMETERS ; USE DATUM, ONLY : PDF_STOCKRECRUIT, IGNORE_RECRUIT, T_PROPORTION
      IMPLICIT NONE
      INTEGER :: IDUMMY,I,LINE=0
      REAL (KIND=8) :: DUMMY
      CHARACTER (LEN=1) :: CH,CH2,CH3
      COMPUTE_COVAR=1
      CALL GETCL(CONFILE)
      IF(CONFILE=='' .OR. CONFILE==' ') THEN
        WRITE(*,*) 'ENTER THE NAME OF THE CONTROL FILE: ' ; READ(*,*) CONFILE
      ENDIF
1     OPEN(10,FILE=CONFILE,STATUS='OLD',ERR=10)
      DO I=0,6
2       READ(10,'(3A1)')  CH,CH2,CH3 ; LINE=LINE+1
        IF(CH.NE.'#'.AND.CH.NE.'!'.AND.CH.NE.'_'.AND.CH.NE.'@') THEN
            IF(CH2=='#'.OR.CH2=='!'.OR.CH2=='@'.OR.CH3=='#'.OR.CH3=='!'.OR.CH3=='@') THEN ;
              WRITE(*,'(1x,a23,a1,a4,a1,a8)') 'INPUT ERROR: Misplaced ',CH2,'and ',CH3,' symbols'
              WRITE(*,'(1x,a21,i3,a8,a120)')  '             in line ',line,'of file ',CONFILE
              WRITE(12,'(1x,a23,a1,a4,a1,a8)') 'INPUT ERROR: Misplaced ',CH2,'and ',CH3,' symbols'
              WRITE(12,'(1x,a21,i3,a8,a120)')  '             in line ',line,'of file ',CONFILE ; STOP
            ENDIF
            BACKSPACE(10) ; READ(10,*) INFILE(I)  
          ELSE
            GOTO 2
        ENDIF
      END DO
      CALL READVAR('INTR',DUMMY,OPTION(1),LINE)
      CALL READVAR('INTR',DUMMY,MODEL_TYPE,LINE)
      CALL READVAR('INTR',DUMMY,PDF_TAG,LINE)
      IF(PDF_TAG /= 0) THEN ; BACKSPACE(10) ; READ(10,*) IDUMMY,SIGMA_TAG_MOD,T_PROPORTION; ENDIF
      CALL READVAR('INTR',DUMMY,SEED,LINE)
      CALL READVAR('INTR',DUMMY,MAXITER,LINE)
      CALL READVAR('INTR',DUMMY,CHECKFLAG,LINE)
      CALL READVAR('REAL',PDEV,IDUMMY,LINE)
      CALL READVAR('INTR',DUMMY,SCALES,LINE)
      CALL READVAR('REAL',CV_OVERIDE,IDUMMY,LINE) ! over-ride index cv's to this value
      CALL READVAR('INTR',DUMMY,ADD_VAR,LINE) ! (1) additive variance (0) multiplicative variance
      CALL READVAR('INTR',DUMMY,LINK_FT,LINE) ! link selectivities in last n years
      IF(LINK_FT.GT.1) THEN ; BACKSPACE(10) ; READ(10,*) IDUMMY,SIGMA_FT,LINK_YOUNGEST,LINK_OLDEST ; ENDIF
      CALL READVAR('INTR',DUMMY,LINK_REC,LINE) ! link recruitments in last n years
      IF(LINK_REC.GT.0) THEN ; BACKSPACE(10) ; READ(10,*) IDUMMY,SIGMA_REC ; ENDIF
      CALL READVAR('INTR',DUMMY,LINK_STOCK,LINE) ! link recruitments of two stocks
      IF(LINK_STOCK.GT.0) THEN ; BACKSPACE(10) ; READ(10,*) IDUMMY,SIGMA_STOCK,RATIO_STOCK ; ENDIF
      CALL READVAR('INTR',DUMMY,PDF_STOCKRECRUIT,LINE) ! Impose a Beverton and Holt stock recruitment penalty over certain years
      IF(PDF_STOCKRECRUIT /= 0) THEN ; BACKSPACE(10) ; READ(10,*) IDUMMY,IGNORE_RECRUIT ; ENDIF
      CALL READVAR('INTR',DUMMY,OPTION(2),LINE) ! estimate F's or N's
      CALL READVAR('INTR',DUMMY,OPTION(4),LINE) ! estimate q via mle or as part of search
      CALL READVAR('INTR',DUMMY,OPTION(3),LINE) ! bootstraps
      IF(OPTION(3)>0) THEN ; BACKSPACE(10) ; READ(10,*) IDUMMY,STINE_CORR,BOOT_FILE_TYPE ; ENDIF
      CALL READVAR('INTR',DUMMY,N_RETRO,LINE) ! retrospective analyses
      N_RETRO=MIN(20,N_RETRO)
!     read covariance estimation switch, if it exists
8     READ(10,'(3A1)',END=9)  CH,CH2,CH3 ; LINE=LINE+1
      IF(CH.NE.'#'.AND.CH.NE.'!'.AND.CH.NE.'@') THEN
          IF(CH2=='#'.OR.CH2=='!'.OR.CH2=='@'.OR.CH3=='#'.OR.CH3=='!'.OR.CH3=='@') THEN ;
            WRITE(*,'(1x,a23,a1,a4,a1,a8)') 'INPUT ERROR: Misplaced ',CH2,'and ',CH3,' symbols'
            WRITE(*,'(1x,a21,i3,a8,a120)')  '             in line ',line,'of file ',CONFILE
            WRITE(12,'(1x,a23,a1,a4,a1,a8)') 'INPUT ERROR: Misplaced ',CH2,'and ',CH3,' symbols'
            WRITE(12,'(1x,a21,i3,a8,a120)')  '             in line ',line,'of file ',CONFILE ; STOP
          ENDIF
          BACKSPACE(10) ; READ(10,*) COMPUTE_COVAR
        ELSE
          IF(CH=='@') GOTO 9; GOTO 8
      ENDIF
9     IF(ABS(OPTION(3))>0 .AND. N_RETRO>0) THEN
        WRITE(*,*) 'ERROR: You cannot do retrospective analyses and bootstraps at the same time'
        WRITE(12,*) 'ERROR: You cannot do retrospective analyses and bootstraps at the same time' ; STOP
      ENDIF
      RETURN
10    WRITE(*,'(1X,A24,A120)') 'ERROR: The control file ',CONFILE
      WRITE(*,'(8X,A78)') 'does not exist or is being used by another application. Use CONTROL-C to abort'
      WRITE(*,'(8X,A26)') 'or enter a new file name: ' ; READ(*,*) CONFILE
      WRITE(12,'(1X,A24,A120)') 'ERROR: The control file ',CONFILE
      WRITE(12,'(8X,A78)') 'did not exist or was being used by another application.                       '
      WRITE(12,'(8X,A30)') 'You re-entered the file named: ' ; WRITE(12,*) CONFILE
      GOTO 1
      END

!_________________________________________________________________      
      SUBROUTINE READ_PARAMETERS  
! reads the parameters to be used in the objective function
!
!    Table of variables
!    -------------------
!    PARM_SPECS(I,J=1-5)    specifications for I'th parameter of model (1: lower bound,
!                           2: best estimate, 3: upper bound, 4: estimation indicator, 5: log-scale std. error)
!             PARM_SPECS(i,4) = 0      set to fixed value
!                               3 (.1) random deviation from previous parameter (random walk)
!                               4 (.2) random deviation from STAR_VARIABLE
!                               2 (.3) random deviation from prior "best" estimate
!                               1      independent parameter without Bayesian prior imposed
!                              -0.1    set equal to value of previous estimated parameter (0.3 or 1 designation)
!                              -n      set equal to value of parameter number n
!    PARM_EST(K)            estimate of k'th estimable parameter
!    PARM_KEY(K)            index of parameter corresponding to k'th estimable parameter
!_________________________________________________________________      
      USE PARAMETERS
      IMPLICIT NONE
      CHARACTER (LEN=1) :: CH,CH2,CH3
      REAL (KIND=8) :: DEFAULTS(5)
      INTEGER :: I,K,J,REPETITIONS,L,LINE=0,METHOD
      OPEN(11,FILE=INFILE(2),STATUS='OLD',ERR=100)
      I=1 ; K=0 ; ICALL=0   
1     READ(11,'(3A1)',END=10)  CH,CH2,CH3 ; LINE=LINE+1
!       wade through comments and read the parameter specifications
        IF(CH.EQ.'@') THEN
            GOTO 10 ! End of file
          ELSEIF(CH.NE.'#'.AND.CH.NE.'!'.AND.CH.NE.'_') THEN
            IF(CH2=='#'.OR.CH2=='!'.OR.CH2=='@'.OR.CH3=='#'.OR.CH3=='!'.OR.CH3=='@') THEN ;
              WRITE(*,'(1x,a23,a1,a4,a1,a8)') 'INPUT ERROR: Misplaced ',CH2,'and ',CH3,' symbols'
              WRITE(*,'(1x,a21,i3,a8,a120)')   '             in line ',line,'of file ',INFILE(2)
              WRITE(12,'(1x,a23,a1,a4,a1,a8)') 'INPUT ERROR: Misplaced ',CH2,'and ',CH3,' symbols'
              WRITE(12,'(1x,a21,i3,a8,a120)')   '             in line ',line,'of file ',INFILE(2) ; STOP
            ENDIF
            BACKSPACE(11)
            IF(CH.EQ.'$') THEN ! Multiple parameters are being specifed with the same form
                READ(11,*) CH,REPETITIONS,DEFAULTS
              ELSE ! One lone parameter is being specified
                READ(11,*) DEFAULTS ; REPETITIONS=1
            ENDIF
            METHOD=INT(DEFAULTS(4))
            SELECT CASE (METHOD)  ! Convert new Bayesian designation to old style
                 CASE (2) ; DEFAULTS(4)=0.3
                 CASE (3) ; DEFAULTS(4)=0.1
                 CASE (4) ; DEFAULTS(4)=0.2
                 CASE DEFAULT
            END SELECT
            IF(DEFAULTS(4)>=1.0) THEN
!             identify improperly-specified boundary conditions of free parameters
              IF(DEFAULTS(1).GT.DEFAULTS(2)) THEN
                   WRITE(*,'(1X,A59,I3)')          'ERROR: Lower bound greater than best estimate of parameter ',I
                   WRITE(*,'(1x,a25,i4,a8,a50)')   '       specified on line ',line,'of file ',INFILE(2)
                   WRITE(12,'(1X,A59,I3)')          'ERROR: Lower bound greater than best estimate of parameter ',I
                   WRITE(12,'(1x,a25,i4,a8,a50)')   '       specified on line ',line,'of file ',INFILE(2) ; STOP
                ELSEIF(DEFAULTS(3).LT.DEFAULTS(2)) THEN
                   WRITE(*,'(1X,A56,I3)')          'ERROR: Upper bound less than best estimate of parameter ',I
                   WRITE(*,'(1x,a25,i4,a8,a50)')   '       specified on line ',line,'of file ',INFILE(2)
                   WRITE(12,'(1X,A59,I3)')          'ERROR: Lower bound greater than best estimate of parameter ',I
                   WRITE(12,'(1x,a25,i4,a8,a50)')   '       specified on line ',line,'of file ',INFILE(2) ; STOP
              ENDIF
            ENDIF
            DO 2 L=1,REPETITIONS
              DO 3 J=1,5 ; PARM_SPECS(I,J)=DEFAULTS(J) ; 3 CONTINUE
!             identify estimable parameters and set them initial to the best guess
              IF(PARM_SPECS(I,4).GT.0) THEN
                K=K+1 ; PARM_KEY(K)=I
                PARM_EST(K)=PARM_SPECS(I,2)
              ENDIF
              I=I+1
2           CONTINUE
        ENDIF
      GOTO 1
10    PARM_KEY(0)=K     ! total number of parameters to be estimated 
      PARM_KEY(NPAR+1)=I-1 ! total number of parameters in model
      WRITE(*,'(1x,a31,i4)') 'Number of estimated parameters: ',PARM_KEY(0)
      WRITE(*,'(1x,a31,i4)') 'Total number of parameters:     ',PARM_KEY(NPAR+1)
      RETURN 
100   WRITE(*,'(1X,A26,A120)') 'ERROR: The parameter file ',INFILE(2)
      WRITE(*,'(8X,A55)') 'does not exist or is being used by another application.'
      WRITE(12,'(1X,A26,A120)') 'ERROR: The parameter file ',INFILE(2)
      WRITE(12,'(8X,A55)') 'does not exist or is being used by another application.'  ; PAUSE ; STOP

      END


!----------------------------------------------------------------------------        
      SUBROUTINE GETPAR(I,K,X,NPARM,VARIABLE,ETA,SIGMA)
! tells how the parameter VARIABLE should be treated (estimable, equal to previous parameter, etc...)
!----------------------------------------------------------------------------
      USE PARAMETERS, ONLY : PARM_SPECS, PARM_KEY, PI, ICALL
      USE LOGLIKELIHOODS, ONLY : PENALTY

      IMPLICIT NONE
      INTEGER I,K,J,L,NPARM
      REAL (KIND=8) :: VARIABLE,X(NPARM),ETA,SIGMA,OLD_VARIABLE,STAR_VARIABLE
      
      SIGMA=-1 ; ETA=0
      IF(ICALL==0) THEN ! Or if parm_specs(i,4)=0
          ! Equate this parameter to the fixed value specified by parm_specs(i,2)
          VARIABLE=PARM_SPECS(I,2)
        ELSEIF(PARM_SPECS(I,4) > 0.01) THEN
          ! Equate this parameter to the corresponding estimate from the search routine 
          K=K+1 ; VARIABLE=X(K)
          IF(VARIABLE>PARM_SPECS(I,3)) THEN
              PENALTY=PENALTY+10+1000*((VARIABLE-PARM_SPECS(I,3))/PARM_SPECS(I,2))**2
            ELSEIF(VARIABLE<PARM_SPECS(I,1)) THEN
              PENALTY=PENALTY+10+1000*((VARIABLE-PARM_SPECS(I,1))/PARM_SPECS(I,2))**2
          ENDIF
        ELSEIF(PARM_SPECS(I,4).LT.-.01 .AND. PARM_SPECS(I,4).GT.-1) THEN
          ! Equate this parameter to the value of the immediately preceding ESTIMABLE parameter
          VARIABLE=X(K)
        ELSEIF(PARM_SPECS(I,4)<=-1) THEN 
          ! Equate this parameter to the value of the parameter # given by parm_specs(i,4)
          J=INT(ABS(parm_specs(i,4)))
          IF(PARM_SPECS(J,4).EQ.0) THEN; VARIABLE=PARM_SPECS(J,2)
            ELSEIF(PARM_SPECS(J,4)>0 .and. J<I) THEN;
              DO L=1,NPARM; IF(PARM_KEY(L).EQ.J) EXIT ; END DO  ; VARIABLE=X(L)
            ELSE
              WRITE(*,'(1x,a33,i4,a14,i3)')  'ERROR: You cannot link parameter ',i,' to parameter ',J
              WRITE(12,'(1x,a33,i4,a14,i3)')  'ERROR: You cannot link parameter ',i,' to parameter ',J
              IF(J>=I) THEN
                  WRITE(*,'(1x,a25,i4,a50)') '       because parameter ',J,' has not been defined yet.'
                  WRITE(*,*) '       (Linked parameters must be specified after the parameter they are linked to.)'
                  WRITE(12,'(1x,a25,i4,a50)') '       because parameter ',J,' has not been defined yet.'
                  WRITE(12,*) '       (Linked parameters must be specified after the parameter they are linked to.)'
                ELSE
                  WRITE(*,*) '       because it is itself linked to another parameter.'
                  WRITE(12,*) '       because it is itself linked to another parameter.'
              ENDIF
              STOP
          ENDIF
      ENDIF
      IF(PARM_SPECS(I,4)>0.01 .AND. PARM_SPECS(I,4)<0.99) THEN ! Bayesian structures
          SIGMA = PARM_SPECS(I,5)
          IF(PARM_SPECS(I,4)<0.11) THEN ! Correlated random deviation structure (random walk)
              IF(VARIABLE/OLD_VARIABLE>0) THEN; ETA=LOG(VARIABLE/OLD_VARIABLE); ELSE; ETA=-20; ENDIF
              STAR_VARIABLE=VARIABLE ; OLD_VARIABLE=VARIABLE
            ELSEIF(PARM_SPECS(I,4)<0.21) THEN ! Uncorrelated random deviation structure
              IF(VARIABLE/STAR_VARIABLE>0) THEN; ETA=LOG(VARIABLE/STAR_VARIABLE); ELSE; ETA=-20; ENDIF
              OLD_VARIABLE=VARIABLE
            ELSEIF(PARM_SPECS(I,4)<0.31) THEN ! random deviations from prior expectation of present parameter
              IF(VARIABLE/PARM_SPECS(I,2)>0) THEN; ETA=LOG(VARIABLE/PARM_SPECS(I,2)); ELSE; ETA=-20; ENDIF
              STAR_VARIABLE=VARIABLE ; OLD_VARIABLE=VARIABLE
          ENDIF
        ELSEIF(PARM_SPECS(I,4) > 0.99) THEN ! Frequentist estimation
          STAR_VARIABLE=VARIABLE ; OLD_VARIABLE=VARIABLE
        ELSEIF(PARM_SPECS(I,4) < 0.01 .AND. PARM_SPECS(I,4) > -0.01) THEN ! fixed parameter
          STAR_VARIABLE=VARIABLE ; OLD_VARIABLE=VARIABLE
        ELSE
          STAR_VARIABLE=VARIABLE ; OLD_VARIABLE=VARIABLE
      ENDIF
      I=I+1
      RETURN ; END

!----------------------------------------------------------------------------        
      SUBROUTINE WRITEPAR(I,K)
! Writes the values of the parameters to a file
!
! Table of variables
!     PARM_SPECS(i,4) = 0      set to fixed value
!                       1      independent parameter without Bayesian prior imposed
!                       2 (.3) random deviation from prior "best" estimate
!                       3 (.1) random deviation from previous parameter
!                       4 (.2) random deviation from reference parameter
!                      -0.1    set equal to value of previous estimated parameter (0.3 or 1 designation)
!                      -n      set equal to value of parameter number n
!
!----------------------------------------------------------------------------
      USE PARAMETERS ; USE STATISTICS, ONLY : FIRSTYEAR, LASTYEAR, DISPLAYYEAR, FIRSTAGE, LASTAGE
      IMPLICIT NONE
      INTEGER :: I,K,J,L,WRITE_YRAGE=0,A,Y
      REAL (KIND=8) :: VARIABLE,CV,VARIANCE,METHOD
      CHARACTER (LEN=5) :: BOUND
      IF(I>=0) THEN
          WRITE_YRAGE=0
        ELSE
          IF(WRITE_YRAGE==0) THEN ; Y=FIRSTYEAR+DISPLAYYEAR-1 ; A=FIRSTAGE ; WRITE_YRAGE=1 ; ENDIF
          I=ABS(I)
      ENDIF
      CV=-1
      SELECT CASE (INT(PARM_SPECS(I,4)*10))
          CASE (3)     ; METHOD=2
          CASE (1)     ; METHOD=3
          CASE (2)     ; METHOD=4
          CASE DEFAULT ; METHOD=PARM_SPECS(I,4)
      END SELECT
      IF(PARM_SPECS(I,4)==0) THEN ! parameter set to fixed value
          VARIABLE=PARM_SPECS(I,2) 
        ELSEIF(PARM_SPECS(I,4)>0) THEN 
          K=K+1 ; VARIABLE=PARM_EST(K)
          VARIANCE=PARM_VAR(K)
          IF(VARIANCE>0) CV=100*ABS(SQRT(VARIANCE)/VARIABLE)
        ELSEIF(PARM_SPECS(I,4).LT.0 .AND. PARM_SPECS(I,4).GT.-1) THEN ! Equate to most recently estimated parameter value
          VARIABLE=PARM_EST(K)
        ELSEIF(PARM_SPECS(I,4).LE.-1) THEN ! Equate to value of parameter n
          J=INT(ABS(parm_specs(i,4)))
          IF(PARM_SPECS(J,4).EQ.0) THEN; VARIABLE=PARM_SPECS(J,2) 
            ELSEIF(PARM_SPECS(J,4)>0) THEN;
              DO L=1,PARM_KEY(0); IF(PARM_KEY(L).EQ.J) EXIT ; END DO ; VARIABLE=PARM_EST(L)
          ENDIF
      ENDIF
      IF((VARIABLE<=1.01*PARM_SPECS(I,1) .OR. VARIABLE>=0.99*PARM_SPECS(I,3)) .AND. PARM_SPECS(I,4)>=1) THEN 
        BOUND='BOUND' ; ELSE ; BOUND=' '
      ENDIF
      IF(WRITE_YRAGE==0) THEN
          IF(CV<=0) THEN ; WRITE(50,100) PARM_SPECS(I,1),VARIABLE,PARM_SPECS(I,3),METHOD,PARM_SPECS(I,5),I,BOUND
            ELSE ; WRITE(50,101) PARM_SPECS(I,1),VARIABLE,PARM_SPECS(I,3),METHOD,PARM_SPECS(I,5),I,k,CV,BOUND
          ENDIF
        ELSE
          IF(CV<=0) THEN ; WRITE(50,102) PARM_SPECS(I,1),VARIABLE,PARM_SPECS(I,3),METHOD,PARM_SPECS(I,5),I,Y,A,BOUND
            ELSE ; WRITE(50,103) PARM_SPECS(I,1),VARIABLE,PARM_SPECS(I,3),METHOD,PARM_SPECS(I,5),I,Y,A,k,CV,BOUND
          ENDIF
          IF(A<LASTAGE) THEN ; A=A+1 ; ELSE ; A=FIRSTAGE ; Y=Y+1 ; ENDIF
      ENDIF
      I=I+1
100   FORMAT(1X,3(D11.4,2X),F7.1,2X,D11.4,2X,I4,14X,A5)
101   FORMAT(1X,3(D11.4,2X),F7.1,2X,D11.4,2X,I4,1x,i4,2X,F5.0,2x,A5)
102   FORMAT(1X,3(D11.4,2X),F7.1,2X,D11.4,2X,3(I4,1X),13X,A5)
103   FORMAT(1X,3(D11.4,2X),F7.1,2X,D11.4,2X,4(I4,1X),1X,F5.0,2x,A5)
      RETURN ; END

!------------------------------------------------------------------------------
      SUBROUTINE ESTIMATE(DO_PRINT)
! Calls the search routine for estimating the parameters      
!------------------------------------------------------------------------------
      USE PARAMETERS, ONLY : SEED, MAXITER, OPTION, CHECKFLAG, PARM_KEY, PDEV, NPAR
      INTEGER ITER,DO_PRINT
      REAL (KIND=8) :: DUMMY1(501,500),DUMMY2(501),DUMMY3(501),DUMMY4(500)
      IF(PARM_KEY(0).GT.NPAR) THEN
        WRITE(*,999) PARM_KEY(0) 
        WRITE(12,999) PARM_KEY(0)
        WRITE(*,*) '        EITHER REFORMULATE THE MODEL OR INCREASE THE DIMENSIONS IN THE DUMMY ARRAYS'
        WRITE(12,*) '        EITHER REFORMULATE THE MODEL OR INCREASE THE DIMENSIONS IN THE DUMMY ARRAYS' ; STOP
      ENDIF
      ITER=MAXITER

! {conduct search for best estimates of parameters}      
      CALL SEARCH(DUMMY1,DUMMY2,DUMMY3,DUMMY4,PARM_KEY(0),SEED,ITER,CHECKFLAG,PDEV,DO_PRINT)
! {notify user of poor performance}      
      IF(ITER.EQ.-1) THEN    
          WRITE(*,*) 'MAXIMUM NUMBER (',MAXITER,') OF ITERATIONS EXCEEDED'
          WRITE(12,*) 'MAXIMUM NUMBER (',MAXITER,') OF ITERATIONS EXCEEDED'
      ENDIF
999   FORMAT('ERROR:  You are trying to estimate too many (',I4,') parameters')
      RETURN ; END

!------------------------------------------------------------------------------
      SUBROUTINE PROFILES
! Generates likelihood profiles for selected parameters      
!------------------------------------------------------------------------------
      USE PARAMETERS
      IMPLICIT NONE
      INTEGER I,J,K,KK,IT,N,CHECK,DO_PRINT
      PARAMETER(IT=20, CHECK=2, DO_PRINT=0)
      REAL (KIND=8) :: DUMMY1(501,500),DUMMY2(501),DUMMY3(501),DUMMY4(500),PD,PROFILE(0:500),OLD_4, OLD_2,OLD_3,OLDPAR(NPAR)
      PD=0.5
      OPEN(13,FILE=INFILE(8),STATUS='UNKNOWN')
      WRITE(*,*) 'WRITING TO FILE ',INFILE(8)
      OLDPAR=PARM_EST  ; OLD_3=PARM_KEY(0) ; K=0
      DO 1 I=1,PARM_KEY(NPAR+1)
        IF(PARM_SPECS(I,4).GT.0) K=K+1
        IF(ABS(PARM_SPECS(I,5)).GT.0) THEN
          WRITE(13,*) '  PARAMETER                VALUE              LIKELIHOOD'
          N=INT(ABS(PARM_SPECS(I,5)))-1  
          OLD_4=PARM_SPECS(I,4) ; OLD_2=PARM_SPECS(I,2) 
          IF(PARM_SPECS(I,4).GT.0) THEN
            PARM_SPECS(I,4)=0 ! Fix new value of parameter to redefined best estimate (above)
            DO KK=K,PARM_KEY(0)-1
              PARM_EST(KK)=PARM_EST(KK+1) ; PARM_KEY(KK)=PARM_KEY(KK+1)
            END DO
            PARM_KEY(0)=PARM_KEY(0)-1
          ENDIF
          DO 2 J=0,N
            ICALL=0 
            IF(PARM_SPECS(I,5).GT.0) THEN 
                PARM_SPECS(I,2)=PARM_SPECS(I,1)+(PARM_SPECS(I,3)-PARM_SPECS(I,1))*DBLE(J)/DBLE(N) 
              ELSE
                PARM_SPECS(I,2)=DEXP(DLOG(PARM_SPECS(I,1))+(DLOG(PARM_SPECS(I,3))-DLOG(PARM_SPECS(I,1)))*DBLE(J)/DBLE(N)) 
            ENDIF
            CALL SEARCH(DUMMY1,DUMMY2,DUMMY3,DUMMY4,PARM_KEY(0),SEED,IT,CHECK,PD,DO_PRINT) ! Conduct search to determine new likelihood
            PROFILE(J)=PARM_EST(PARM_KEY(0)+1)
            WRITE(*,*) I,PARM_SPECS(I,2),parm_key(0),PROFILE(J)
2           WRITE(13,*) I,PARM_SPECS(I,2),PROFILE(J)
          PARM_EST=OLDPAR ; PARM_SPECS(I,4)=OLD_4 ; PARM_SPECS(I,2)=OLD_2 ; PARM_KEY(0)=OLD_3
        ENDIF
1     CONTINUE
      CLOSE(13)
      RETURN
      END

!------------------------------------------------------------------------------
      SUBROUTINE FORMULAE(OBS,PRED,SIGMATEMP,SCALE_VAR,CHOICE,ANSWER,CONSTANTS)
! computes log-likelihoods of different distributions (not negative log-likelihood)
! CAUTION: Be careful not to modify any variable other than answer
!
! SIGMATEMP = coefficient of variation if positive and standard deviation if negative.
!             Depending on the pdf, sigmatemp is converted to either the cv or the std. dev.
!             by use of the observed data; i.e., cv= |std. dev. / observed value| or
!             std. dev. = cv*obs      
!------------------------------------------------------------------------------
      USE PARAMETERS, ONLY : ADD_VAR
      REAL (KIND=8) :: OBS,PRED,SIGMATEMP,SIGMA,SCALE_VAR,SCALE_VAR2,ANSWER,HALFLOG2PI,LOGFACTORIAL,A,B,LOGSQRT2
      INTEGER CHOICE
      CHARACTER*1 CONSTANTS
      PARAMETER(HALFLOG2PI=0.91893853D0, LOGSQRT2=0.34657359D0)
      SIGMA=SIGMATEMP ; SCALE_VAR2=SCALE_VAR*SCALE_VAR; IF(ADD_VAR<=0 .AND. SCALE_VAR2<=0) SCALE_VAR2=1.0;
      IF(SIGMA.EQ.0) SIGMA=0.1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The special case missing data (neg. values)
! Note that negative PRED values should be penalized in subroutine LIKELIHOOD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(OBS.LT.0 .AND. CHOICE /= 12) THEN ; ANSWER=0.D0 ; RETURN ; ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute likelihoods of CHOICE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SELECT CASE (CHOICE)

        CASE (0) ! ignore
          RETURN
        CASE (1) ! lognormal
          IF(OBS.EQ.0) THEN
            WRITE(*,*) ' ERROR: The lognormal density was selected for the likelihood, but one of the'
            WRITE(*,*) '        observations is zero-- LOG(0) is undefined'
            WRITE(12,*) ' ERROR: The lognormal density was selected for the likelihood, but one of the'
            WRITE(12,*) '        observations is zero-- LOG(0) is undefined' ; PAUSE ; STOP
          ENDIF
          IF(SIGMA<0) THEN  ! convert std. dev. to CV
            IF(PRED>0) THEN ; SIGMA = ABS(SIGMA/PRED) ; ELSE ; SIGMA = ABS(SIGMA/OBS) ; ENDIF
          ENDIF
          SIGMA = DLOG(SIGMA**2.0 + 1.0d0) ! Convert CV to log-scale variance
          IF(ADD_VAR>0) THEN ; SIGMA=SQRT(SIGMA+SCALE_VAR2) ; ELSE ; SIGMA=SQRT(SIGMA*SCALE_VAR2) ; ENDIF
          IF(PRED>0) THEN
              ANSWER = -0.5d0*((DLOG(OBS)-DLOG(PRED))/SIGMA)**2.0 -DLOG(SIGMA)
            ELSE
              ANSWER = -0.5d0*((DLOG(OBS)+14.0)/SIGMA)**2.0 -DLOG(SIGMA)
          ENDIF
          IF(CONSTANTS.EQ.'Y') ANSWER=ANSWER-HALFLOG2PI-DLOG(OBS)

        CASE (2) ! normal
          IF(SIGMA>0) THEN  ! {convert cv's to variances}
              IF(PRED>0) THEN
                  SIGMA = (SIGMA*PRED)**2.0
                ELSEIF(OBS>0) THEN
                  SIGMA = (SIGMA*OBS)**2.0
                ELSE
                  SIGMA = 1.0 ! hopefully gives moderate weight to a zero obs. (treats it as though it were 1.0 with a cv of 1.0)
              ENDIF
          ENDIF
          SIGMA=SIGMA*SIGMA ! variance
          IF(ADD_VAR>0) THEN ; SIGMA=SQRT(SIGMA+SCALE_VAR2) ; ELSE ; SIGMA=SQRT(SIGMA*SCALE_VAR2) ; ENDIF
          ANSWER = -0.5d0 * ((OBS - PRED)/SIGMA)**2.0 - DLOG(ABS(SIGMA))
          IF(CONSTANTS.EQ.'Y') ANSWER=ANSWER-HALFLOG2PI

        CASE (3) ! multinomial
	  !     {obs = sample size * observed proportion, pred = estimated proportion}
          IF(OBS.GE.0 .AND. (PRED.GT.0 .AND. PRED.LE.1.D0)) THEN
              ANSWER = OBS*DLOG(PRED)
            ELSEIF(PRED.LE.0) THEN
!             {zero predictions are feasible, but not in the log-likelihood formula
              ANSWER = -OBS*74.0 ! log(1.d-32)=74
            ELSE
!             {predictions > 1 are infeasible}
              ANSWER = -OBS*74.0 ! log(1.d-32)=74
          ENDIF
          IF(CONSTANTS.EQ.'Y') ANSWER=ANSWER-LOGFACTORIAL(OBS)

        CASE (4) ! Poisson
          IF(OBS.GE.0 .AND. PRED.GT.0) THEN
              ANSWER = -PRED + OBS*DLOG(PRED)
            ELSEIF(PRED.EQ.0) THEN
!             {zero predictions are feasible, but not in loglikelihood formula}
              ANSWER = OBS*DLOG(1.D-06)
          ENDIF
          IF(CONSTANTS.EQ.'Y') ANSWER=ANSWER-LOGFACTORIAL(OBS)

        CASE (5) ! Chi-square criterion (robustified)
          IF(PRED.GE.0) THEN
              ANSWER = -(((OBS - PRED)/SIGMA)**2.0)/(PRED + 1.d0) ! sigma here is just a means of downweighting, not true variance
            ELSE
!             {negative predictions are infeasible}
              ANSWER = -2.0*(OBS - PRED)**2.0
          ENDIF

        CASE (6) ! Robust (double exponential)
          IF(SIGMA>0) THEN  ! {convert cv's to std. dev.'s}
              IF(PRED>0) THEN
                  SIGMA = ABS(SIGMA*PRED)
                ELSEIF(OBS>0) THEN
                  SIGMA = ABS(SIGMA*OBS)
                ELSE
                  SIGMA = 1.0 ! hopefully gives moderate weight to a zero obs. (treats it as though it were 1.0 with a cv of 1.0)
              ENDIF
          ENDIF
          SIGMA=SIGMA*SIGMA ! variance
          IF(ADD_VAR>0) THEN ; SIGMA=SQRT(SIGMA+SCALE_VAR2) ; ELSE ; SIGMA=SQRT(SIGMA*SCALE_VAR2) ; ENDIF
          ANSWER = -LOGSQRT2*ABS((OBS - PRED)/SIGMA) - LOG(ABS(SIGMA))
          IF(CONSTANTS.EQ.'Y') ANSWER=ANSWER-LOGSQRT2

        CASE (7) ! Gamma (using Stirling's approximation for the gamma function)
          IF(OBS.EQ.0) THEN
            WRITE(*,*) ' ERROR: The gamma density was selected, but one or more of the'
            WRITE(*,*) '        observations are zero-- LOG(0) is undefined'
            WRITE(12,*) ' ERROR: The gamma density was selected, but one or more of the'
            WRITE(12,*) '        observations are zero-- LOG(0) is undefined' ; STOP
          ENDIF
          IF(SIGMA>0) THEN  ! {convert cv's to std. dev.'s}
              IF(PRED>0) THEN
                  SIGMA = ABS(SIGMA*PRED)
                ELSEIF(OBS>0) THEN
                  SIGMA = ABS(SIGMA*OBS)
                ELSE
                  SIGMA = 1.0 ! hopefully gives moderate weight to a zero obs. (treats it as though it were 1.0 with a cv of 1.0)
              ENDIF
          ENDIF
          SIGMA=SIGMA*SIGMA ! variance
          IF(ADD_VAR>0) THEN ; SIGMA=SIGMA+SCALE_VAR2 ; ELSE ; SIGMA=SIGMA*SCALE_VAR2 ; ENDIF
          A=PRED*PRED/SIGMA ; B=SIGMA/PRED ! get a and b from variances
          ANSWER = -A*LOG(B) + (A-1.)*LOG(OBS) - OBS/B - LOGFACTORIAL(a-1.0)

        CASE (12) ! Standard normal deviates (index - index mean)
          ANSWER = -0.5d0 * ((OBS - PRED)/SIGMA)**2.0 - DLOG(ABS(SIGMA))
          IF(CONSTANTS.EQ.'Y') ANSWER=ANSWER-HALFLOG2PI

        CASE DEFAULT
          WRITE(*,'(1X,A9,I20,A43)') 'PDF TYPE ',CHOICE,' NOT AVAILABLE FOR LIKELIHOOD COMPUTATIONS' ; STOP

      END SELECT

      RETURN ; END
!--------------------------------------------------------------------------
       SUBROUTINE GENDATA(xbar,SIGMATEMP,CHOICE,ANSWER,ISEED)   
!  This subroutine generates a random variate that follows one of
!  several possible types of distributions:
!      1=LOGNORMAL, 2=NORMAL, 3=MULTINOMIAL (NOT AVAILABLE), 4=POISSON
!      5=CHI-SQUARE (NOT AVAILABLE), 6=DOUBLE-EXPONENTIAL (NOT AVAILABLE)
!      7=UNIFORM, 8=TRIANGULAR CENTERED AT XBAR.
!  Note: all inputs must be on arithmetic scale
!--------------------------------------------------------------------------

      implicit none
      integer CHOICE,iseed,ib,i,ICHOICE
      REAL (KIND=8) ::  a,b,xmu,c,bot,r,xbar,answer,RNORM,RAN1,SIGMATEMP,sigma
 
      external rnorm,ran1
      SIGMA=SIGMATEMP ; ICHOICE=CHOICE
      IF(SIGMA.EQ.0) THEN ; ANSWER=XBAR ; RETURN ; ENDIF
      IF(CHOICE==4 .AND. XBAR>80) THEN ; ICHOICE=2 ; SIGMA=-SQRT(XBAR) ; ENDIF ! {use normal approximation to Poisson for lg XBAR}

      SELECT CASE (ICHOICE)

        CASE (1) ! lognormal
          IF(XBAR.LE.0) THEN
            WRITE(*,*) ' ERROR: The lognormal density was selected for the likelihood, but one of the'
            WRITE(*,*) '        observations is zero-- LOG(0) is undefined'
            WRITE(12,*) ' ERROR: The lognormal density was selected for the likelihood, but one of the'
            WRITE(12,*) '        observations is zero-- LOG(0) is undefined' ; STOP
          ENDIF
          IF(SIGMA.LT.0) SIGMA = ABS(SIGMA/XBAR) ! convert std. dev. to CV
          SIGMA = LOG(SIGMA**2 + 1.0d0) ; xmu=LOG(XBAR)-SIGMA/2
          ANSWER=DEXP(XMU+SQRT(SIGMA)*RNORM(iseed))

        CASE (2) ! normal
          IF(SIGMA.GT.0) THEN ; SIGMA = ABS(SIGMA*XBAR) ; ELSEIF(SIGMA.LT.0) THEN ; SIGMA = ABS(SIGMA) ; ENDIF
          ANSWER=XBAR+SIGMA*RNORM(iseed)

        CASE (3) ! multinomial
          WRITE(*,*) 'ERROR: Multinomial density not avaialble for random variable generation'
          WRITE(12,*) 'ERROR: Multinomial density not avaialble for random variable generation' ; STOP

        CASE (4) ! Poisson
          IF(XBAR<0) THEN
              WRITE(*,*) 'ERROR: The Poisson density was selected, but the input expectation is < 0'
              WRITE(12,*) 'ERROR: The Poisson density was selected, but the input expectation is < 0' ; RETURN
            ELSE
              a=EXP(-XBAR) ; b=1 ; ib=0
              do 20 i=1,1000
                b=b*RAN1(iseed)
                IF(B.LT.A) THEN ; ANSWER=IB ; RETURN ; ELSE ; IB=IB+1 ; ENDIF
20            CONTINUE
          ENDIF

        CASE (5) ! chi-square
          WRITE(*,*) 'CHI-SQUARE density not avaialble for random variable generation'
          WRITE(12,*) 'CHI-SQUARE density not avaialble for random variable generation' ; STOP

        CASE (6) ! double-exponential
          WRITE(*,*) 'DOUBLE EXPONENTIAL DENSITY density not avaialble for random variable generation'
          WRITE(12,*) 'DOUBLE EXPONENTIAL DENSITY density not avaialble for random variable generation' ; STOP

        CASE (7) ! uniform:
          IF(SIGMA<0) SIGMA=ABS(SIGMA/XBAR)
          A=XBAR*(1-SIGMA*DSQRT(3.D0)) ; B=2.*XBAR - A ; ANSWER=A+(B-A)*RAN1(iseed)

        CASE (8) ! triangular centered at XBAR (i.e., C parameter = (a+b)/2)
          IF(SIGMA<0) SIGMA=ABS(SIGMA/XBAR)
          b=XBAR*(1.+DSQRT(6.D0)*SIGMA) ; a=2.*XBAR-b ; bot=b-a ; c=xbar ; r=RAN1(iseed)
          IF(r.le.c) THEN ; answer=a+bot*SQRT(r*c) ; ELSE ; answer=a+bot*(1.-SQRT((1.-c)*(1.-r))) ; ENDIF

        CASE DEFAULT ! no random deviations
          ANSWER=XBAR

      END SELECT

      RETURN ; END

!------------------------------------------------------------------------------
      SUBROUTINE GETVARIANCE(OBS,PRED,SIGMATEMP,SCALE_VAR,CHOICE,SIGMA)
! computes variances of various data types
! CAUTION: Be careful not to modify any variable other than answer
!------------------------------------------------------------------------------
      USE PARAMETERS, ONLY : ADD_VAR
      REAL (KIND=8) :: OBS,PRED,SIGMATEMP,SIGMA,SCALE_VAR,SCALE_VAR2
      INTEGER CHOICE
      SIGMA=SIGMATEMP ; SCALE_VAR2=SCALE_VAR*SCALE_VAR
      IF(SIGMA.EQ.0) SIGMA=0.1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The special cases of missing data (neg. values)
! Note that negative PRED values should be penalized in subroutine LIKELIHOOD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(OBS.LT.0 .AND. CHOICE /=12) THEN ; WRITE(*,*) 'OBS<0' ; SIGMA=0.D0 ; RETURN ; ENDIF

      SELECT CASE (CHOICE)

        CASE (0) ! IGNORE
          RETURN
        CASE (1) ! lognormal
          IF(SIGMA<0) THEN  ! convert std. dev. to CV
            IF(PRED>0) THEN ; SIGMA = ABS(SIGMA/PRED) ; ELSE ; SIGMA = ABS(SIGMA/OBS) ; ENDIF
          ENDIF
          SIGMA = DLOG(SIGMA**2.0 + 1.0d0) ! Convert CV to log-scale variance
          IF(ADD_VAR>0) THEN ; SIGMA=SIGMA+SCALE_VAR2 ; ELSE ; SIGMA=SIGMA*SCALE_VAR2 ; ENDIF

        CASE (2,6,7) ! normal, laplace and gamma
          IF(SIGMA>0) THEN  ! {convert cv's to std. dev.'s}
              IF(PRED>0) THEN
                  SIGMA = ABS(SIGMA*PRED)
                ELSEIF(OBS>0) THEN
                  SIGMA = ABS(SIGMA*OBS)
                ELSE
                  SIGMA = 1.0 ! hopefully gives moderate weight to a zero obs. (treats it as though it were 1.0 with a cv of 1.0)
              ENDIF
          ENDIF
          SIGMA=SIGMA*SIGMA ! variance
          IF(ADD_VAR>0) THEN ; SIGMA=SIGMA+SCALE_VAR2 ; ELSE ; SIGMA=SIGMA*SCALE_VAR2 ; ENDIF

        CASE (3) ! multinomial
	  !     {obs = sample size * observed proportion, pred = estimated proportion}

        CASE (4) ! Poisson
          IF(PRED>0) THEN ; SIGMA=PRED ; ELSE ; SIGMA=OBS ; ENDIF

        CASE (5) ! Chi-square criterion (robustified)
          ! NA

        CASE (12) ! Standard normal deviates (index - index mean)
          SIGMA=SIGMA**2

        CASE DEFAULT
          WRITE(*,'(1X,A9,I3,A41)') 'PDF TYPE ',CHOICE,' NOT AVAILABLE FOR VARIANCE COMPUTATIONS' ; STOP

      END SELECT

      RETURN ; END

!------------------------------------------------------------------------------
      SUBROUTINE GETDISCREPANCY(OBS,PRED,SIGMA2,CHOICE,DISCREPANCY)
! computes Chi-squared discrepancy statistic
! CAUTION: Be careful not to modify any variable other than answer
!------------------------------------------------------------------------------
      implicit none
      REAL (KIND=8) :: OBS,PRED,SIGMA2,DISCREPANCY,EXPECTATION,VARIANCE
      INTEGER CHOICE

      IF(CHOICE==1) THEN ; ! Convert log-scale variances to arithmetic scale and median to mean
          EXPECTATION=PRED*EXP(SIGMA2/2.0)
          IF(PRED>0) THEN ; VARIANCE = (EXP(SIGMA2)-1)*EXPECTATION*EXPECTATION ; ELSE ; VARIANCE=(EXP(SIGMA2)-1.0)*OBS*OBS ; ENDIF
        ELSE
          EXPECTATION=PRED; VARIANCE=SIGMA2
      ENDIF
      DISCREPANCY=(OBS-EXPECTATION)*(OBS-EXPECTATION)/VARIANCE
      RETURN ; END

!_________________________________________________________________
      SUBROUTINE BLANK(IVAR,RECORD)
!_________________________________________________________________
      IMPLICIT NONE ; INTEGER IVAR,RECORD ; CHARACTER CH*1, CH2*1, CH3*1

1     READ(IVAR,'(3A1)')  CH,CH2,CH3
      RECORD=RECORD+1
      IF(CH.NE.'#'.AND.CH.NE.'!'.AND.CH.NE.'_'.AND.CH.NE.'@') THEN
        IF(CH2=='#'.OR.CH2=='!'.OR.CH2=='_'.OR.CH3=='#'.OR.CH3=='!'.OR.CH3=='_') THEN ;
          WRITE(*,*) 'INPUT ERROR: Misplaced #, ! or _ symbol in line ',record
          WRITE(12,*) 'INPUT ERROR: Misplaced #, ! or _ symbol in line ',record
         ELSE
          BACKSPACE(IVAR) ; RETURN
        ENDIF
       ELSE
          GOTO 1
      ENDIF
      END

!_________________________________________________________________
      SUBROUTINE READVAR(TYPE,VAR,IVAR,LINE)
!_________________________________________________________________
      IMPLICIT NONE
      REAL (KIND=8) :: VAR
      INTEGER :: IVAR,LINE
      CHARACTER (LEN=1) :: CH,CH2,CH3
      CHARACTER (LEN=4) :: TYPE

1     READ(10,'(3A1)')  CH,CH2,CH3 ; LINE=LINE+1
        IF(CH.NE.'#'.AND.CH.NE.'!'.AND.CH.NE.'_'.AND.CH.NE.'@') THEN
            IF(CH2=='#'.OR.CH2=='!'.OR.CH2=='@'.OR.CH3=='#'.OR.CH3=='!'.OR.CH3=='@') THEN ;
              WRITE(*,'(1x,a23,a1,a4,a1,a8)') 'INPUT ERROR: Misplaced ',CH2,'and ',CH3,' symbols'
              WRITE(*,'(1x,a21,i3,a8,a20)')   '             in line ',line,'of the control file.'
              WRITE(12,'(1x,a23,a1,a4,a1,a8)') 'INPUT ERROR: Misplaced ',CH2,'and ',CH3,' symbols'
              WRITE(12,'(1x,a21,i3,a8,a20)')   '             in line ',line,'of the control file.' ; STOP
            ENDIF
            BACKSPACE(10)
            IF(TYPE.EQ.'INTR') THEN
                READ(10,*) IVAR
              ELSEIF(TYPE.EQ.'REAL') THEN
                READ(10,*) VAR
            ENDIF
            RETURN
          ELSEIF(CH.NE.'@') THEN
            GOTO 1
          ELSE
            WRITE(*,*) 'INPUT ERROR: improperly placed @'
            WRITE(12,*) 'INPUT ERROR: improperly placed @' ; STOP
        ENDIF
      END

!_________________________________________________________________
      SUBROUTINE WRITE_REAL(IFILE,WIDTH,VALUE,ADV)
!_________________________________________________________________
      IMPLICIT NONE
      REAL (KIND=8) :: VALUE
      INTEGER :: IFILE, WIDTH
      CHARACTER (LEN=2) :: ADV
      CHARACTER (LEN=2) FFORM(7:14), DFORM(7:14)
      CHARACTER (LEN=2) TEMPCH(6)
      DATA FFORM/'7' ,'8' ,'9' ,'10','11','12','13','14'/
      DATA DFORM/'.2','.2','.3','.3','.4','.5','.5','.5'/
      DATA TEMPCH/'(','','','','',')'/
      TEMPCH(2)=FFORM(WIDTH); TEMPCH(4)='' ;TEMPCH(5)=''
      IF(ABS(VALUE)<10.0**(WIDTH-2)) THEN
          TEMPCH(1)='(F' ; TEMPCH(3)='.2'
        ELSE
          TEMPCH(1)='(E' ; TEMPCH(3)=DFORM(WIDTH);
      ENDIF
      IF(ADV/='EN') THEN
          IF(ADV/='0') THEN ; TEMPCH(4)=',' ; TEMPCH(5)=ADV ; ENDIF
          WRITE(IFILE,TEMPCH,ADVANCE='NO') VALUE
        ELSE
          WRITE(IFILE,TEMPCH,ADVANCE='YES') VALUE
      ENDIF
      RETURN ; END

!_________________________________________________________________
      SUBROUTINE WRITE_CHAR(IFILE,WIDTH,CH,ADV)
!_________________________________________________________________
      IMPLICIT NONE
      INTEGER :: IFILE, WIDTH
      CHARACTER (LEN=WIDTH) :: CH
      CHARACTER (LEN=2) :: ADV
      CHARACTER (LEN=2) AFORM(1:30)
      CHARACTER (LEN=2) TEMPCH(5)
      DATA AFORM/'1' ,'2' ,'3' ,'4' ,'5' ,'6' ,'7' ,'8' ,'9' ,'10', &
                 '11','12','13','14','15','16','17','18','19','20', &
                 '21','22','23','24','25','26','27','28','29','30'/
      DATA TEMPCH/'(A','1','','',')'/
      TEMPCH(2)=AFORM(WIDTH); TEMPCH(3)='' ; TEMPCH(4)=''
      IF(ADV/='EN') THEN
          IF(ADV/='0') THEN ; TEMPCH(3)=',' ; TEMPCH(4)=ADV ; ENDIF
          WRITE(IFILE,TEMPCH,ADVANCE='NO') CH
        ELSE
          WRITE(IFILE,TEMPCH,ADVANCE='YES') CH
      ENDIF
      RETURN ; END

!_________________________________________________________________
      SUBROUTINE WRITE_INT(IFILE,WIDTH,INTGR,ADV)
!_________________________________________________________________
      IMPLICIT NONE
      INTEGER :: IFILE, WIDTH,INTGR
      CHARACTER (LEN=2) :: ADV
      CHARACTER (LEN=2) IFORM(1:30)
      CHARACTER (LEN=2) TEMPCH(5)
      DATA IFORM/'1' ,'2' ,'3' ,'4' ,'5' ,'6' ,'7' ,'8' ,'9' ,'10', &
                 '11','12','13','14','15','16','17','18','19','20', &
                 '21','22','23','24','25','26','27','28','29','30'/
      DATA TEMPCH/'(I','1','','',')'/
      TEMPCH(2)=IFORM(WIDTH); TEMPCH(3)='' ; TEMPCH(4)=''
      IF(ADV/='EN') THEN
          IF(ADV/='0') THEN ; TEMPCH(3)=',' ; TEMPCH(4)=ADV ; ENDIF
          WRITE(IFILE,TEMPCH,ADVANCE='NO') INTGR
        ELSE
          WRITE(IFILE,TEMPCH,ADVANCE='YES') INTGR
      ENDIF
      RETURN ; END

!_________________________________________________________________
      SUBROUTINE CHECKDIMENSION(I,IDIM,IMAX,CH,RECORD)
!_________________________________________________________________
      IMPLICIT NONE
      INTEGER I,IMAX, IDIM , RECORD
      CHARACTER*10 CH
      IF(I.GT.IDIM) THEN
        WRITE(*,'(1X,A21,A10,A1,I3)') 'ERROR: The index for ',CH,' ',I
        WRITE(*,*) '       exceeds the maximum dimension allowed by the program (',IDIM,')'
        WRITE(12,'(1X,A21,A10,A1,I3)') 'ERROR: The index for ',CH,' ',I
        WRITE(12,*) '       exceeds the maximum dimension allowed by the program (',IDIM,')'
       ELSEIF(I.GT.IMAX) THEN
        WRITE(*,'(1X,A21,A10,A1,I3)') 'ERROR: The index for ',CH,' ',I
        WRITE(*,*) '       exceeds the maximum number indicated in the data file(',IMAX,')'
        WRITE(12,'(1X,A21,A10,A1,I3)') 'ERROR: The index for ',CH,' ',I
        WRITE(12,*) '       exceeds the maximum number indicated in the data file(',IMAX,')'
       ELSEIF(I.LT.0) THEN
        WRITE(*,'(1X,A21,A10,A1,I3)') 'ERROR: The index for ',CH,' ',I
        WRITE(*,*) '       is less than zero'
        WRITE(12,'(1X,A21,A10,A1,I3)') 'ERROR: The index for ',CH,' ',I
        WRITE(12,*) '       is less than zero'
       ELSE
        GOTO 1
      ENDIF
      WRITE(*,'(1X,A9,I5)') 'Record # ',RECORD
      WRITE(12,'(1X,A9,I5)') 'Record # ',RECORD ; STOP
1     RETURN ; END

!------------------------------------------------------------------------------
      REAL (KIND=8) FUNCTION LOGFACTORIAL(OBS)
!     Returns (obs)! or, equivalently, GAMMA(obs+1)
!     uses STIRLINGs approximation to GAMMA(obs+1) for non-integer arguments
!        exp(-a)*(a**(a-.5))*2.506628*(1.+1./(12.*a)+1./(288.*a*a)-139./(51840.*a**3.)-571./(2488320.*a**4.))
!------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NUM,I
      REAL (KIND=8) :: OBS,REMAINDER,A
      NUM=NINT(OBS); REMAINDER=OBS-DBLE(NUM)
      IF(REMAINDER<1.0D-16) THEN
          ! factorial of integer may be computed exactly
     	  LOGFACTORIAL=0.D0
          DO I=NUM,2,-1 ; LOGFACTORIAL=LOGFACTORIAL+DLOG(DBLE(I)); END DO
        ELSE
          A=OBS+1.0
          LOGFACTORIAL= -a+(a-.5)*LOG(a)+0.91893853+LOG(1.+1./(12.*a)+1./(288.*a*a)-139./(51840.*a**3.)-571./(2488320.*a**4.))
      ENDIF
      RETURN
      END


!_________________________________________________________________
      REAL (KIND=8) FUNCTION RNORM(IDUM)
!_________________________________________________________________
! generates normally distributed random deviates N(0,1)
      REAL (KIND=8) RAN1
      EXTERNAL RAN1
      
      integer iset,idum
      real (KIND=8) v1,v2,r,fac,gset

      DATA ISET/0/
      IF (ISET.EQ.0) THEN
1       V1=2.d0*RAN1(IDUM)-1.d0
        V2=2.d0*RAN1(IDUM)-1.d0
        R=V1*V1+V2*V2
        IF(R.GE.1.)GO TO 1
        FAC=DSQRT(-2.D0*DLOG(R)/R)
        GSET=V1*FAC
        RNORM=V2*FAC
        ISET=1
      ELSE
        RNORM=GSET
        ISET=0
      ENDIF
      RETURN
      END

!_________________________________________________________________      
      REAL (KIND=8) FUNCTION RAN1(IDUM)
!_________________________________________________________________      
! generates uniform random numbers U(0,1)
! uses 3 linear congruential generators and "shuffles" #s
! It's portable between machines with different precision

      integer m1,m2,m3,ia1,ia2,ia3,ic1,ic2,ic3,iff,ix1,ix2,ix3,j,idum
      real (KIND=8)  rm1,rm2,r(97)
      PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247d-6)
      PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773d-6)
      PARAMETER (M3=243000,IA3=4561,IC3=51349)
      DATA IFF /0/
      IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
        IFF=1
        IX1=MOD(IC1-IDUM,M1)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX2=MOD(IX1,M2)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX3=MOD(IX1,M3)
        DO 11 J=1,97
          IX1=MOD(IA1*IX1+IC1,M1)
          IX2=MOD(IA2*IX2+IC2,M2)
          R(J)=(DBLE(IX1)+DBLE(IX2)*RM2)*RM1
11      CONTINUE
        IDUM=1
      ENDIF
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      IX3=MOD(IA3*IX3+IC3,M3)
      J=1+(97*IX3)/M3
       IF(J.GT.97.OR.J.LT.1) then
       write(*,*) ' ERROR ENCOUNTERED IN RANDOM NUMBER GENERATOR' 
       write(12,*) ' ERROR ENCOUNTERED IN RANDOM NUMBER GENERATOR'
       STOP
       endif
      RAN1=R(J)
      R(J)=(DBLE(IX1)+DBLE(IX2)*RM2)*RM1
      RETURN ; END

!--------------------------------------------------------------------------------
      SUBROUTINE MISTAKE(I,RECORD)
!--------------------------------------------------------------------------------
      INTEGER I,RECORD
      IF(I>0) THEN
          WRITE(*,'(1x,a39,i5)') 'ERROR:  There is a bad field in record ',RECORD
          WRITE(12,'(1x,a39,i5)') 'ERROR:  There is a bad field in record ',RECORD
        ELSEIF(I<0) THEN
          WRITE(*,'(1x,a53,i5)') 'ERROR:  Either there are not enough fields in record ',RECORD
          WRITE(*,'(1x,a42,i5)') '        or the end of the file was reached'
          WRITE(12,'(1x,a53,i5)') 'ERROR:  Either there are not enough fields in record ',RECORD
          WRITE(12,'(1x,a42,i5)') '        or the end of the file was reached'
      ENDIF
      STOP ; END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! FUNCTION MINIMIZATION SUBROUTINES !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!--------------------------------------------------------------------------------
      SUBROUTINE SEARCH(P,X,O,PB,NPARM,SEED,MAX_ITER,CHECKFLAG,PDEV,DO_PRINT)
!--------------------------------------------------------------------------------
      USE PARAMETERS, ONLY : PARM_EST, PARM_SPECS, PARM_KEY, PARM_VAR, SPACE, ENDL
      INTEGER :: I,II,J,K,MP,NDIM,NPARM,NPP,NTRY,GOOD,SEED,ITER,CHECKFLAG,MAX_ITER,DO_PRINT,JMIN,JMAX,IMAX,HMIN,HMAX,THE_END=0, &
                 FOOTNOTE
      REAL (KIND=8) :: p(NPARM+1,NPARM),x(NPARM+1),O(NPARM+5),ftol,RSS,RNORM,PDEV,PB(NPARM),AA(500,500),BB(500,500),H,H2,RAN1, &
                       COVARIANCE(500,500),delta
      PARAMETER (NPP=500)
      EXTERNAL OBJECTIVE,RNORM,RAN1
      MP=NPARM+1 ; ndim=NPARM ; ntry=0 ; FTOL=1.d-6

! Initialize routine. If some parameters are set to the value of a
! preceding estimable parameter this step is needed to assign their
! proper values (otherwise routine starts off with the best guess
! specified for that parameter).
      CALL OBJECTIVE(PARM_EST,PARM_EST(NPARM+1),NPARM)

! Get initial value of the objective function with redefined parameters (this
! will be the same as above call if none of the parameters are linked)
      CALL OBJECTIVE(PARM_EST,PARM_EST(NPARM+1),NPARM)
      IF(DO_PRINT>0) THEN
        IF(PARM_EST(NPARM+1)<1.0E+7) THEN
            WRITE(*,'(1x,a42,f12.4)')  'Initial value of the OBJECTIVE function = ',PARM_EST(NPARM+1)
            WRITE(12,'(1x,a42,f12.4)')  'Initial value of the OBJECTIVE function = ',PARM_EST(NPARM+1)
          ELSE
            WRITE(*,'(1x,a42,e12.4)')  'Initial value of the OBJECTIVE function = ',PARM_EST(NPARM+1)
            WRITE(12,'(1x,a42,e12.4)')  'Initial value of the OBJECTIVE function = ',PARM_EST(NPARM+1)
        ENDIF
        CALL RESULTS(THE_END) ! output for initial guesses (for error checking)
      ENDIF
      IF(MAX_ITER.LE.0) RETURN
! {Belur's Optimized Step-Size Random Search}
      CALL OSSRS(PARM_EST,RSS,NPARM,SEED,DO_PRINT)

! {The downhill simplex method must be started with NPARM+1 (MP) points--
!  defining an initial simplex (a total of NP*NP+NP values).  If you
!  think of one of these points as the initial point (P[1,j], j=1,np), the
!  second point (P[2,j], j=1,np) would be obtained adding to it the quantity 
!  lambda[1] (which is the approximate length scale associate with parameter 1.  
!  The third point would be obtained by adding lambda[2] to P[1,2] and so on.  

210   DO 1 I=1,NPARM+1
        DO 1 J=1,NPARM
          IF(j.NE.(i-1)) THEN
              P(i,j)=PARM_EST(j)
            ELSE
!             {do random search in vicinity of current parameter estimates}
              P(i,j)=PARM_EST(j)*dEXP(PDEV*RNORM(SEED))
              K=PARM_KEY(j)
              IF(P(i,j)>PARM_SPECS(k,3) .OR. P(i,j)<PARM_SPECS(k,1)) THEN
                P(i,j)=PARM_SPECS(k,1)+(PARM_SPECS(k,3)-PARM_SPECS(k,1))*RAN1(SEED) !vertex can't have value outside parameter range
              ENDIF
          ENDIF
1     CONTINUE
      DO 2 I=1,MP
          do 3 j=1,NPARM
              x(j)=p(i,j)
3         continue
          CALL OBJECTIVE(X,O(i),NPARM)
2     continue
      RSS=O(1)


! {let amoeba do its thing}
211   CALL AMOEBA(p,O,mp,NPARM,ndim,ftol,OBJECTIVE,iter)

! {find the vertex with the lowest (best) residual sum of squares}
      do 4 i=1,mp
        do 5 j=1,NPARM
          x(j)=p(i,j)
5       continue
        CALL OBJECTIVE(X,O(i),NPARM)
        if(O(i).lt.rss) then
          rss=O(i)
          DO 290 J=1,NPARM
290         PARM_EST(J)=x(J)
          PARM_EST(NPARM+1)=RSS
        endif
4     continue

! {check for convergence}
      NTRY=NTRY+1
      CALL CHECK(parm_est,ntry,good,NPARM,CHECKFLAG)
      IF(DO_PRINT>0) THEN
        IF(RSS<1.0E+7) THEN
            WRITE(*,'(1x,a17,i5,a23,f12.4)')  'Amoeba restart # ',NTRY,', Objective function = ',RSS
            WRITE(12,'(1x,a17,i5,a23,f12.4)')  'Amoeba restart # ',NTRY,', Objective function = ',RSS
          ELSE
            WRITE(*,'(1x,a17,i5,a23,e12.4)')  'Amoeba restart # ',NTRY,', Objective function = ',RSS
            WRITE(12,'(1x,a17,i5,a23,e12.4)')  'Amoeba restart # ',NTRY,', Objective function = ',RSS
        ENDIF
        IF(NPARM<7) WRITE(*,*) '         ',(PARM_EST(J),J=1,NPARM)
      ENDIF
      IF(NTRY.EQ.MAX_ITER) MAX_ITER=-1
      IF(DO_PRINT>1 .AND. (GOOD==1 .OR. MAX_ITER<0)) THEN ! get a numerical approximation to the variance_covariance matrix
        WRITE(12,'(/)')
        WRITE(12,*) 'FIRST DERIVATIVE TEST'
        write(12,'(22(A13))') ('=============',I=1,3),'=======      '
        WRITE(12,'(8X,3A12)') '     -h     ','  central   ','     +h     '
        WRITE(12,'(8X,3A12)') ' ---------- ',' ---------- ',' ---------- '
        ! get central difference approximation to Hessian (lower triangular form)
        ! note: while in one dimension the point is a minima if the 2nd deriv. is
        ! continuous on the interval and its 2nd derivative is positive, this is not necessarily
        ! the case in multiple dimensions. Consider quadratic x*x - 2xy + y = 0. Has minimum of zero wherever x=y=1/2
        ! dQ/dx = 2x - 2y , dQ/dy = 1 - 2x , dq/dxdx= 2 , dQ/dydy = 0 , dQ/dxdy = -2
        ! note that, by the 1d definition, we would call this a minimum, but in fact it is
        ! a saddle point as attested by the proper 2nd derivative test:
        ! let G = h*h*Qxx + 2hkQxy + k*k*Qyy:
        !           if G > 0 for all small h,k then minimum
        !           if G < 0 for all small h,k then maximum
        !           if G < 0 for some and > 0 for others then saddle point
        ! in this case G=h*h*2 + 2hk*-2 + k*k*0, so sign depends on relative size of h and k
        ! Additionally, often with noisy functions we get a bifurcation that is infact
        ! the lowest point but has nonzero derivatives and meaningless Hessian elements.
        !
        ! Another issue is that the magnitude of the covariances depends on the value
        ! assumed for the weights given to the data.
        !
        ! The point? Don't trust the variance-covariance matrix much, especially the
        ! absolute magnitude.
        DO I=1,NPARM ;  X(I)=PARM_EST(I) ; END DO
        DO I=1,NPARM
          DO J=1,I
            IF(I==J) THEN
                DELTA=1.0D-2
555             H=ABS(PARM_EST(I))*DELTA
                X(I)=PARM_EST(I)+H ; CALL DONOTHING(X(I)) ; H=X(I)-PARM_EST(I) ! trick to reduce finite precision error
                X(I)=PARM_EST(I) + H ; CALL OBJECTIVE(X,O(1),NPARM)
                X(I)=X(I)+H ; CALL OBJECTIVE(X,O(4),NPARM)
                X(I)=PARM_EST(I) - H ; CALL OBJECTIVE(X,O(3),NPARM)
                X(I)=X(I)-H ; CALL OBJECTIVE(X,O(5),NPARM)
                X(I)=PARM_EST(I)   ; CALL OBJECTIVE(X,O(2),NPARM) ! This is computed last to re-establish x(i)=parm_est(i)
                IF(DELTA<=1.0D-2) THEN
                  ! first derivative test for minimum (should be near zero)
                  PB(I)=(O(4)-O(5)+8.0*(O(1)-O(3)))/20.0/H
                  WRITE(12,'(1X,I4,A3,20(d11.4,1X))') I,' : ',(O(2)-O(5))/(2.0*H),PB(i),(O(4)-O(2))/2.0/H
                ENDIF
                ! second derivative (should be positive)
                AA(I,I)=(O(1)-2.0*O(2)+O(3))/H/H
                IF(AA(I,I)<=0 .and. DELTA<=1.0E-2) THEN
                  ! If point is a well-behaved minimum then 2nd deriv. should all > 0, otherwise negative variance will result
                  DELTA=DELTA*2.0 ; GO TO 555
                END IF
              ELSE
                ! cross deriative (can be positive, negative or zero)
                DELTA=1.0D-2
                H=ABS(PARM_EST(I))*DELTA ; H2=ABS(PARM_EST(J))*DELTA
                X(I)=PARM_EST(I)+H ; CALL DONOTHING(X(I)) ; H=X(I)-PARM_EST(I) ! trick to reduce finite precision error
                X(J)=PARM_EST(J)+H2; CALL DONOTHING(X(J)) ; H2=X(J)-PARM_EST(J) ! trick to reduce finite precision error
                X(I)=PARM_EST(I)+H ; X(J)=PARM_EST(J)+H2 ; CALL OBJECTIVE(X,O(1),NPARM)
                X(J)=PARM_EST(J)-H2 ; CALL OBJECTIVE(X,O(2),NPARM)
                X(I)=PARM_EST(I)-H ; X(J)=PARM_EST(J)+H2 ; CALL OBJECTIVE(X,O(3),NPARM)
                X(J)=PARM_EST(J)-H2 ; CALL OBJECTIVE(X,O(4),NPARM)
                AA(I,J)=(O(1)-O(2)-O(3)+O(4))/4/H/H2
                X(j)=PARM_EST(J) ;  X(i)=PARM_EST(i) ! re-establish x=parm_est
            ENDIF
          END DO
        END DO

        WRITE(12,'(/)')
        WRITE(12,*) 'HESSIAN MATRIX'
        write(12,'(22(A13))') ('=============',I=1,MIN(PARM_KEY(0),20)),'=======      '
       IMAX=0 ; JMAX=0
        DO WHILE (IMAX<NPARM)
           IMAX=MIN(NPARM,20+IMAX) ; JMIN=JMAX+1
           DO II=0,(IMAX-1)/20
             HMIN=20*II+1 ; HMAX=HMIN+19  ; JMAX=JMIN+19 ; IF(JMAX>IMAX) JMAX=IMAX
             DO I=JMIN,JMAX
               IF(I==JMIN) WRITE(12,'(7X,20(I8,5X))') (J,J=HMIN,JMAX)
               WRITE(12,'(1X,I4,A3,20(d12.4,1X))') I,' : ',(AA(I,J),J=HMIN,MIN(HMAX,I))
             END DO
             WRITE(12,*)
           END DO
           WRITE(12,*)
        END DO

! {invert Hessian to get variance-covariance matrix for parameters}
        bb=0 ; DO I=1,NPARM ; bb(i,i)=1 ; ENDDO
        CALL GAUSSJ(AA,NPARM,NPP,BB,NPARM,NPP)

! {get covariance matrix (in lower triangular form)}
        DO I=1,NPARM
          DO J=1,I
            IF(I==J) THEN
                PARM_VAR(J)=AA(J,J) ; COVARIANCE(J,J)=AA(J,J)
              ELSE
                COVARIANCE(i,j)=AA(i,j)
            ENDIF
          END DO
        END DO
        WRITE(12,*) 'CORRELATION MATRIX OF PARAMETERS'
        write(12,'(22(A5))') ('=====',I=1,MIN(PARM_KEY(0),20)+1),'===  '
        IMAX=0 ; JMAX=0
        DO WHILE (IMAX<NPARM)
           IMAX=MIN(NPARM,20+IMAX) ; JMIN=JMAX+1
           DO II=0,(IMAX-1)/20
             HMIN=20*II+1 ; HMAX=HMIN+19  ; JMAX=JMIN+19 ; IF(JMAX>IMAX) JMAX=IMAX ; FOOTNOTE=0
             DO I=JMIN,JMAX
               DO J=HMIN,I
                 IF(AA(I,I)<0.0 .OR. AA(J,J)<0.0) THEN
                     X(J)=9.0 ; FOOTNOTE=1 ! negative variances make no sense
                   ELSEIF(AA(I,J)==0) THEN ; X(J)=0.0 ! if covariance really zero then correlation is zero no matter what
                   ELSEIF(AA(I,I)==0 .OR. AA(J,J)==0) THEN
                     X(J)=8.0 ; FOOTNOTE=1 ! zero variances with nonzero covariance makes no sense
                   ELSE ; X(J)=AA(I,J)/SQRT(AA(I,I))/SQRT(AA(J,J)) ; IF(ABS(X(J))>1.0) X(J)=X(J)/ABS(X(J))
                 ENDIF
               END DO
               IF(I==JMIN) WRITE(12,'(7X,20(I4,1X))') (J,J=HMIN,JMAX)
               WRITE(12,'(1X,I4,A3,20(F4.1,1X))') I,' : ',(X(J),J=HMIN,MIN(HMAX,I))
             END DO
             IF(FOOTNOTE>0) THEN
               WRITE(12,*) '    Values of 9.0 indicate that there were negative variances'
               WRITE(12,*) '    Values of 8.0 indicate that the covariance was nonzero but'
               WRITE(12,*) '      one or both variances were zero'
             ENDIF
             WRITE(12,*)
           END DO
           WRITE(12,*)
        END DO
        WRITE(12,*) 'COVARIANCE MATRIX OF PARAMETERS'
        write(12,'(22(A13))') ('=============',I=1,MIN(PARM_KEY(0),20)),'=======      '
       IMAX=0 ; JMAX=0
        DO WHILE (IMAX<NPARM)
           IMAX=MIN(NPARM,20+IMAX) ; JMIN=JMAX+1
           DO II=0,(IMAX-1)/20
             HMIN=20*II+1 ; HMAX=HMIN+19  ; JMAX=JMIN+19 ; IF(JMAX>IMAX) JMAX=IMAX
             DO I=JMIN,JMAX
               IF(I==JMIN) WRITE(12,'(7X,20(I8,5X))') (J,J=HMIN,JMAX)
               WRITE(12,'(1X,I4,A3,20(d12.4,1X))') I,' : ',(COVARIANCE(I,J),J=HMIN,MIN(HMAX,I))
             END DO
             WRITE(12,*)
           END DO
           WRITE(12,*)
        END DO
        RETURN
       ELSEIF(GOOD==1 .OR. MAX_ITER<0) THEN
        RETURN
      ENDIF

! {Update the results report every 5 iterations}     
      IF(DO_PRINT>0 .AND. MOD(NTRY,5).EQ.0) CALL RESULTS(THE_END)      ! A user-defined subroutine

! {restart the minimization routine at the point where it claims to have
!  found a minimum by setting P(1,J)=P(minimum RSS) and then reinitializing
!  the remaining NPARM number of simplex vertices}
      GO TO 210
      END
      

! ----------------------------------------------------      
      SUBROUTINE DONOTHING(TEMP)
      REAL (KIND=8) :: TEMP
      RETURN; END

! ----------------------------------------------------      
      SUBROUTINE CHECK(parmest,ntry,good,NP,NFLAG)
! check for convergence by comparing the results of the last n iterations (up to 10)
      USE PARAMETERS , ONLY : NPAR
      IMPLICIT NONE      
      INTEGER I,J,K,NP,NTRY,GOOD,NFLAG
      REAL (KIND=8) :: parmest(NP+1),parmold(10,NPAR),r
      IF(NFLAG.GT.10) NFLAG=10
      GOOD=0
      
      IF(ABS(PARMEST(NP+1)).LT.1.d-8) GOTO 300

      if(ntry.le.NFLAG) then
        DO 1 J=1,NP+1
1         parmold(ntry,J)=parmest(J)
        return
      endif

      do 2 J=1,NP+1
        DO 3 I=1,NFLAG-1
3         parmold(I,J)=parmOLD(I+1,J)
        PARMOLD(NFLAG,J)=PARMEST(J)
2     continue

!  Determine whether the value of the objective function varies by less
!  than 0.5% among the last NFLAG iterations.  If not, then return.  If so, then
!  apply the same criteria to the parameter estimates.
      DO 200 K=NP+1,1,-1
        do 100 i=1,NFLAG
          do 100 J=i,NFLAG
              IF(PARMOLD(J,K).NE.0.) THEN
                  r=parmold(i,K)/parmold(J,K)
                  IF(R.LT.0.995d0.or.R.GT.1.005d0) return
              ENDIF
100       continue
200   CONTINUE
300   GOOD=1
      return
      end
      
!--------------------------------------------------------------------------------
      SUBROUTINE NEWTON(XVEC,NV,CHECK,FUNC)
! Solves any real system of NV (one or two) equations with NV unknowns 
! (X,Y). Calls a function FUNC(NV,XVEC,FVEC,DF,I) where FVEC is the vector
! function whose zeroes are being sought and DF is the Jacobian (derivatives 
! with respect to XVEC).  The routine terminates when XVEC changes by  
! less than the TOLERANCE variable or after 20 iterations.
! Note the counting variable I in the list of function routine call variables.  
! This can be used in an IF-THEN statement to intialize certain constants 
! during the first iteration rather than recomputing them for all 20 calls.
! The logical variable CHECK returns a value of .TRUE. if the routine
! converges within 20 iterations.
!--------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: I,NV
      REAL (KIND=8) :: XVEC(NV),FVEC(2),DF(2,2),DELTA_X,DELTA_Y,TOLERANCE,DIVISOR,XOLD,YOLD,HI,LO,XINP
      EXTERNAL FUNC
      LOGICAL CHECK
      CHECK=.TRUE. ; TOLERANCE=0.00001
      IF(NV==1) THEN  ! use bisection
          ! starting values for bisection
          I=1; XINP=XVEC(1) ; XVEC(1)=1.0D-7 ; CALL FUNC(NV,XVEC,FVEC,DF,I)
          IF(FVEC(1)<0) THEN
              LO=XVEC(1) ; XVEC(1)=5.0 ; CALL FUNC(NV,XVEC,FVEC,DF,I) ; HI=XVEC(1)
            ELSE
              HI=XVEC(1) ; XVEC(1)=5.0 ; CALL FUNC(NV,XVEC,FVEC,DF,I) ; LO=XVEC(1)
          ENDIF
          XVEC(1)=XINP ; CALL FUNC(NV,XVEC,FVEC,DF,I)
          IF(FVEC(1)<0.0) THEN ; LO=XVEC(1) ; ELSE ; HI=XVEC(1) ; ENDIF
          ! bisection iterations
          DO I=1,100
             XVEC(1)=(LO+HI)/2.0 ; CALL FUNC(NV,XVEC,FVEC,DF,I)
             IF(FVEC(1)==0.0 .OR. ABS(HI-LO)<TOLERANCE) THEN ; RETURN
               ELSEIF(FVEC(1)<0)                        THEN ; LO=XVEC(1)
               ELSE                                          ; HI=XVEC(1)
             ENDIF
!             XOLD = XVEC(1) ; XVEC(1) = XVEC(1)-FVEC(1)/DF(1,1)
!             IF(ABS(XVEC(1)-XOLD) < TOLERANCE) RETURN
          END DO
        ELSE
          DO 2 I=1,20
            CALL FUNC(NV,XVEC,FVEC,DF,I) ; DIVISOR=DF(1,1)*DF(2,2)-DF(1,2)*DF(2,1)
            DELTA_X=(FVEC(2)*DF(1,2)-FVEC(1)*DF(2,2))/DIVISOR ; DELTA_Y=(FVEC(1)*DF(2,1)-FVEC(2)*DF(1,1))/DIVISOR
            XOLD = XVEC(1) ; YOLD=XVEC(2)
            XVEC(1) = XVEC(1)+DELTA_X ; XVEC(2) = XVEC(2)+DELTA_Y
            IF(ABS(XVEC(1)-XOLD) < TOLERANCE .AND. ABS(XVEC(2)-YOLD) < TOLERANCE) RETURN
2         CONTINUE
      ENDIF
      CHECK=.FALSE.
      RETURN ; END


!--------------------------------------------------------------------------------
      SUBROUTINE CRAMERS_RULE(x,y,a,b,c,d,k,j)
! Solves two linear equations, ax+by=k and cx+dy=j , for x and y using Cramer's rule
!--------------------------------------------------------------------------------
      REAL (KIND=8) :: X,Y,A,B,C,D,K,J,DIVISOR
      DIVISOR=A*D-B*C ; X = (D*K-B*J)/DIVISOR ; Y = (A*J-C*K)/DIVISOR
      RETURN ; END
      
!--------------------------------------------------------------------------------
      SUBROUTINE VARIANCE_RECURSION(NPRIME,S,SIGMA2,X)
! A stable algorithm for computing the variance and sum of a series if observations
! in a single pass      
!--------------------------------------------------------------------------------
      INTEGER NPRIME,N
      REAL (KIND=8) :: MU,SIGMA2,X,S
      N=NPRIME-1 ; S=S+X ; MU=S/DBLE(NPRIME)
      IF(N==0) THEN
          SIGMA2=0
        ELSE
          SIGMA2=N*SIGMA2/NPRIME+(MU-X)**2/N
      ENDIF
      RETURN ; END

!---------------------------------------------------------------------------------
      SUBROUTINE OSSRS(X,FUNC,N,SEED,DO_PRINT)
! Finds the orbital elements corresponding to the minimum of the function 'func'
! using random search in a larger interval determined by allowed ranges.
! This particular version was adapted from Fortran 77 code written by
! Sheela Belur.  The algorith was published as
!        Belur S.V., An Optimized Step-Size Random Search,
!        Computer Methods in Appl.Mech and Engg. Vol 19, pp 99-106
!---------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: NX,N,NFE,NOS,J,SEED,NRPT,MN,INP,MAXEVAL,DO_PRINT
      REAL (KIND=8) :: X(N),Y(N),AN(5,N),F(5),SS,A,AMDA,EPS,R(N),DIFF,FUNC,RNORM
      EXTERNAL RNORM
      NFE=0 ; NOS=0 ; EPS=0.0001 ; NRPT=10 ; SS=.2 ; Y=X ; MAXEVAL=1000
      DO J=1,N
        AN(1,J)=X(J)
      END DO
      INP=1 ; CALL FUN(X,AN,F,N,NFE,INP)
      IF(DO_PRINT>0) THEN ;
        IF(F(INP)>1.0E+10) THEN ; WRITE(*,'(1x,A20,2X,E12.5)') 'OSSRS initial search',F(INP)
          ELSE ; WRITE(*,*) 'OSSRS initial search',F(INP); ENDIF
      ENDIF
1     CONTINUE
      NX=3
      DO J=1,N
!       R(J)=SS*RNORM(SEED) ! Use when parameter range spans 1 unit (as with transformations)
        R(J)=SS*RNORM(SEED)*Y(J) ! Use when parameter range depends on scale of parameter (no transformation)
        AN(2,J)=AN(1,J)-R(J) ; AN(3,J)=AN(1,J)+R(J) ! When parameter range spans 1 unit (as with transformations)
      END DO
      CALL FUN(X,AN,F,N,NFE,2) ; CALL FUN(X,AN,F,N,NFE,3)
      A=F(2)-2*F(1)+F(3) ; IF(A<=0.0) GOTO 2
      AMDA=-0.5*(F(3)-F(2))/A ; NX=4
      DO J=1,N
        AN(4,J)=AN(1,J)+AMDA*R(J)
      END DO
      CALL FUN(X,AN,F,N,NFE,4)
2     CONTINUE
      MN=1
      DO J=1,NX
        IF(F(J)<F(MN)) MN=J
      END DO
      IF(F(MN)==0) GOTO 10
      DIFF=ABS(F(1)-F(MN)) ; IF(DIFF==0.0) GOTO 5
      NOS=0
      IF(DO_PRINT>0) THEN
        IF(F(INP)>1.0E+10) THEN
            WRITE(*,'(10x,A20,2X,E12.5,2X,I3)') 'OSSRS function value',f(mn),nfe
          ELSE
            WRITE(*,*) 'OSSRS function value',f(mn),nfe
         ENDIF
      ENDIF
      IF(DIFF <= EPS .OR. NFE > MAXEVAL) GOTO 10
      DO J=1,N
        AN(1,J)=AN(MN,J)
      END DO
      F(1)=F(MN)
      GOTO 1
5     NOS=NOS+1
      IF(NOS.GT.NRPT) GOTO 10
      GOTO 1
10    IF(DO_PRINT>0) THEN
        IF(F(INP)>1.0E+10) THEN
            WRITE(*,'(1x,A29,2X,E12.5,2X,I3)') 'Final OSSRS function value = ',f(mn),nfe
          ELSE
            WRITE(*,*) 'Final OSSRS function value = ',f(mn),nfe ; WRITE(*,'(/)')
         ENDIF
         FUNC=F(MN)
      ENDIF
      CALL FUN(X,AN,F,N,NFE,MN)
      DO J=1,N
        X(J)=AN(MN,J)
      END DO

      RETURN ; END SUBROUTINE

      SUBROUTINE FUN(X,AN,F,N,NFE,INP)
      INTEGER N,NFE,J,INP
      REAL (KIND=8) :: X(N),F(5),AN(5,N)
      DO J=1,N
        X(J)=AN(INP,J)
      END DO
      CALL OBJECTIVE(X,F(INP),N) ; NFE=NFE+1
      RETURN ; END SUBROUTINE
