!     Last change:  CEP  29 Mar 2017    1:55 pm

!//////////////////////////////////////////////////////////////////////////////////
! *************program PRO-2BOX*********
!                version 2.01
! multi-fleet, two-stock(sex) projection software
! by Clay E. Porch  9/07/2000
!
! changes for 2.01: 1) use AMOEBA routine for computing fishing mortality rates corresponding to a given TAC
!                   2) modify penalty constaining asymptote of stock-recruitment relation not to exceed average
!		       observed recruitment for a range of years
!                   3) add Ricker spawner-recruit function to suite of MSY calculations
! changes for 3.0   1) computes MSY with two-line model when Fmax is unsustainable (i.e., when SSBmax < hinge point)
!   5/20/2016       2) outputs selectivity for each year
!                   3) includes option for steepness parameterizations of spawner-recruit curves
!                   4) implements bias correction for spawner-recruit curves with lognormal error structure
!                   5) writes ASC-ii files with bootstrap results
!                   6) Generates uncertainty in N, F and M using input standard deviations (if bootstraps unavailable)
!                   7) Adds capability to estimate autocorrelated recruitment
!                   8) Allows user to specify target SPR values
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! VARIABLE MODULE DEFINITIONS !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE CNTRL
        TYPE GRWTH
          INTEGER :: CURVE
          REAL (KIND=8) :: LINF,K,T0,M,WA,WB,DATE
        END TYPE GRWTH
        TYPE FLS ; INTEGER :: FILETYPE ; CHARACTER (LEN=50) :: FILENAME ; END TYPE FLS
        TYPE YRT ; INTEGER :: FIRST, LAST, PRJCT, DISPLAY, STORE ; END TYPE YRT
        TYPE AGES ; INTEGER :: FIRST, LAST ; END TYPE AGES
        TYPE (GRWTH), SAVE :: GROWTH(2,0:100)
        TYPE (YRT), SAVE :: YEAR, BOX
        TYPE (AGES), SAVE :: AGE(2), SEL_TYPE, WEIGHT_TYPE
        TYPE (FLS), SAVE :: N_FILE,F_FILE,W_FILE,M_FILE,S_FILE,R_FILE,X_FILE,C_FILE,D_FILE
        TYPE (FLS), SAVE :: Q_FILE,T_FILE,Nsig_FILE,Fsig_FILE,Msig_FILE,Gsig_FILE,Rmod_FILE
        INTEGER :: PATCH, SEED, N_LOOPS, N_boots, SR_BOX, PR_TYPE, N_BOX, N_SCENARIOS,YR_BOX,N_GEARS(2),MODEL_TYPE, &
                   SSB_YEAR(2),REC_YEAR(2),N_REC,SR_TYPE(2),BAD(:),OPENSRFILE=0,ERRORTYPE,SR_PENALTY(2),REC_PENALTY(2), &
                   DO_STEEPNESS(2),F_REPORT
        REAL (KIND=8) :: CI,AGE_PLUS(2),var_R(2),SEX_FRACTION(2),BIAS_COR,SPR_TARGET(3),ESTIMATE_AUTOCORRELATION,VARY_M,F_MULT
        CHARACTER (LEN=5) :: CLASS(8)
        ALLOCATABLE BAD
        DATA CLASS/'SStot','YIELD','BIO_f','BIO_t','SSnum','RECRT','Fapex','Fcurr'/
      END MODULE CNTRL

      MODULE STATISTICS
        TYPE BENCH ; REAL (KIND=8) :: FMSY, MSY, YPRMSY, SPRMSY, SPRATIOMSY,SSBMSY,FMAX,YMAX,YPRMAX,SPRMAX,SPRATIOMAX,SSBMAX,F01, &
                                      Y01,YPR01,SPR01, SPRATIO01, SSBF01, &
                                      F20, Y20,YPR20, SPR20, SSB20, F30,Y30, YPR30, SPR30, SSB30, F40, Y40, YPR40, SPR40, SSB40, &
                                      F90,YF90,YPRF90,SPRF90,SSBF90,F75FMAX,Y75FMAX,YPR75FMAX,SPR75FMAX,SSB75FMAX
        END TYPE BENCH
        TYPE (BENCH), SAVE :: BENCHMARKS(:,:)
        REAL(KIND=8), SAVE ::  R(:,:),N(:,:,:),Fp(:,:,:,:),Cp(:,:,:,:),YLD(:,:,:),SSB(:,:),SSB_MALE(:,:),SSN(:,:),SSN_MALE(:,:),   &
                               BIOMASS(:,:),F_APICAL(:,:,:),SELECTIVITY(:,:,:),FBIOMASS(:,:),TAC(:,:,:),C(:,:,:),F(:,:,:),T(:,:),  &
                               U(:,:,:),UN(:,:,:),UF(:,:,:),UM(:,:,:),UG(:,:,:),SR(1000,2),PLUSAGE(2),SEL_MOD(:,:,:,:),W(:,:,:,:), &
                               M(:,:,:),SSBa(:,:,:),MATURITY(2,0:100),F100(:,:),F_CURRENT(:,:,:),YRTEMP(2),SRTEMP(2),SSBMSYTEMP(2),&
                               MSYTEMP(2),SUM_YLD(:,:),QUOTA(:,:,:),FMAX(:,:,:),LANDED_FRACTION(:,:,:,:),DISCARD_MOD(:,:,:),       &
                               P_MALE(:,:,:),RAND(:,:,:),RANDM(:,:,:),RANDG(:,:,:),SR_mse_minimum(2),REC_MOD(:,:),Nsig(:,:),       &
                               Fsig(:,:),Msig(:,:),Gsig(:,:),BIAS_ADJ(2),SR_GUESS(5,2),F_CURRENT_RATIO
        INTEGER :: SIGNAL(2,100),SUMSIGNAL(2)
        ALLOCATABLE R,N,F,C,YLD,SSB,SSB_MALE,SSN,SSN_MALE,BIOMASS,FBIOMASS,SELECTIVITY,F_APICAL,U,RAND,BENCHMARKS,SEL_MOD,W,M,F100,&
                    F_CURRENT,Fp,Cp,TAC,T,SUM_YLD,QUOTA,FMAX,LANDED_FRACTION,DISCARD_MOD,P_MALE,REC_MOD,UN,UF,UM,UG,Nsig,Fsig,MSig,&
                    GSig,SSBa,RANDM,RANDG
      END MODULE STATISTICS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!       MAIN PROGRAM        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      PROGRAM PRO2BOX
      WRITE(*,'(7(/))')
      WRITE(*,*) '+-------------------------------------------------------------------------+'
      WRITE(*,*) '| Program PRO-2BOX, Version 3.0 (May 15, 2016)                      |'
      WRITE(*,*) '|                                                                         |'
      WRITE(*,*) '| A companion tool for VPA-2BOX that projects the abundance and mortality |'
      WRITE(*,*) '| of one or two populations. Accomodates the overlap and diffusion models |'
      WRITE(*,*) '| of stock intermixing and sex-specific analyses.                         |'
      WRITE(*,*) '|                                                                         |'
      WRITE(*,*) '| programmed by                                                           |'
      WRITE(*,*) '|   Clay E. Porch                                                         |'
      WRITE(*,*) '|   NOAA Fisheries                                                        |'
      WRITE(*,*) '|   75 Virginia Beach Drive                                               |'
      WRITE(*,*) '|   Miami, Fl 33149 (USA)                                                 |'
      WRITE(*,*) '+-------------------------------------------------------------------------+'
      WRITE(*,'(/,/)')

      CALL CONTROL    ! options controlling projections
      CALL PROJECT    ! project original point estimates and bootstrap results
      CALL RESULTS    ! summarize results for point estimates and bootstraps
      STOP
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!       SUBROUTINES         !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      INCLUDE 'PRO-2BOX.INC'

!-----------------------------------------------------------------------------
      SUBROUTINE CONTROL
! wade through comments and read the specifications
!-----------------------------------------------------------------------------
      USE CNTRL ; USE STATISTICS ; IMPLICIT NONE
      INTEGER :: INTEGERS(100),Nage(2),a,B,G,i
      REAL (KIND=8) :: REALS(100),MINSPR
      CHARACTER (LEN=30) :: CHARACTERS(100)
      SR_PENALTY=0;
      CALL GETARG(1,CHARACTERS(1))
      IF(CHARACTERS(1)=='' .OR. CHARACTERS(1)==' ') THEN
        WRITE(*,*) 'ENTER THE NAME OF THE CONTROL FILE: ' ; READ(*,*) CHARACTERS(1) ; WRITE(*,'(3/)')
      ENDIF

1     OPEN(11,FILE=CHARACTERS(1),STATUS='OLD',ERR=10)

      WRITE(*,*) 'READING CONTROL FILE'
      CALL READIT(11,1,0,0,INTEGERS,REALS,CHARACTERS) ; MODEL_TYPE=INTEGERS(1)
      CALL READIT(11,1,0,0,INTEGERS,REALS,CHARACTERS) ; N_BOX=INTEGERS(1)
           IF(N_BOX<1 .OR. N_BOX>2) THEN ; WRITE(*,*) 'ERROR: inappropriate number of areas/stocks/sexes' ; STOP ; ENDIF
      CALL READIT(11,1,0,0,INTEGERS,REALS,CHARACTERS) ; N_LOOPS=INTEGERS(1)+1 ! add deterministic
      CALL READIT(11,0,1,0,INTEGERS,REALS,CHARACTERS) ; CI=REALS(1)
      CALL READIT(11,1,0,0,INTEGERS,REALS,CHARACTERS) ; SEED=INTEGERS(1)
      CALL READIT(11,1,0,0,INTEGERS,REALS,CHARACTERS) ; PATCH=INTEGERS(1)
      CALL READIT(11,1,0,0,INTEGERS,REALS,CHARACTERS) ; YR_BOX=INTEGERS(1)
      CALL READIT(11,1,0,0,INTEGERS,REALS,CHARACTERS) ; SR_BOX=INTEGERS(1)   ! SR_BOX > 0 indicates a sex-specific run
      CALL READIT(11,0,1,0,INTEGERS,REALS,CHARACTERS) ; SEX_FRACTION(1)=REALS(1)
           IF(N_BOX==1) THEN
               MODEL_TYPE=1 ; SR_BOX=0 ; YR_BOX=1 ; SEX_FRACTION=1.0
             ELSEIF(MODEL_TYPE/=3) THEN
               SR_BOX=0 ; SEX_FRACTION=1.0 ! Not sex-specific
             ELSEIF(SR_BOX>0) THEN ! Sex-specific
               IF(SR_BOX>3) SR_BOX=3
               WRITE(*,'(/,/,/)')
               WRITE(*,'(1X,A58,A30)') 'NOTE: You have set the model type to a value of 3 in file ',CHARACTERS(1)//','
               WRITE(*,'(7x,a31)') 'making this a sex-specific run.' ; WRITE(*,*)
               WRITE(*,*) '      Press ENTER key to continue or CONTROL-C to abort.' ; READ(*,*) ; WRITE(*,*)
               IF(SEX_FRACTION(1)>=1.0d0) THEN
                  WRITE(*,*) 'ERROR: There must be a finite fraction of sex 2' ; STOP
                 ELSEIF(SEX_FRACTION(1)<=0.0d0) THEN
                  WRITE(*,*) 'ERROR: There must be a finite fraction of sex 1' ; STOP
               ENDIF
               SEX_FRACTION(2)=1.0-SEX_FRACTION(1) ; MODEL_TYPE=1
             ELSE
               WRITE(*,*) 'ERROR: A sex-specific analysis is indicated, but a non-positive sex type was entered' ; STOP
           ENDIF
      CALL READIT(11,1,0,0,INTEGERS,REALS,CHARACTERS) ; WEIGHT_TYPE%FIRST=INTEGERS(1)
      CALL READIT(11,3,0,0,INTEGERS,REALS,CHARACTERS) ; YEAR%FIRST=INTEGERS(1);YEAR%LAST=INTEGERS(2);YEAR%PRJCT=INTEGERS(3)
           YEAR%DISPLAY=YEAR%FIRST ; YEAR%FIRST=1 ; YEAR%LAST=YEAR%LAST-YEAR%DISPLAY+1 ; YEAR%PRJCT=YEAR%PRJCT-YEAR%DISPLAY+1
      IF(N_BOX==2) THEN
          CALL READIT(11,3,0,0,INTEGERS,REALS,CHARACTERS) ; AGE%FIRST=INTEGERS(1); AGE(1)%LAST=INTEGERS(2); AGE(2)%LAST=INTEGERS(3)
          IF(AGE(1)%LAST/=AGE(2)%LAST .AND. SR_BOX/=3) THEN
             WRITE(*,'(/,/,/)')
             WRITE(*,'(1X,A61)') 'Error: This is not a sex-specific run, therefore the last age'
             WRITE(*,'(7x,A61)') '       must be the same for each of the two stocks.          '
             STOP
          ENDIF
        ELSE
          CALL READIT(11,2,0,0,INTEGERS,REALS,CHARACTERS) ; AGE%FIRST=INTEGERS(1); AGE(1)%LAST=INTEGERS(2)
      ENDIF
      CALL READIT(11,2,0,0,INTEGERS,REALS,CHARACTERS) ; SEL_TYPE%FIRST=INTEGERS(1);SEL_TYPE%LAST=INTEGERS(2)
           IF(SEL_TYPE%FIRST>0 .AND. SEL_TYPE%FIRST<YEAR%DISPLAY) THEN
             WRITE(*,*) 'ERROR: First year given for selectivity computations precedes first year in data' ; STOP
           ENDIF
           SEL_TYPE%FIRST=SEL_TYPE%FIRST-YEAR%DISPLAY+1 ; SEL_TYPE%LAST=SEL_TYPE%LAST-YEAR%DISPLAY+1
      CALL READIT(11,4,0,0,INTEGERS,REALS,CHARACTERS)
           REC_YEAR(1)=INTEGERS(1) ; REC_YEAR(2)=INTEGERS(2) ; SSB_YEAR(1)=INTEGERS(3) ; SSB_YEAR(2)=INTEGERS(4)
           REC_YEAR=REC_YEAR-YEAR%DISPLAY+1 ; N_REC=REC_YEAR(2)-REC_YEAR(1)+1 ; SSB_YEAR=SSB_YEAR-YEAR%DISPLAY+1
      CALL READIT(11,1,0,0,INTEGERS,REALS,CHARACTERS) ; ERRORTYPE=INTEGERS(1)
           IF(ERRORTYPE>2 .AND. ERRORTYPE<1) ERRORTYPE=1
      WRITE(*,*) 'RECRUITMENT YEARS',(INTEGERS(I),I=1,4)
      DO B=1,N_BOX
        Nage(b)=AGE(b)%LAST-AGE(b)%FIRST+1
        CALL READIT(11,1,0,0,INTEGERS,REALS,CHARACTERS) ; N_GEARS(b)=INTEGERS(1)
        WRITE(*,*) 'N_gears',B,N_GEARS(b)
           IF(SR_BOX>0 .AND. B==2 .AND. N_GEARS(1)/=N_GEARS(2)) THEN
             WRITE(*,*) 'ERROR: Different number of fisheries for each sex.' ; STOP
           ENDIF
        CALL READIT(11,0,1,0,INTEGERS,REALS,CHARACTERS) ; AGE_PLUS(B)=REALS(1)
        WRITE(*,*) 'Plus-age',B,age_plus(b)
        CALL READIT(11,0,Nage(b),0,INTEGERS,REALS,CHARACTERS)
        DO A=AGE(b)%first,AGE(b)%last ; MATURITY(B,a)=REALS(A-age(b)%first+1) ; END DO
        DO G=0,N_GEARS(B)
          CALL READIT(11,1,8,0,INTEGERS,REALS,CHARACTERS)
          GROWTH(B,g)%CURVE=INTEGERS(1);GROWTH(B,g)%LINF=REALS(1);GROWTH(B,g)%K=REALS(2)
          GROWTH(B,g)%M=REALS(4);GROWTH(B,g)%WA=REALS(5)
          GROWTH(B,g)%WB=REALS(6);GROWTH(B,g)%DATE=REALS(7)
          GROWTH(B,g)%T0=REALS(3)+REALS(8)/12.0
        END DO
      END DO
      CALL READIT(11,1,0,1,INTEGERS,REALS,CHARACTERS) ; Q_FILE%FILETYPE=INTEGERS(1);Q_FILE%FILENAME=CHARACTERS(1)
      WRITE(*,*) 'Q_File type',Q_FILE%FILETYPE
      CALL READIT(11,1,0,1,INTEGERS,REALS,CHARACTERS) ; N_FILE%FILETYPE=INTEGERS(1);N_FILE%FILENAME=CHARACTERS(1)
      CALL READIT(11,1,0,1,INTEGERS,REALS,CHARACTERS) ; F_FILE%FILETYPE=INTEGERS(1);F_FILE%FILENAME=CHARACTERS(1)
      CALL READIT(11,1,0,1,INTEGERS,REALS,CHARACTERS) ; W_FILE%FILETYPE=INTEGERS(1);W_FILE%FILENAME=CHARACTERS(1)
      CALL READIT(11,1,0,1,INTEGERS,REALS,CHARACTERS) ; C_FILE%FILETYPE=INTEGERS(1);C_FILE%FILENAME=CHARACTERS(1)
      CALL READIT(11,1,0,1,INTEGERS,REALS,CHARACTERS) ; M_FILE%FILETYPE=INTEGERS(1);M_FILE%FILENAME=CHARACTERS(1)
      WRITE(*,*) 'M_File type',M_FILE%FILETYPE
      INTEGERS(1)=-9; CHARACTERS(1)="NULL"
      CALL READITORNOT(11,1,0,1,INTEGERS,REALS,CHARACTERS) ; S_FILE%FILETYPE=INTEGERS(1);S_FILE%FILENAME=CHARACTERS(1)
      WRITE(*,*) 'S_File type',S_FILE%FILETYPE
      INTEGERS(1)=-9; CHARACTERS(1)="NULL"
      CALL READITORNOT(11,1,0,1,INTEGERS,REALS,CHARACTERS) ; R_FILE%FILETYPE=INTEGERS(1);R_FILE%FILENAME=CHARACTERS(1)
      WRITE(*,*) 'R_File type',R_FILE%FILETYPE
      INTEGERS(1)=-9; CHARACTERS(1)="NULL"
      CALL READITORNOT(11,1,0,1,INTEGERS,REALS,CHARACTERS) ; T_FILE%FILETYPE=INTEGERS(1);T_FILE%FILENAME=CHARACTERS(1)
      WRITE(*,*) 'T_File type',T_FILE%FILETYPE
      INTEGERS(1)=-9; CHARACTERS(1)="NULL"
      CALL READITORNOT(11,1,0,1,INTEGERS,REALS,CHARACTERS) ; D_FILE%FILETYPE=INTEGERS(1);D_FILE%FILENAME=CHARACTERS(1)
      WRITE(*,*) 'D_File type',D_FILE%FILETYPE
      INTEGERS(1)=-9; CHARACTERS(1)="NULL"
      CALL READITORNOT(11,1,0,1,INTEGERS,REALS,CHARACTERS) ; Rmod_FILE%FILETYPE=INTEGERS(1);Rmod_FILE%FILENAME=CHARACTERS(1)
      WRITE(*,*) 'Rmod_File type',Rmod_FILE%FILETYPE
      INTEGERS(1)=-9; CHARACTERS(1)="NULL"
      CALL READITORNOT(11,1,0,1,INTEGERS,REALS,CHARACTERS) ; Nsig_FILE%FILETYPE=INTEGERS(1);Nsig_FILE%FILENAME=CHARACTERS(1)
      WRITE(*,*) 'Nsig_File type',Nsig_FILE%FILETYPE
      INTEGERS(1)=-9; CHARACTERS(1)="NULL"
      CALL READITORNOT(11,1,0,1,INTEGERS,REALS,CHARACTERS) ; Fsig_FILE%FILETYPE=INTEGERS(1);Fsig_FILE%FILENAME=CHARACTERS(1)
      WRITE(*,*) 'Fsig_File type',Fsig_FILE%FILETYPE
      REALS(1)=-1.
      CALL READITORNOT(11,0,1,0,INTEGERS,REALS,CHARACTERS) ; VARY_M=REALS(1)
      IF(VARY_M>1.0) THEN; WRITE(*,*) 'ERROR: The lower limit of the uniformly distributed M cannot be less than 0'
                           WRITE(*,*) '       The fraction you input is > 1.0, hence (1-input)*M is negative.' ; STOP
      END IF
      REALS(1)=-1. ;  REALS(2)=-1
      CALL READITORNOT(11,0,2,0,INTEGERS,REALS,CHARACTERS) ; BIAS_COR=REALS(1); ESTIMATE_AUTOCORRELATION=REALS(2);
      REALS(1)=0.2;  REALS(2)=0.3; REALS(3)=0.4;
      CALL READITORNOT(11,0,3,0,INTEGERS,REALS,CHARACTERS) ;
      SPR_TARGET(1)=MIN(REALS(1),REALS(2),REALS(3)); SPR_TARGET(3)=MAX(REALS(1),REALS(2),REALS(3)); SPR_TARGET(2)=SPR_TARGET(1)
      DO b=1,3; IF(REALS(b)>SPR_TARGET(1) .AND. REALS(b)<SPR_TARGET(3) ) SPR_TARGET(2)=REALS(B); END DO
      IF(SPR_TARGET(3)>1.0) THEN; WRITE(*,*) 'ERROR: SPR targets must be a fraction between 0 and 1.0' ; STOP; ENDIF
      REALS(1)=-999;
      CALL READITORNOT(11,0,1,0,INTEGERS,REALS,CHARACTERS) ;  ! F multiplier for multi-stock calculations (-999 use default, -value = fixed multiple, value = multiple of Fcurrent(1)/Fcurrent(2))
      F_MULT=REALS(1);
      IF(N_BOX<2 .OR. SR_BOX>0) F_MULT=-999; ! F_MULT does not apply when only one stock
      INTEGERS(1)=0;
      CALL READITORNOT(11,1,0,0,INTEGERS,REALS,CHARACTERS) ;  ! F_report = first age to compute F=C/N (0 means use apical)
      F_REPORT=INTEGERS(1);
      IF(SR_BOX<1 .AND. N_BOX==2 .AND. (AGE(1)%LAST /= AGE(2)%LAST)) THEN
          WRITE(*,*) 'ERROR: Oldest age must be the same for both stocks/areas' ; STOP
        ELSEIF(SEX_FRACTION(2)<0) THEN
          WRITE(*,*) 'ERROR: Sex fraction > 1' ; STOP
        ELSEIF(SEL_TYPE%LAST>YEAR%LAST) THEN
          WRITE(*,*) 'ERROR: Last year given for selectivity computations greater than last year in data' ; STOP
        ELSE
          DO B=1,N_BOX
            IF(REC_YEAR(B)<1 .OR. REC_YEAR(B)>YEAR%LAST) THEN
                WRITE(*,*) 'ERROR: Years specified for recruitment computations outside range' ; STOP
              ELSEIF(SSB_YEAR(B)<1 .OR. SSB_YEAR(B)>YEAR%LAST) THEN
                WRITE(*,*) 'ERROR: Years specified for spawning biomass computations outside range' ; STOP
            ENDIF
          END DO
      ENDIF
      ALLOCATE(BAD(N_LOOPS)) ; BAD=0
      IF(N_LOOPS<0) THEN
        OPEN(12,FILE='BAD.OUT',STATUS='UNKNOWN')
2       READ(12,*,END=22) INTEGERS(1)
          BAD(INTEGERS(1))=1
        GOTO 2
      ENDIF
22    CLOSE(11) ; CLOSE(12) ; RETURN
10    WRITE(*,'(3/)')
      WRITE(*,'(1X,A24,A30)') 'ERROR: The control file ',CHARACTERS(1)
      WRITE(*,*) 'is not in the current directory or is being used by another process.'
      WRITE(*,*) 'Use CONTROL-C to abort or enter a new name: ' ; READ(*,*) CHARACTERS(1)
      GOTO 1
      END

!-----------------------------------------------------------------------------
      SUBROUTINE READIT(ID,N_INTEGERS,N_REALS,N_CHARACTERS,INTEGERS,REALS,CHARACTERS)
!-----------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ID,R,I,C,N_INTEGERS,N_REALS,N_CHARACTERS,INTEGERS(1:N_INTEGERS)
      REAL (KIND=8) :: REALS(1:N_REALS)
      CHARACTER (LEN=1) :: CH
      CHARACTER (LEN=30) :: CHARACTERS(1:N_CHARACTERS)
1     READ(ID,'(A1)',END=10) CH
      SELECT CASE (CH)
        CASE ('#','*','!')
          GOTO 1
        CASE DEFAULT
          BACKSPACE(11)
          READ (11,*) (INTEGERS(I),I=1,N_INTEGERS),(REALS(R),R=1,N_REALS),(CHARACTERS(C),C=1,N_CHARACTERS)
          RETURN
      END SELECT
10    WRITE(*,*) 'ERROR: End of file reached before all required inputs have'
      WRITE(*,*) '       been read (you are missing something in this file).'
      STOP ; END

!-----------------------------------------------------------------------------
      SUBROUTINE READITORNOT(ID,N_INTEGERS,N_REALS,N_CHARACTERS,INTEGERS,REALS,CHARACTERS)
!-----------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ID,R,I,C,N_INTEGERS,N_REALS,N_CHARACTERS,INTEGERS(1:N_INTEGERS),TESTIT
      REAL (KIND=8) :: REALS(1:N_REALS)
      CHARACTER (LEN=1) :: CH
      CHARACTER (LEN=30) :: CHARACTERS(1:N_CHARACTERS)
1     READ(ID,'(A1)',END=10) CH
      SELECT CASE (CH)
        CASE ('#','*','!')
          GOTO 1
        CASE DEFAULT
          BACKSPACE(11)
          READ (11,*) (INTEGERS(I),I=1,N_INTEGERS),(REALS(R),R=1,N_REALS),(CHARACTERS(C),C=1,N_CHARACTERS)
          RETURN
      END SELECT
10    WRITE(*,*) 'WARNING: This is an old file that does not specify the SPR targets, bias correction and other controls.'
      WRITE(*,*) '         Defaults will be used consistent with earlier versions of this program.'
      WRITE(*,*)
      END

!-----------------------------------------------------------------------------
	SUBROUTINE PROJECT
!-----------------------------------------------------------------------------
      USE CNTRL ; USE STATISTICS ; IMPLICIT NONE
      REAL (KIND=8) :: PARM(100),p(100,99),x(100),O(100),PB(100),func
      LOGICAL :: CHECK
      INTEGER :: NPARM,B,G,GG,GGG,Y,Y1,SCENARIO,LOOP,I,IC,N_BADLOOPS=0,FCHOICE,LUP,MAX_ITER=100
      CHARACTER (LEN=1) :: CH
      EXTERNAL OBJECTIVE_F

      CALL OPENFILES

      ALLOCATE ( U(1:N_LOOPS,2,1:YEAR%PRJCT+1),UN(1:N_LOOPS,2,0:100),UF(1:N_LOOPS,2,0:100),RAND(1:N_LOOPS,2,1:YEAR%PRJCT+1), &
                 BENCHMARKS(1:N_LOOPS,2),QUOTA(N_BOX,1:100,1:YEAR%PRJCT),FMAX(N_BOX,1:100,1:YEAR%PRJCT), &
                 RANDM(1:N_LOOPS,2,1:YEAR%PRJCT+1),UM(1:N_LOOPS,2,1:YEAR%PRJCT+1), &
                 RANDG(1:N_LOOPS,2,1:YEAR%PRJCT+1),UG(1:N_LOOPS,2,1:YEAR%PRJCT+1) )
      CALL GETRANDOM(U,RAND,N_LOOPS,YEAR%PRJCT,SEED)
      CALL GETRANDOM(UG,RANDG,N_LOOPS,YEAR%PRJCT,SEED)
      IF(Nsig_FILE%FILETYPE>-8) CALL GETRANDOM2(UN,N_LOOPS,SEED)
      IF(Fsig_FILE%FILETYPE>-8) CALL GETRANDOM2(UF,N_LOOPS,SEED)
      IF(VARY_M>0) CALL GETRANDOM(UM,RANDM,N_LOOPS,YEAR%PRJCT,SEED)
      BIAS_ADJ=1.0;

      DO SCENARIO=1,N_SCENARIOS

!       read catch and F quotas
        QUOTA=0 ; FMAX=0 ! Default is no fishing
1       READ(9,'(A1)') CH
          IF(CH/='#' .AND. CH/='!' .AND. CH/='*') THEN
            BACKSPACE(9) ; READ(9,*) B
            IF(B>0) THEN ; BACKSPACE(9) ; READ (9,*) b,g,(QUOTA(b,g,y),Y=1,(YEAR%PRJCT-YEAR%LAST)) ; ELSE ; GOTO 2 ; ENDIF
          ENDIF
        GOTO 1
2       READ(9,'(A1)') CH
          IF(CH/='#' .AND. CH/='!' .AND. CH/='*') THEN
            BACKSPACE(9) ; READ(9,*) B
            IF(B>0) THEN ; BACKSPACE(9) ; READ (9,*) b,g,(FMAX(b,g,y),Y=1,(YEAR%PRJCT-YEAR%LAST)) ; ELSE ; GOTO 3 ; ENDIF
          ENDIF
        GOTO 2

!       begin loop-specific projections
3       LOOP=0
        DO LUP=1,N_LOOPS   ! N_LOOPS includes the deterministic run
          IF(LUP>1 .AND. MOD(LUP-1,100)==0) THEN
              WRITE(*,*) 'SCENARIO ',SCENARIO,',  BOOTSTRAP RUN #: ',lUP-1
          ENDIF
          IF(BAD(LUP)>0) THEN
              N_BADLOOPS=N_BADLOOPS+1 ; WRITE(*,*) 'BAD' ; CYCLE
            ELSE
              LOOP=LOOP+1
          ENDIF
          CALL ALLOCATION
          CALL GETDATA(LOOP,SCENARIO)

          IF(SCENARIO.EQ.1) then  ; CALL GETYR(LOOP) ;  endif

          DO Y=YEAR%LAST+1,YEAR%PRJCT
            Y1=Y-YEAR%LAST

!           weight at age for gears and ssb
            CALL GETWAA(Y,LOOP)

!           projected fishing mortality rate with F limits
            IF(LOOP==1 .AND. SCENARIO==1 .AND. Y==YEAR%LAST+1) WRITE(*,*) 'SETTING F LIMITS'
            DO B=1,N_BOX
              DO G=1,N_GEARS(B)
                FCHOICE=INT(100.0*Fmax(b,g,y1))
                SELECT CASE (FCHOICE)
                  CASE (-10)     ; F_APICAL(b,g,y)=BENCHMARKS(loop,b)%F01
                  CASE (-2000)   ; F_APICAL(b,g,y)=BENCHMARKS(loop,b)%F20
                  CASE (-3000)   ; F_APICAL(b,g,y)=BENCHMARKS(loop,b)%F30
                  CASE (-4000)   ; F_APICAL(b,g,y)=BENCHMARKS(loop,b)%F40
                  CASE (-100)    ; F_APICAL(b,g,y)=BENCHMARKS(loop,b)%FMAX
                  CASE (-200)    ; F_APICAL(b,g,y)=BENCHMARKS(loop,b)%FMSY
                  CASE (-99900)  ; F_APICAL(b,g,y)=F_CURRENT(loop,b,g)
                  CASE (-75)     ; F_APICAL(b,g,y)=BENCHMARKS(loop,b)%F75Fmax
                  CASE (-90)     ; F_APICAL(b,g,y)=BENCHMARKS(loop,b)%F90
                  CASE DEFAULT ; IF(SR_BOX<1) THEN ; F_APICAL(b,g,y)=Fmax(b,g,y1) ; ELSE ; F_APICAL(b,g,y)=Fmax(1,g,y1) ; ENDIF
                END SELECT
                IF(G>1 .AND. FCHOICE<0 .AND. FCHOICE/=-99900) THEN
                  WRITE(*,*) 'ERROR: CANNOT APPLY AGGREGATE BENCHMARKS TO MULTIFLEET PROJECTIONS' ; STOP
                ENDIF
              END DO ! G
            END DO ! B
            CALL GETFAA(Y) ! Partial F's at age with above F limits
            CALL GETYIELD(y) ! yield under F limits

            IF(LOOP==1 .AND. SCENARIO==1 .AND. Y==YEAR%LAST+1) WRITE(*,*) 'APPLYING CATCH LIMITS'
!           determine if limit on fishing mortality rate or catch quota is more restrictive
            SUMSIGNAL=0 ; SIGNAL=0
            DO B=1,N_BOX ; DO G=1,N_GEARS(B)
               TAC(b,g,y)=QUOTA(b,g,Y1)*1000.0  ; IF(SR_BOX>0 .AND. B==2) TAC(2,g,y)=TAC(1,g,y)
               IF(SR_BOX<1) THEN ! Not sex-specific
                   IF(TAC(b,g,y)<YLD(b,g,Y)) SIGNAL(b,g)=1 ! F limit less restrictive than TAC
                 ELSEIF(B==1) THEN ! Sex-specific
                   IF(TAC(b,g,y)<(YLD(1,g,Y)+YLD(2,g,Y))) SIGNAL(1,g)=1 ! F limit less restrictive than TAC
               ENDIF
               SUMSIGNAL(b)=SUMSIGNAL(b)+SIGNAL(b,g)
            END DO ; END DO

            DO B=1,N_BOX
              IF(SUMSIGNAL(b)>0) THEN ! find F vector that yields fleet-specific TAC's
                YEAR%STORE=Y  ; BOX%STORE=B
                ! estimate parameters
                NPARM=0
                DO G=1,N_GEARS(B)
                  IF(SIGNAL(b,g)>0) THEN ! get starting values
                      NPARM=NPARM+1
                      IF(FBIOMASS(b,y-1)>1.0) THEN ; PARM(NPARM)=TAC(b,g,y)/FBIOMASS(b,y-1)/1000.0
                        ELSEIF(TAC(b,g,y)<=0) THEN ; PARM(NPARM)=0
                        ELSE                       ; PARM(NPARM)=TAC(b,g,y)/BIOMASS(b,y-1)/1000.0
                      ENDIF
                  ENDIF
                  PARM(NPARM)=MIN(PARM(NPARM),2.0)
                END DO
                CALL ESTIMATE(OBJECTIVE_F,PARM,P,X,O,PB,NPARM,SEED,MAX_ITER)
                CALL OBJECTIVE_F(PARM,func,NPARM)
                ! assign estimated values
                I=0 ; DO G=1,N_GEARS(B) ; IF(SIGNAL(b,g)>0) THEN ; I=I+1 ; F_APICAL(B,G,Y)=PARM(I) ; ENDIF ; END DO
              ENDIF
            END DO
            IF(SR_BOX>0) THEN ; DO G=1,N_GEARS(1) ; F_APICAL(2,G,Y)=F_APICAL(1,G,Y) ; END DO ; ENDIF
            CALL GETFAA(Y)   ! new F at age matrix
            CALL GETYIELD(y) ! yield under TAC's

!           projected population abundance at beginning of next year etc...
            CALL GETSSB(Y)
            CALL GETBIOMASS(Y)
            IF(Y<YEAR%PRJCT) CALL GETN(Y,YEAR%PRJCT-YEAR%LAST,LOOP)

          END DO ! Y (years)
          
!         Write projections to files
          IF(LOOP==1 .AND. SCENARIO==1) WRITE(*,*) 'WRITING RESULTS TO FILES'
          DO B=1,N_BOX
            IC=B*20
            WRITE(1+ic,1000) SCENARIO,LOOP-1,(SSB(B,Y),Y=1,YEAR%PRJCT)
            WRITE(2+ic,1000) SCENARIO,LOOP-1,(SUM_YLD(B,Y),Y=1,YEAR%PRJCT)
            WRITE(3+ic,1000) SCENARIO,LOOP-1,(FBIOMASS(B,Y),Y=1,YEAR%PRJCT)
            WRITE(4+ic,1000) SCENARIO,LOOP-1,(BIOMASS(B,Y),Y=1,YEAR%PRJCT)
            WRITE(5+ic,1000) SCENARIO,LOOP-1,(SSN(B,Y),Y=1,YEAR%PRJCT)
            WRITE(6+ic,1000) SCENARIO,LOOP-1,(N(B,AGE(B)%FIRST,Y),Y=1,YEAR%PRJCT)
            WRITE(7+ic,1000) SCENARIO,LOOP-1,(F100(B,Y),Y=1,YEAR%PRJCT)
            IF(LOOP==1) THEN
              DO Y=1,YEAR%PRJCT; WRITE(98,1003) B,SCENARIO,y,(SSBa(B,i,Y),i=AGE(B)%FIRST,AGE(B)%LAST); END DO;
            END IF
            IF(SCENARIO==1) THEN
              WRITE(8+IC,1000) SCENARIO,LOOP-1,(F_CURRENT(LOOP,B,G),g=1,n_gears(B))
              IF(LOOP==1) THEN
                  G=INT(100*SPR_TARGET(1));GG=INT(100*SPR_TARGET(2));GGG=INT(100*SPR_TARGET(3))
                  WRITE(31+IC,1002) 'RUN #','F_MSY   ','MSY    ','Y/R_MSY  ','S/R_MSY  ','SPR_MSY  ','SSB_MSY  ',      &
                                    'F_MAX   ','Y_MAX    ','Y/R_MAX  ','S/R_MAX  ','SPR_MAX  ','SSB_MAX  ',            &
                                    'F_0.1   ','Y_0.1    ','Y/R_0.1  ','S/R_0.1  ','SPR_0.1  ','SSB_0.1  ',            &
                                    ' F',G,'%   ',  '  Y',G,'%   ',  'Y/R',G,'%   ',  'S/R',G,'%   ',  ' SS',G,'%   ',      &
                                    ' F',GG,'%   ', '  Y',GG,'%   ', 'Y/R',GG,'%   ', 'S/R',GG,'%   ', ' SS',GG,'%   ',     &
                                    ' F',GGG,'%   ','  Y',GGG,'%   ','Y/R',GGG,'%   ','S/R',GGG,'%   ',' SS',GGG,'%   ',    &
                                    ' F_90%Y/Rmax',' Y_90%Y/Rmax',' 90%_Y/Rmax ','S/R90%Y/Rmax ','SSB90%Y/Rmax',           &
                                    ' F_75%Fmax  ',' Y_75%Fmax  ',' Y/R_75%Fmax',' S/R_75%Fmax ','SSB_75%Fmax '
              ENDIF
              WRITE(31+IC,1001) LOOP-1,BENCHMARKS(LOOP,B)
            ENDIF
          END DO ! B

          DEALLOCATE (R,N,C,F,W,M,T,TAC,Fp,Cp,YLD,F_APICAL,SSB,SSB_MALE,SSN,SSN_MALE,SELECTIVITY,BIOMASS,FBIOMASS,SEL_MOD,SUM_YLD)
          DEALLOCATE (SSBa,Nsig,Fsig,Msig,F_CURRENT,F100,LANDED_FRACTION,DISCARD_MOD,P_MALE,REC_MOD)

          IF(LOOP==1 .AND. SCENARIO==1 .AND. Y==YEAR%LAST+1) WRITE(*,*) 'END FIRST SET OF PROJECTIONS'

        END DO ! Loops (bootstraps)
        IF(SCENARIO==1) OPENSRFILE=1 ;
        CLOSE(60); CLOSE(101); CLOSE(102); CLOSE(103); CLOSE(104); CLOSE(105);

      END DO ! Scenarios

      N_boots=N_LOOPS-N_BADLOOPS-1 ! first loop not a bootstrap

      write(*,*) 'The number of usuable bootstrap runs is ',N_boots
1000  FORMAT(1X,I5,1X,I5,200(1X,E12.4))
1001  FORMAT(1X,I5,200(1X,E12.4))
1002  FORMAT(1X,A5,18(1X,A12),15(3x,a3,i3,a4),16(1X,A12))
1003  FORMAT(1X,I5,1X,I5,1X,I5,200(1X,E12.4))


      RETURN ; END SUBROUTINE

!-----------------------------------------------------------------------------
      SUBROUTINE OPENFILES
! opens output files
!-----------------------------------------------------------------------------
      USE CNTRL
      INTEGER :: IC,B,I
      CHARACTER (LEN=1) :: CH
      CHARACTER (LEN=4) :: NOBOX
      CHARACTER (LEN=6) :: IDBOX(2)
      CHARACTER (LEN=11) :: FILENAME
      DATA IDBOX/'-1.OUT','-2.OUT'/
      DATA NOBOX/'.OUT'/
      WRITE(*,*) 'OPENING OUTPUT FILES' ; WRITE(*,*)
      DO B=1,N_BOX
        IC=B*20
        DO I=1,8
          IF(N_BOX==1) THEN; FILENAME=CLASS(I)//NOBOX; ELSE; FILENAME=CLASS(I)//IDBOX(B); ENDIF
          open(I+IC,file=FILENAME,status='unknown')
        END DO
        IF(N_BOX==1) THEN; FILENAME='BENCH'//NOBOX; ELSE; FILENAME='BENCH'//IDBOX(B); ENDIF
        OPEN(31+IC,FILE=FILENAME,STATUS='unknown')
      END DO
      OPEN(9,FILE=Q_FILE%FILENAME,STATUS='OLD')
1     READ(9,'(A1)') CH ; SELECT CASE (CH) ; CASE ('#','*','!') ; GOTO 1 ; CASE DEFAULT ; BACKSPACE(9)
      READ (9,*) N_SCENARIOS ; END SELECT

      OPEN(98,FILE='SSbyage.OUT',STATUS='UNKNOWN')
      OPEN(99,FILE='SELECTIVITY.OUT',STATUS='UNKNOWN')
      WRITE(98,*) "Stock Scenario year Spawning stock by age"
      WRITE(99,*) " RUN# Stock Fleet Selectivity by age"

      RETURN ; END SUBROUTINE

!-----------------------------------------------------------------------------
      SUBROUTINE GETRANDOM(U,RAND,NLOOPS,NYEARS,SEED)
!-----------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: NLOOPS,SEED,NYEARS,I,J,B
      REAL (KIND=8) :: U(NLOOPS,2,NYEARS),RAND(NLOOPS,2,NYEARS),RNORM,ran
      EXTERNAL RNORM,ran

      DO I=1,NLOOPS ; DO J=1,NYEARS ; DO B=1,2
            U(I,B,J)=RNORM(SEED)
            RAND(I,B,J)=ran(SEED)
      END DO ; END DO ; END DO

      RETURN ; END

!-----------------------------------------------------------------------------
      SUBROUTINE GETRANDOM2(U2,NLOOPS,SEED)
!-----------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: NLOOPS,SEED,I,J,B
      REAL (KIND=8) :: U2(NLOOPS,2,0:100),RNORM
      EXTERNAL RNORM,ran

      DO I=1,NLOOPS ; DO J=0,100 ; DO B=1,2
            U2(I,B,J)=RNORM(SEED)
      END DO ; END DO ; END DO

      RETURN ; END

!-----------------------------------------------------------------------------
      SUBROUTINE ALLOCATION
!-----------------------------------------------------------------------------
      USE CNTRL ; USE STATISTICS ; IMPLICIT NONE
      INTEGER :: G,Y,A
      G=MAXVAL(N_GEARS) ; Y=YEAR%PRJCT+1 ; A=MAX(AGE(1)%LAST,AGE(2)%LAST)+1
      ALLOCATE (R(2,Y),REC_MOD(2,Y),N(2,AGE(1)%FIRST:A,Y),Cp(2,G,AGE(1)%FIRST:A,Y))
      ALLOCATE (Nsig(2,AGE(1)%FIRST:A),Fsig(2,AGE(1)%FIRST:A),Msig(2,AGE(1)%FIRST:A))
      ALLOCATE (Fp(2,G,AGE(1)%FIRST:A,Y),W(2,0:G,AGE(1)%FIRST:A,Y),YLD(2,G,Y),SUM_YLD(2,Y))
      ALLOCATE (F_APICAL(2,G,Y),SSB(2,Y),SSB_MALE(2,Y),SSN_MALE(2,Y),SSN(2,Y),M(2,AGE(1)%FIRST:A,Y),P_MALE(2,AGE(1)%FIRST:A,Y))
      ALLOCATE (SSBa(2,AGE(1)%FIRST:A,Y),DISCARD_MOD(2,0:G,AGE(1)%FIRST:A))
      ALLOCATE (BIOMASS(2,Y),SELECTIVITY(2,0:G,AGE(1)%FIRST:A),SEL_MOD(2,0:G,AGE(1)%FIRST:A,YEAR%LAST+1:Y))
      ALLOCATE (fBIOMASS(2,Y),T(2,AGE(1)%FIRST:A),TAC(2,G,Y),LANDED_FRACTION(2,0:g,AGE(1)%FIRST:A,Y))
      ALLOCATE (C(2,AGE(1)%FIRST:A,Y),F(2,AGE(1)%FIRST:A,Y),F_current(n_loops,2,g),F100(2,y))
      N=0.0 ; F=0.0 ; FP=0.0 ; T=0.0 ; LANDED_FRACTION=0 ; SEL_MOD=0; REC_MOD=1.0; Fsig=0.0; Msig=0.0; Nsig=0.0;
      RETURN ; END SUBROUTINE


!-----------------------------------------------------------------------------
      SUBROUTINE GETDATA(loop,scenario)
!-----------------------------------------------------------------------------
      USE CNTRL ; USE STATISTICS ; IMPLICIT NONE
      INTEGER (KIND=4) :: Y,A,S,IRL1,loop,YR,G,O,i,draw,first,last,j,NPS,SCENARIO,READ_RECRUIT_PENALTY=0
      REAL (KIND=8) :: U1,REC,XMU,GETF,ZMAX(2,0:100),temp,rnorm,ftemp(0:200),siggy(2),rho(2),stddev(2),TL,SEL_MOD_MAX(2), &
                       R_VIRGIN(2),STEEPNESS(2),X(2)
      REAL (KIND=4) :: FT(2,0:100,1:200),NT(2,0:100,1:200),ADD(2)
      CHARACTER (LEN=1) :: CONT
      EXTERNAL GETF,RNORM,OBJECTIVE_YR

      N=0; F=0; SSB=0; W=0; M=0; SEL_MOD_MAX=0; DO_STEEPNESS=0; SR_mse_minimum=0.0001;
!
! Read historical N, F, C, W and M values as well as recruitment parameters
!
      IRL1=0 ; DO S=1,N_BOX ; IRL1=IRL1+(AGE(S)%LAST-AGE(S)%FIRST+1) ; END DO ; IRL1=IRL1*4*YEAR%LAST
      IF(LOOP==1 .and. SCENARIO==1) WRITE(*,*) 'READING HISTORICAL ABNDANCE ESTIMATES'
      SELECT CASE (N_FILE%FILETYPE)
        CASE (0) ! ASC-II
          OPEN(10,FILE=N_FILE%FILENAME,STATUS='OLD')
          DO S=1,N_BOX ; DO Y=1,YEAR%LAST
            READ(10,*) i,yr,(N(S,a,y),a=AGE(s)%first,AGE(s)%last)
            IF(i/=s .or. y/=yr-year%display+1) THEN
              WRITE(*,*) 'ERROR: Inputs in abundance file out of proper sequence'; STOP
            ENDIF
          END DO ; END DO
        CASE (-1) ! VPA-2BOX asc-ii format
          IF(LOOP==1) OPEN(101,FILE=N_FILE%FILENAME,STATUS='OLD')
          DO S=1,N_BOX ; DO a=AGE(s)%first,AGE(s)%last
            READ(101,*) j,i,yr,(N(S,a,y),Y=1,YEAR%LAST)
          END DO ; END DO
        CASE DEFAULT ! VPA-2BOX binary format
          OPEN(10,FILE=N_FILE%FILENAME,STATUS='old',ACCESS='direct',RECL=irl1)
          READ(10,REC=loop) (((NT(S,a,y),y=1,year%last),a=AGE(s)%first,AGE(s)%last),S=1,N_BOX)
          DO S=1,N_BOX ; DO A=AGE(s)%first,AGE(s)%last ; DO Y=1,year%last ; N(s,a,y)=DBLE(NT(s,a,y)) ; END DO ; END DO ;END DO
      END SELECT
      CLOSE(10)
      SELECT CASE (Nsig_FILE%FILETYPE)
        CASE (0) ! ASC-II
          IF(LOOP==1 .and. SCENARIO==1) WRITE(*,*) 'READING HISTORICAL ABUNDANCE STANDARD ERRORS'
          OPEN(10,FILE=Nsig_FILE%FILENAME,STATUS='OLD')
          DO S=1,N_BOX ;
            DO a=AGE(s)%first,AGE(s)%last
              READ(10,*) i,Nsig(S,a)
              IF(loop>1 .and. Nsig(S,a)>0) THEN ! generate normal distributed random variable
                temp=N(s,a,year%last)+Nsig(S,a)*UN(LOOP,s,a)
                IF( temp < 0.1*N(s,a,year%last) ) THEN ;        N(s,a,year%last)=0.1*N(s,a,year%last);
                  ELSEIF ( temp > 1.9*N(s,a,year%last) ) THEN ; N(s,a,year%last)=1.9*N(s,a,year%last);
                  ELSE; N(s,a,year%last)=temp;
                ENDIF
              ENDIF
            END DO
          END DO
          CLOSE(10)
        CASE DEFAULT ! do nothing
      END SELECT

      IF(LOOP==1 .and. SCENARIO==1) WRITE(*,*) 'READING HISTORICAL FISHING MORTALITY RATE ESTIMATES'
      SELECT CASE (F_FILE%FILETYPE)
        CASE (0) ! ASC-II
          OPEN(10,FILE=F_FILE%FILENAME,STATUS='OLD')
          DO S=1,N_BOX ; DO Y=1,YEAR%LAST
            READ(10,*) i,yr,(F(S,a,y),a=AGE(s)%first,AGE(s)%last)
           IF(i/=s .or. y/=yr-year%display+1) THEN
              WRITE(*,*) 'ERROR: Inputs in fishing mortality rate file out of proper sequence'; STOP
            ENDIF
          END DO ; END DO
        CASE (-1) ! VPA-2BOX asc-ii format
          IF(LOOP==1) OPEN(102,FILE=F_FILE%FILENAME,STATUS='OLD')
          DO S=1,N_BOX ; DO a=AGE(s)%first,AGE(s)%last
            READ(102,*) j,i,yr,(F(S,a,y),Y=1,YEAR%LAST)
          END DO ; END DO
        CASE DEFAULT ! VPA-2BOX binary format
          OPEN(10,FILE=F_FILE%FILENAME,STATUS='old',ACCESS='direct',RECL=irl1)
          READ(10,REC=loop) (((FT(S,a,y),y=1,year%last),a=age(s)%first,age(s)%last),S=1,N_BOX)
          DO S=1,N_BOX ; DO A=age(s)%first,age(s)%last ; DO Y=1,year%last ; F(s,a,y)=DBLE(FT(s,a,y)) ; END DO ; END DO ;END DO
      END SELECT
      CLOSE(10)
      SELECT CASE (Fsig_FILE%FILETYPE)
        CASE (0) ! ASC-II
          IF(LOOP==1 .and. SCENARIO==1) WRITE(*,*) 'READING HISTORICAL FISHING MORTALITY RATE STANDARD ERRORS'
          OPEN(10,FILE=Fsig_FILE%FILENAME,STATUS='OLD')
          DO S=1,N_BOX ;
            DO a=AGE(s)%first,AGE(s)%last
              READ(10,*) i,Fsig(S,a)
              IF(loop>1 .and. Fsig(S,a)>0) THEN ! generate normal distributed random variable
                temp=F(s,a,year%last)+Fsig(S,a)*UF(LOOP,s,a)
                IF( temp < 0.1*F(s,a,year%last) ) THEN ;        F(s,a,year%last)=0.1*F(s,a,year%last);
                  ELSEIF ( temp > 1.9*F(s,a,year%last) ) THEN ; F(s,a,year%last)=1.9*F(s,a,year%last);
                  ELSE; F(s,a,year%last)=temp;
                ENDIF
              ENDIF
            END DO
          END DO
          CLOSE(10)
        CASE DEFAULT ! do nothing
      END SELECT

      IF(LOOP==1 .and. SCENARIO==1) WRITE(*,*) 'READING HISTORICAL WEIGHTS'
      SELECT CASE (W_FILE%FILETYPE)
          CASE (0) ! ASC-II
            OPEN(10,FILE=W_FILE%FILENAME,STATUS='OLD')
            DO S=1,N_BOX ; DO G=0,N_GEARS(S) ; DO Y=1,YEAR%LAST
              READ(10,*) i,j,yR,(W(s,g,a,y),a=age(s)%first,age(s)%last)
              IF(i/=s .or. j/=g .or. y/=yr-year%display+1) THEN
                WRITE(*,*) 'ERROR: Inputs in weight file out of proper sequence'; STOP
              ENDIF
            END DO ; END DO ; END DO
          CASE DEFAULT ! VPA-2box format
            WRITE(*,*) 'ERROR: Weights at age may only be provided through an'
            write(*,*) '       ASC-ii file' ; STOP
      END SELECT ; CLOSE(10)

      IF(LOOP==1 .and. SCENARIO==1) WRITE(*,*) 'READING HISTORICAL CATCHES'
      C=0
      SELECT CASE (C_FILE%FILETYPE)
        CASE (0) ! ASC-II
          OPEN(10,FILE=C_FILE%FILENAME,STATUS='OLD')
          C=0
          DO S=1,N_BOX ; DO G=1,N_GEARS(S) ; DO Y=1,YEAR%LAST
            READ(10,*) i,j,yR,(Cp(s,g,a,y),a=age(s)%first,age(s)%last)
            IF(i/=s .or. j/=g .or. y/=yr-year%display+1) THEN
              WRITE(*,*) 'ERROR: Inputs in catch file out of proper sequence'; STOP
            ENDIF
            DO a=age(s)%first,age(s)%last
              IF(Cp(s,g,a,y)<1.0D-10) Cp(s,g,a,y)=0.1
              C(s,a,y)=C(s,a,y)+Cp(s,g,a,y)
            END DO
          END DO ; END DO ; END DO
        CASE (-1) ! VPA-2BOX asc-ii format
          IF(N_GEARS(1)==1.AND.N_GEARS(N_BOX)==1) THEN
              IF(LOOP==1) OPEN(103,FILE=C_FILE%FILENAME,STATUS='OLD')
              DO S=1,N_BOX ; DO a=AGE(s)%first,AGE(s)%last
                READ(103,*) j,i,yr,(C(S,a,y),Y=1,YEAR%LAST)
              END DO ; END DO
              DO S=1,N_BOX ; DO A=age(s)%first,age(s)%last ; DO Y=1,year%last
                Cp(s,1,a,y)=C(s,a,y)
              END DO ; END DO; END DO
            ELSE
              WRITE(*,*) 'ERROR: Only total catches at age are provided through the VPA binary'
              write(*,*) '       files, fleet-specific calculations require an ASC-ii file with'
              write(*,*) '       fleet-specific catches at age'
              STOP
          ENDIF
        CASE DEFAULT ! VPA-2BOX binary format
          IF(N_GEARS(1)==1.AND.N_GEARS(N_BOX)==1) THEN
              OPEN(10,FILE=C_FILE%FILENAME,STATUS='old',ACCESS='direct',RECL=irl1)
              READ(10,REC=loop) (((NT(S,a,y),y=1,year%last),a=age(s)%first,age(s)%last),S=1,N_BOX)
              DO S=1,N_BOX ; DO A=age(s)%first,age(s)%last ; DO Y=1,year%last
                C(s,a,y)=MAX(DBLE(NT(s,a,y)),1.0d-1) ; Cp(s,1,a,y)=C(s,a,y)
              END DO ; END DO ; END DO
            ELSE
              WRITE(*,*) 'ERROR: Only total catches at age are provided through the VPA binary'
              write(*,*) '       files, fleet-specific calculations require an ASC-ii file with'
              write(*,*) '       fleet-specific catches at age'
              STOP
          ENDIF
      END SELECT ; CLOSE(10)

      IF(LOOP==1 .and. SCENARIO==1) WRITE(*,*) 'READING HISTORICAL DISCARD FRACTIONS'
      SELECT CASE (D_FILE%FILETYPE)
          CASE (0) ! ASC-II
            OPEN(10,FILE=D_FILE%FILENAME,STATUS='OLD')
            DO S=1,N_BOX ; DO G=1,N_GEARS(S) ; DO Y=1,YEAR%LAST   ! G=0 corresponds to FISHERY WIDE
              READ(10,*) i,j,yR,(LANDED_FRACTION(s,g,a,y),a=age(s)%first,age(s)%last)
              IF(i/=s .or. j/=g .or. y/=yr-year%display+1) THEN
                WRITE(*,*) 'ERROR: Inputs in discard file out of proper sequence'; STOP
              ENDIF
              DO a=age(s)%first,age(s)%last
                 LANDED_FRACTION(s,0,a,y)=LANDED_FRACTION(s,0,a,y)+LANDED_FRACTION(s,g,a,y)*Cp(s,g,a,y) ! total number discarded
                 LANDED_FRACTION(s,g,a,y)=1.0-LANDED_FRACTION(s,g,a,y)
              END DO
            END DO ; END DO ; END DO
            DO S=1,N_BOX ; DO Y=1,YEAR%LAST   ! G=0 corresponds to FISHERY WIDE
              DO a=age(s)%first,age(s)%last
                 LANDED_FRACTION(s,0,a,y)=1.0-LANDED_FRACTION(s,0,a,y)/C(s,a,y) ! total fraction discarded
              END DO
            END DO ; END DO
          CASE DEFAULT
            WRITE(*,*) 'WARNING: Discard fractions at age were not specified in a file and are assumed to be 0.'
            LANDED_FRACTION=1.0
      END SELECT ; CLOSE(10)

      IF(LOOP==1 .and. SCENARIO==1) WRITE(*,*) 'READING NATURAL MORTALITY RATES'
      SELECT CASE (M_FILE%FILETYPE)
        CASE (0) ! ASC-II
          OPEN(10,FILE=M_FILE%FILENAME,STATUS='OLD')
          DO S=1,N_BOX
            READ(10,*) i,(M(s,a,year%last),a=age(s)%first,age(s)%last)
            DO a=AGE(s)%first,AGE(s)%last
              DO y=1,year%last
                M(s,a,y)=M(s,a,year%last)
              END DO
            END DO
            IF(i/=s) THEN ; WRITE(*,*) 'ERROR: Inputs in natural mortality rate file out of proper sequence'; STOP ; ENDIF
          END DO
        CASE (-1) ! VPA-2BOX asc-ii format
          IF(LOOP==1) OPEN(104,FILE=M_FILE%FILENAME,STATUS='OLD')
          DO S=1,N_BOX ; DO a=AGE(s)%first,AGE(s)%last
            READ(104,*) j,i,yr,(M(S,a,y),Y=1,YEAR%LAST)
          END DO ; END DO
        CASE DEFAULT ! VPA-2box format
          OPEN(10,FILE=M_FILE%FILENAME,STATUS='old',ACCESS='direct',RECL=irl1)
          READ(10,REC=loop) (((FT(S,a,y),y=1,year%last),a=age(s)%first,age(s)%last),S=1,N_BOX)
          DO S=1,N_BOX ; DO A=age(s)%first,age(s)%last ; DO Y=1,year%last ; M(s,a,y)=DBLE(FT(s,a,y)) ; END DO ; END DO ;END DO
          CLOSE(10)
      END SELECT
      IF(LOOP==1) THEN
        OPEN(70,FILE='M_values.OUT',STATUS='UNKNOWN')
      ENDIF
      DO S=1,N_BOX ;
        DO a=AGE(s)%first,AGE(s)%last
          DO y=year%last+1,year%prjct
            M(s,a,y)=M(s,a,year%last)
            IF(VARY_M>0 .AND. LOOP>1) M(s,a,y)=M(s,a,y)*(1.0-VARY_M+2.0*VARY_M*RANDM(LOOP,s,y));! projected M varies uniformly between interval M(1-VARY_M) and M(1+VARY_M)
          END DO ! y
        END DO ! a
        WRITE(70,'(1x,i4,1x,I4,1000(1X,f8.3))') LOOP,S,VARY_M,(M(s,AGE(s)%last,y),y=year%last+1,year%prjct)
        !WRITE(*,'(1x,i4,2x,1000(f8.3,1x))') S,vary_m,(M(s,AGE(s)%last,y),y=year%last+1,year%prjct)
      END DO ! s

      IF(N_box>1 .AND. SR_BOX<1) THEN
        IF(LOOP==1 .and. SCENARIO==1) WRITE(*,*) 'READING TRANSFER COEFFICIENTS'
        SELECT CASE (T_FILE%FILETYPE)
          CASE (0) ! ASC-II
            OPEN(10,FILE=T_FILE%FILENAME,STATUS='OLD')
            DO S=1,N_BOX
              READ(10,*) i,(T(s,a),a=age(s)%first,age(s)%last)
              IF(i/=s) THEN ; WRITE(*,*) 'ERROR: Inputs in transfer coefficient file out of proper sequence'; STOP ; ENDIF
            END DO
          CASE (-1) ! VPA-2BOX asc-ii format
            IF(LOOP==1) OPEN(105,FILE=T_FILE%FILENAME,STATUS='OLD')
            DO S=1,N_BOX ; DO a=AGE(s)%first,AGE(s)%last
              READ(105,*) i,yr,(FT(S,a,y),Y=1,YEAR%LAST)
              IF(i/=s .or. a/=yr) THEN
                WRITE(*,*) 'ERROR: Inputs in fishing mortality rate file out of proper sequence'; STOP
              ENDIF
              T(s,a)=DBLE(FT(s,a,year%last))
            END DO ; END DO
          CASE (1)  ! VPA-2box binary format
            OPEN(10,FILE=T_FILE%FILENAME,STATUS='old',ACCESS='direct',RECL=irl1)
            READ(10,REC=loop) (((FT(S,a,y),y=1,year%last),a=age(s)%first,age(s)%last),S=1,N_BOX)
            DO S=1,N_BOX ; DO A=age(s)%first,age(s)%last ; T(s,a)=DBLE(FT(s,a,year%last)) ; END DO ; END DO
            CLOSE(10)
          CASE DEFAULT ! no mix
            WRITE(*,*) 'WARNING: Transfer coefficients were not specified in a file and are assemed to be 0.'
            T=0.0
         END SELECT
       ELSE
        T=0
      ENDIF

      IF(LOOP==1 .and. SCENARIO==1) WRITE(*,*) 'READING SELECTIVITY MODIFIERS'
      SELECT CASE (S_FILE%FILETYPE)
          CASE (0) ! ASC-II
            OPEN(10,FILE=S_FILE%FILENAME,STATUS='OLD')
            DO y=year%last+1,year%prjct
              DO S=1,N_BOX ; DO G=1,N_GEARS(S)
                READ(10,*,END=11) i,j,(Sel_mod(s,g,a,y),a=age(s)%first,age(s)%last)
                IF(i/=s .or. j/=g) THEN ; WRITE(*,*) 'ERROR: Inputs in vulnerability file out of proper sequence'; STOP ; ENDIF
              END DO ; END DO
            END DO
11          DO I=y,year%prjct   ! set all remaining years = to last year in file (recalling that y is incremented in the above loop)
              DO S=1,N_BOX ; DO G=1,N_GEARS(S); DO a=age(s)%first,age(s)%last;
                Sel_mod(s,g,a,i)=Sel_mod(s,g,a,y-1)
              END DO ; END DO; END DO
	    END DO
          CASE (1) ! VPA-2box format
            WRITE(*,*) 'ERROR: Selectivity modifiers may only be provided through an'
            write(*,*) '       ASC-ii file' ; STOP
          CASE DEFAULT ! modifications
            WRITE(*,*) 'WARNING: Selectivity modifiers were not specified in a file and are assemed to be 1.'
            Sel_mod=1.0
      END SELECT ; CLOSE(10)

      IF(LOOP==1 .and. SCENARIO==1) WRITE(*,*) 'READING SPAWNER-RECRUIT PARAMETERS'
      SELECT CASE (R_FILE%FILETYPE)
          CASE (0) ! ASC-II
            OPEN(10,FILE=R_FILE%FILENAME,STATUS='OLD')
            DO S=1,N_BOX
              READ(10,*) i,(SR(G,s),G=1,6)
              IF(i/=s) THEN ; WRITE(*,*) 'ERROR: Inputs in spawner-recruit file out of proper sequence'; STOP ; ENDIF
              IF(SR(4,S)<0) SR_mse_minimum(s)=ABS(SR(4,S)) ! impose a minimum level of variance even if you estimate S/R relationship
              IF(SR_BOX>0) THEN ; SR(1,2)=SR(1,1) ; EXIT ; ENDIF
           END DO
          CASE (-1) ! somebody else's binary format
          CASE (1) ! VPA-2box format
            IRL1=4*N_BOX*(5+YEAR%LAST)
            OPEN(10,FILE=R_FILE%FILENAME,STATUS='old',ACCESS='direct',RECL=irl1)
            READ(10,REC=loop) ((NT(S,1,y),Y=1,5+YEAR%LAST),S=1,N_BOX)
            DO S=1,N_BOX
              DO Y=1,5+year%last ; sr(y,s)=DBLE(NT(s,1,y)) ; END DO
              SR(1,S)=-1.0*SR(1,S) ! already estimated by VPA-2box
            END DO
          CASE DEFAULT
            WRITE(*,*) 'WARNING: The spawner-recruit relationship was not specified and option 6 has been invoked.'
            SR=0; SR(1,1)=6; SR(1,2)=6;
      END SELECT ; CLOSE(10)

      IF(Rmod_FILE%FILETYPE==0) THEN
        IF(LOOP==1 .and. SCENARIO==1) WRITE(*,*) 'READING RECRUITMENT MODIFIERS'
        OPEN(10,FILE=Rmod_FILE%FILENAME,STATUS='OLD',ERR=2)
        DO S=1,N_BOX; READ(10,*) i,(Rec_mod(s,y),y=YEAR%LAST+1,YEAR%PRJCT) ; END DO
2       CLOSE(10)
      ENDIF

      DO S=1,N_BOX
      	SR_TYPE(s)=INT(ABS(SR(1,s))+0.5)
        IF(SR_TYPE(S)>10) THEN ! invoke penalty
          SR_TYPE(S)=SR_TYPE(S)-10 ; SR_PENALTY(S)=1
          IF(READ_RECRUIT_PENALTY==0) THEN
            WRITE(*,'(/,/,/,/,/,/,/,/,/,/)')
            WRITE(*,*) 'ENTER THE FIRST AND LAST YEARS OF THE RANGE OF RECRUITMENT OBSERVATIONS THAT'
            WRITE(*,*) 'WILL BE USED TO DEFINE THE MAXIMUM RECRUITMENT PENALTY (ENTER ANY NEGATIVE '
            WRITE(*,*) 'VALUE TO USE THE MAXIMUM OBSERVED RECRUITMENT): '
            READ(*,*) REC_PENALTY
            REC_PENALTY=REC_PENALTY-YEAR%DISPLAY+1
            IF(REC_PENALTY(1)>0 .AND. REC_PENALTY(1)<REC_YEAR(1))    REC_PENALTY(1)=REC_YEAR(1)
            IF(REC_PENALTY(1)>0 .AND. REC_PENALTY(1)>REC_YEAR(2))    REC_PENALTY(1)=REC_YEAR(2)
            IF(REC_PENALTY(2)>0 .AND. REC_PENALTY(2)<REC_PENALTY(1)) REC_PENALTY(2)=REC_PENALTY(1)
            IF(REC_PENALTY(2)>0 .AND. REC_PENALTY(2)>REC_YEAR(2))    REC_PENALTY(2)=REC_YEAR(2)
            READ_RECRUIT_PENALTY=1
          END IF
        END IF
        IF(SR_TYPE(S)==7 .OR. SR_TYPE(S)==8) THEN ! steepness parameterizations of Bev-Holt and Ricker
          DO_STEEPNESS(S)=6; SR_TYPE(S)=SR_TYPE(S)-DO_STEEPNESS(S)
        END IF

        IF(N_BOX>1 .AND. S==N_BOX .AND. SR_TYPE(N_BOX)/=SR_TYPE(1)) THEN
          WRITE(*,'(/,/,/,/,/,/,/,/,/,/)')
          WRITE(*,*) 'WARNING: Spawner recruit curve specified for area/stock 1 differs from that'
          WRITE(*,*) '         specified for area/stock 2. Do you wish to continue (yes/no)?'
          READ(*,*) cont ; IF(cont=='N' .OR. cont=='n') STOP
        ENDIF
      END DO
!
! fit SR curve to data to get parameters used below
!
      IF(OPENSRFILE==0) THEN ! first scenario
          IF(SR_BOX>0) THEN ; NPS=1 ; ELSE ; NPS=N_BOX ; ENDIF
          IF(SR(1,1)>0) THEN ! Estimate parameters
              CALL SR_FIT(loop)
              IF(LOOP==1) WRITE(*,*) 'ESTIMATING SPAWNER-RECRUIT PARAMETERS'
              IF(DO_STEEPNESS(1)>0) THEN ! Convert estimated alpha and beta into steepness and R0 for reporting
                X=0  ; CALL OBJECTIVE_YR(X,temp,NPS); SRTEMP(1)=SRTEMP(1)/1000. ! spawner per recruit at F=0
                DO S=1,N_BOX
                  IF(SR_TYPE(s)==1) THEN ! Beverton and Holt
                      STEEPNESS(s)=SR(2,s)*SRTEMP(S)/(4.*SR(3,S)+SR(2,S)*SRTEMP(S))
                      R_virgin(S)=(SR(2,S)*SRTEMP(S)-SR(3,S))/SRTEMP(S)
                    ELSE ! Ricker
                      STEEPNESS(S)=0.2*(SR(2,S)*SRTEMP(S))**0.8
                      R_virgin(S)=LOG(SR(2,S)*SRTEMP(S))/(SR(3,S)*SRTEMP(S))
                  ENDIF
                END DO
              ENDIF
            ELSE ! do not estimate parameters; use input values
              IF(DO_STEEPNESS(1)>0) THEN ! Convert input R_virgin and steepness into alpha and beta parameters
                X=0  ; CALL OBJECTIVE_YR(X,temp,NPS) ! spawner per recruit at F=0
                DO S=1,N_BOX
                  R_virgin(s)=SR(2,S); STEEPNESS(S)=SR(3,S); SRTEMP(S)=SRTEMP(S)/1000.0 ! parameters input as R_virgin and steepness
                  IF(SR_TYPE(S)==1) THEN ! Beverton and Holt
                      SR(2,S)=4.0*R_virgin(s)*STEEPNESS(s)/(5.0*STEEPNESS(s)-1.0); SR(3,S)=(SR(2,s)-R_virgin(s))*SRTEMP(S);
                    ELSE ! Ricker
                      SR(2,S)=( (5*STEEPNESS(s))**1.25 )/SRTEMP(S); SR(3,S)=log(SR(2,s)*SRTEMP(S))/SRTEMP(S)/R_virgin(s);
                  ENDIF
                END DO
              ENDIF
          ENDIF

        ELSE ! Don't re-estimate parameters for every scenario

          IF(LOOP==1) THEN; OPEN(60,FILE='SR_PARMS.OUT',STATUS='UNKNOWN') ; READ(60,*) ; ENDIF
          READ(60,*) I,(SR_TYPE(s),SR(2,S),SR(3,S),SR(4,S),SR(5,S),SR(5+YEAR%LAST,S),S=1,N_BOX)
          IF(DO_STEEPNESS(1)>0) THEN ! Convert steepness and R_virgin back into into alpha and beta parameters
            X=0  ; CALL OBJECTIVE_YR(X,temp,NPS) ! spawner per recruit at F=0
            DO S=1,N_BOX
              SRTEMP(S)=SRTEMP(S)/1000.0
              R_virgin(s)=SR(2,S); STEEPNESS(S)=SR(3,S) ! parameters input as R_virgin and steepness
              IF(SR_TYPE(1)==1) THEN ! Beverton and Holt
                  SR(2,S)=4.0*R_virgin(s)*STEEPNESS(s)/(5.0*STEEPNESS(s)-1.0); SR(3,S)=(SR(2,s)-R_virgin(s))*SRTEMP(S);
                ELSE ! Ricker
                  SR(2,S)=( (5*STEEPNESS(s))**1.25 )/SRTEMP(S); SR(3,S)=log(SR(2,s)*SRTEMP(S))/SRTEMP(S)/R_virgin(s);
              ENDIF
            END DO
          ENDIF

      ENDIF
      IF(ERRORTYPE==1 .AND. BIAS_COR>0) THEN; DO S=1,N_BOX; IF(SR_TYPE(s)/=6) BIAS_ADJ(s)=EXP(sr(4,s)*sr(4,s)/2.0); END DO; ENDIF
      IF(SR_BOX>0) THEN ! sex-specific analyses have only one set of spawner/recruit parameters
        DO Y=2,5+year%last ; sr(y,2)=sr(y,1) ; END DO
      ENDIF
!
! Over-ride predicted N's in last year when terminal F specified as random variable
!
      IF(SUM(Fsig)>0.0) THEN
        IF(LOOP==1 .and. SCENARIO==1) WRITE(*,*) 'OVER-RIDING TERMINAL Fs AND RECOMPUTING SUBSEQUENT Ns'
        DO S=1,N_BOX ;
          DO a=AGE(s)%first,AGE(s)%last
            IF(loop>1 .and. Fsig(S,a)>0) THEN
              temp=F(s,a,year%last)+M(s,a,year%last)
              N(s,a,year%last)=C(s,a,year%last)*TEMP/(F(s,a,year%last)*(1-EXP(-TEMP)));
            ENDIF
          END DO
        END DO
      ENDIF
!
! Over-ride predicted N's in bad corner of VPA
!
      IF(PATCH>0) THEN
        IF(LOOP==1 .and. SCENARIO==1) WRITE(*,*) 'OVER-RIDING LAST RECRUITMENTS AND RECOMPUTING SUBSEQUENT Ns AND Fs'

        ! get recruitment replacements
        DO Y=YEAR%LAST-patch+1,YEAR%LAST
          ! Deterministic part
          IF(WEIGHT_TYPE%FIRST==0) CALL GETWAA(Y-AGE(1)%FIRST,LOOP)      ! get weights from growth curve
          CALL GETSSB(y-AGE(1)%FIRST)
          XMU=0
          DO S=1,N_BOX
            rho(s)=sr(5,s) ; stddev(s)=ABS(sr(4,s))
            IF(Y==YEAR%LAST-patch+1) THEN
              IF(R_FILE%FILETYPE==0 .AND. ESTIMATE_AUTOCORRELATION<=0) THEN ; siggy(s)=sr(6,s)  ! Starting residual from single entry in ASC-ii file
                ELSE; siggy(s)=sr(5+YEAR%LAST-PATCH,s)                                          ! Starting residual from binary file entry for year before patch or estimated internally
              ENDIF
            ENDIF
            CALL RECRUITS(REC,S,Y)
            ! get deterministic number of recruits (REC)
            siggy(s)=rho(s)*siggy(s) ;
            N(S,age(s)%first,Y)=MAX(0.0d0,REC)
            XMU=XMU+MAX(0.0d0,REC)
          END DO ! S
          ! Stochastic part
          DO S=1,N_BOX
            IF(SR_BOX>0 .AND. SR_TYPE(s)/=6) THEN ! sex specific
                REC=XMU ;
                IF(SR_BOX>2) THEN ; BOX%STORE=1 ; ELSE ; BOX%STORE=SR_BOX ; ENDIF !s/r relationship assumes sex S is limiting factor
                IF(S>1) EXIT      ! one s/r relationship governs both sexes, so only need one pass through loop
              ELSE ! area/stock specific or SR_TYPE=6
                REC=N(S,age(s)%first,Y) ; BOX%STORE=S
            ENDIF
            IF(REC>0 .AND. LOOP>1) THEN
                IF(SR_TYPE(s)==6) THEN
                    DRAW=REC_YEAR(1)+INT(RAND(LOOP,BOX%STORE,Y)*N_REC)   ! draw at random from historical recruitments
                    REC=N(BOX%STORE,age(s)%first,DRAW)
                  ELSE
                    siggy(BOX%STORE)=siggy(BOX%STORE)+STDDEV(BOX%STORE)*U(LOOP,BOX%STORE,y)
                    ! STDDEV is std dev. of the random component of recruit variation
                    IF(ERRORTYPE==1) THEN
                        xmu=log(rec)                   ! mu = E[ln(R)] = ln(E[R])-mse/2 where E is expectation
                        rec=exp(xmu+siggy(BOX%STORE))  ! stochastic recruitment value with mean = deterministic value
                      ELSE
                        REC=REC+siggy(BOX%STORE)
                    ENDIF
                ENDIF
                N(S,age(s)%first,Y)=REC
                SR(5+Y,s)=SIGGY(s)
              ELSEIF(REC>0 .AND. LOOP==1 .AND. ERRORTYPE==1) THEN
                N(S,age(s)%first,Y)=REC*BIAS_ADJ(s)   ! bias-correct MLE from lognormal distribution
            ENDIF
          END DO
          IF(SR_BOX>0) THEN ; DO S=1,N_BOX ; N(S,age(s)%first,Y)=REC*SEX_FRACTION(S) ; END DO ; ENDIF  ! sex specific
        END DO
        ! determine F's and N's of cohorts with replaced recruitments from catches
        DO Y=YEAR%LAST-patch+1,YEAR%LAST
          ! F's associated with recruitments
          ADD=0
99        DO S=1,N_BOX ; N(S,age(s)%first,Y)=N(S,age(s)%first,Y)+ADD(S) ; END DO
          DO S=1,N_BOX
            IF(S==1) THEN ; O=2 ; ELSE ; O=1 ; ENDIF
            TEMP=N(S,age(s)%FIRST,Y)*(1.0-T(S,age(s)%FIRST))+N(O,age(s)%FIRST,Y)*T(O,age(s)%FIRST) ! abundance in area
            F(S,age(s)%FIRST,Y)=GETF(U1,C(S,age(s)%FIRST,Y),TEMP,M(S,age(s)%FIRST,Y)) ! F by Newton iteration from N and C
            IF(U1>0.001) THEN ! replaced recruitment too small
              ADD(S)=U1*N(S,age(s)%FIRST,Y)*(1.0-T(S,age(s)%FIRST))/TEMP+1.0
              ADD(O)=U1*N(O,age(s)%FIRST,Y)*T(O,age(s)%FIRST)/TEMP+1.0
 !             WRITE(*,*) '****redo1****', add
              GOTO 99
            ENDIF
 !           WRITE(*,*) S,1,Y, F(S,age(s)%FIRST,Y), N(S,age(s)%FIRST,Y)
          END DO
          ! N's and F's in subequent years
          DO I=Y,YEAR%LAST-1
            A=age(1)%FIRST+I-Y
            DO S=1,N_BOX ! abundance of cohort I years after recruitment
              IF(S==1) THEN ; O=2 ; ELSE ; O=1 ; ENDIF
              IF(MODEL_TYPE==1) THEN ! diffusion and no mix
                  N(S,A+1,I+1)=N(S,A,I)*(1.0-T(S,A))*EXP(-M(S,A,I)-F(S,A,I))+N(O,A,I)*T(O,A)*EXP(-M(S,A,I)-F(S,A,I))
                ELSE ! overlap
                  N(S,A+1,I+1)=N(S,A,I)*(1.0-T(S,A))*EXP(-M(S,A,I)-F(S,A,I))+N(S,A,I)*T(S,A)*EXP(-M(O,A,I)-F(O,A,I))
              ENDIF
            END DO ! S
            DO S=1,N_BOX ! F's on cohort I years after recruitment
              IF(S==1) THEN ; O=2 ; ELSE ; O=1 ; ENDIF
              TEMP=N(S,A+1,I+1)*(1.0-T(S,A+1))+N(O,A+1,I+1)*T(O,A+1)
              F(S,A+1,I+1)=GETF(U1,C(S,A+1,I+1),TEMP,M(S,A+1,I))
              IF(U1>0.001) THEN ! replaced recruitment too small
                ADD(S)=U1*N(S,A+1,I+1)*(1.0-T(S,A+1))/TEMP+1.0
                ADD(O)=U1*N(O,A+1,I+1)*T(O,A+1)/TEMP+1.0
!                WRITE(*,*) '****redo2****', add
                GOTO 99
              ENDIF
!              WRITE(*,*) S,A+1,i+1,F(S,A+1,i+1),N(S,A+1,I+1)
            END DO ! S
          END DO ! I
        END DO ! Y
      ENDIF

      IF(OPENSRFILE==0) THEN ! first scenario
        IF(LOOP==1) OPEN(60,FILE='SR_PARMS.OUT',STATUS='UNKNOWN')
        IF(LOOP==1) WRITE(60,'(1x,a5,a10,6(a14))') ' RUN#',' SR_TYPE  ','PARM_1       ','PARM_2      ','  STD_DEV.   ', &
                                                                        '     RHO     ','Last_Residual'
        IF(DO_STEEPNESS(1)>0) THEN
            DO S=1,N_BOX; ftemp(S)=R_virgin(S); ftemp(S+2)=STEEPNESS(S); END DO
          ELSE
            DO S=1,N_BOX; ftemp(S)=SR(2,S); ftemp(S+2)=SR(3,S); END DO
        ENDIF
        IF(ftemp(1)<0.01 .OR. ftemp(3)<0.01) THEN
            WRITE(60,'(1x,i4,2x,2(1X,I4,2(1X,D12.5),3(1x,f14.4)))') &
                  LOOP-1,(SR_TYPE(s),FTEMP(S),FTEMP(S+2),SR(4,S),SR(5,S),SR(5+YEAR%LAST,S),S=1,N_BOX)
          ELSE
            WRITE(60,'(1x,i4,2x,2(1X,I4,5(1x,f14.4)))') &
                  LOOP-1,(SR_TYPE(s),FTEMP(S),FTEMP(S+2),SR(4,S),SR(5,S),SR(5+YEAR%LAST,S),S=1,N_BOX)
        ENDIF
      ENDIF

!
! compute the partial F's
!
      IF(LOOP==1 .AND. (N_GEARS(1)>1 .OR. N_GEARS(2)>1)) WRITE(*,*) 'COMPUTING THE PARTIAL Fs'
      F100=0
      DO S=1,N_BOX ; DO Y=1,YEAR%LAST ; DO A=age(s)%FIRST,age(s)%LAST
        DO G=1,N_GEARS(s)
          Fp(s,g,a,y)=F(s,a,y)*Cp(s,g,a,y)/C(s,a,y)
        END DO
        F100(S,Y)=MAX(F100(S,Y),F(S,A,Y))
      END DO ; END DO ; END DO
!
! compute the relative selectivity vector (SELECTIVITY) for the multi-fleet projections
!
      IF(LOOP==1 .and. SCENARIO==1) WRITE(*,*) 'COMPUTING PROJECTION SELECTIVITIES'
      SELECTIVITY=0 ; ZMAX=0
      DO S=1,N_BOX
        DO G=1,N_GEARS(S)
          DO A=age(s)%FIRST,age(s)%LAST
            IF(SEL_TYPE%FIRST>0) THEN ! Get geometric mean of historical fishing mortalities
                DO Y=SEL_TYPE%FIRST,SEL_TYPE%LAST
                  IF(Fp(S,G,A,Y)*SEL_MOD(s,g,a,YEAR%LAST+1)>0) THEN;
                      SELECTIVITY(S,G,A)=SELECTIVITY(S,G,A)+LOG(Fp(S,G,A,Y)*SEL_MOD(s,g,a,YEAR%LAST+1))
                    ELSE
                      SELECTIVITY(S,G,A)=SELECTIVITY(S,G,A)-10; ! Set selection for this year to tiny non-zero value
                  ENDIF
                END DO ! Y
                SELECTIVITY(S,G,A)=EXP(SELECTIVITY(S,G,A)/DBLE(SEL_TYPE%LAST-SEL_TYPE%FIRST+1)) ! fleet specific F from geometric mean
              ELSE ! Use fishing mortality estimates from last year of data
                SELECTIVITY(S,G,A)=Fp(S,G,A,YEAR%LAST)
            ENDIF
            SELECTIVITY(S,0,A)=SELECTIVITY(S,0,A)+SELECTIVITY(S,G,A) ! total F
            ZMAX(S,G)=MAX(ZMAX(S,G),SELECTIVITY(S,G,A))
          END DO ! A
          F_CURRENT(LOOP,S,G)=ZMAX(S,G)
          DO A=age(s)%FIRST,age(s)%LAST
            SEL_MOD(S,0,A,YEAR%LAST+1)=SEL_MOD(S,0,A,YEAR%LAST+1)+SEL_MOD(S,g,A,YEAR%LAST+1)*F_CURRENT(LOOP,S,G) ! sum F weighted input selectivities
          END DO ! A
        END DO ! G
        DO A=age(s)%FIRST,age(s)%LAST
          ZMAX(S,0)=MAX(ZMAX(S,0),SELECTIVITY(S,0,A))
          SEL_MOD_MAX(S)=MAX(SEL_MOD_MAX(S),SEL_MOD(S,0,A,YEAR%LAST+1))
        END DO
      END DO ! S
      DO S=1,N_BOX
        DO G=0,N_GEARS(S)
          IF(SR_BOX>0.AND.G>0) F_CURRENT(LOOP,S,G)=MAX(F_CURRENT(LOOP,1,G),F_CURRENT(LOOP,2,G)) ! Sex-specfic
          IF(SEL_TYPE%FIRST>0) THEN ! normalize geometric mean F's by maximum value to get relative selectivity
              IF(N_BOX>1 .AND. SR_BOX>0) ZMAX(S,G)=MAX(ZMAX(1,G),ZMAX(2,G)) ! Sex-specfic
              DO A=age(s)%FIRST,age(s)%LAST ; SELECTIVITY(S,G,A)=SELECTIVITY(S,G,A)/ZMAX(S,G) ; END DO
              IF(SCENARIO==1) WRITE(99,1000) LOOP-1,S,G,(SELECTIVITY(S,G,A),A=age(s)%FIRST,age(s)%LAST)
            ELSE ! User-supplied relative selectivity SELMOD
              IF(G==0) THEN
                IF(N_BOX>1 .AND. SR_BOX>0) SEL_MOD_MAX(S)=MAX(SEL_MOD_MAX(1),SEL_MOD_MAX(2)) ! Sex-specfic
                DO A=age(s)%FIRST,age(s)%LAST ; SEL_MOD(S,0,A,YEAR%LAST+1)=SEL_MOD(S,0,A,YEAR%LAST+1)/SEL_MOD_MAX(S) ; END DO
              ENDIF
              DO A=age(s)%FIRST,age(s)%LAST ; SELECTIVITY(S,G,A)=SEL_MOD(S,g,A,YEAR%LAST+1) ; END DO
              IF(SCENARIO==1) WRITE(99,1000) LOOP-1,S,G,(SELECTIVITY(S,G,A),A=age(s)%FIRST,age(s)%LAST)
          ENDIF
        END DO ! G
      END DO ! S
      IF(F_MULT>-999) F_CURRENT_RATIO=F_CURRENT(LOOP,2,0)/F_CURRENT(LOOP,1,0)


!
! compute the discard modification factor (quotas are based on landings, not discards)
!
      IF(LOOP==1 .and. SCENARIO==1) WRITE(*,*) 'COMPUTING DISCARD MODIFICATIONS'
      DO S=1,N_BOX
        DO G=0,N_GEARS(S)
          FTEMP=0
          DO A=age(s)%FIRST,age(s)%LAST
            IF(SEL_TYPE%FIRST>0) THEN ! {geometric mean of historical landed fractions}
                DO Y=SEL_TYPE%FIRST,SEL_TYPE%LAST
                  IF(LANDED_FRACTION(s,g,a,y)>0.0) THEN
                      FTEMP(a)=FTEMP(a)+LOG(LANDED_FRACTION(s,g,a,y))
                    ELSE
                      FTEMP(a)=FTEMP(a)+LOG(0.00000001)
                  ENDIF
                END DO
                DISCARD_MOD(s,g,a)=EXP(FTEMP(a)/DBLE(SEL_TYPE%LAST-SEL_TYPE%FIRST+1))
              ELSE
                DISCARD_MOD(s,g,a)=LANDED_FRACTION(s,g,a,year%last)
            ENDIF
          END DO
        END DO
      END DO

!
! adjust the historical partial F's to discount discards when computing historical yield from catch at age
!
      DO S=1,N_BOX ; DO Y=1,YEAR%LAST ; DO A=age(s)%FIRST,age(s)%LAST
        DO G=1,N_GEARS(s) ; Fp(s,g,a,y)=Fp(s,g,a,y)*LANDED_FRACTION(s,g,a,y) ; END DO
      END DO ; END DO ; END DO


!
! compute historical spawning stock biomasses, fishable biomass, yields, etc..
!
      IF(LOOP==1 .and. SCENARIO==1) WRITE(*,*) 'COMPUTING HISTORICAL SSB, FISHABLE BIOMASS, ETC...'
      DO Y=YEAR%FIRST,YEAR%LAST
        IF(WEIGHT_TYPE%FIRST==0) CALL GETWAA(Y,LOOP)
        CALL GETSSB(y)
        CALL GETYIELD(y)
        CALL GETBIOMASS(y)
!        CALL GETSUM(box,year)
      END DO
!
! compute age of plus groups
      DO S=1,N_BOX
        IF(AGE_PLUS(S)>=REAL(AGE(S)%LAST)) THEN ! Set equal to values in control file
            PLUSAGE(S)=MAX(REAL(age(s)%last),AGE_PLUS(S))
          ELSE ! Compute from average weight of plus-group given for fishery 1 in weight file
            IF(GROWTH(S,1)%CURVE==1) THEN  ! von Bertalanffy/Chapman-Richards growth
                TL=(W(s,1,AGE(S)%LAST,YEAR%LAST)/GROWTH(S,1)%WA)**(1.0/GROWTH(S,1)%WB)
                PLUSAGE(S)=GROWTH(S,1)%T0-LOG((1.0-(TL/GROWTH(S,1)%LINF)**GROWTH(S,1)%M)/GROWTH(S,1)%M)/GROWTH(S,1)%K
              ELSEIF(GROWTH(S,1)%CURVE==2) THEN  ! Gompertz growth
                PLUSAGE(S)=GROWTH(S,1)%T0-LOG(LOG(GROWTH(S,1)%LINF/W(S,1,AGE(S)%LAST,YEAR%LAST)))/GROWTH(S,1)%K
            ENDIF
            PLUSAGE(S)=MAX(REAL(age(s)%last),PLUSAGE(S)-GROWTH(S,1)%DATE/12.0) ! avg age of plus group at beginning of year
        ENDIF
      END DO

!
! compute abundance at beginning of first projection year
!
      CALL GETN(year%last,year%prjct,LOOP)

      IF(LOOP==1 .and. SCENARIO==1) WRITE(*,*) 'END OF SUBROUTINE GETDATA'

1000  FORMAT(1X,I5,1X,I5,1X,I5,200(1X,f6.4))
      RETURN ; END


!-----------------------------------------------------------------------------
      REAL (KIND=8) FUNCTION GETF(OVERAGE,C,N,M)
! Estimates the fishing mortality rate
!-----------------------------------------------------------------------------
      INTEGER :: iter
      REAL(KIND=8) :: c,n,m,overage,F,FOLD,Z,S,FUNCT,DERIVATIVE

      IF(N.LE.0) THEN ! infeasible, probably a bad loop
        WRITE(*,*) 'NEGATIVE OR ZERO N' ; GETF=0 ; RETURN
       ELSEIF(C>3.0*N*(1.0-EXP(-3.0-M))/(3.0+M)) THEN ! infeasible, abundance too low to support obseved catch
        GETF=3.0 ; OVERAGE=C-3.0*N*(1.0-EXP(-3.0-M))/(3.0+M) ; RETURN
       ELSE
        Z=C/N+M ; F=C*Z/N/(1-EXP(-Z)) ! First approximation to F
      ENDIF
      FOLD=-1000.
      DO ITER=1,20
        Z=MAX(M,MIN(709.0,M+F))
        S=EXP(-Z)
        FUNCT=F-C*Z/(N*(1-S))
        DERIVATIVE=1-C*(1-S-Z*S)/(N*(1-S)**2)
        F = F - FUNCT/DERIVATIVE
        IF(ABS(F-FOLD).LE.0.0000000001 .OR. ABS(DERIVATIVE).LT.1.E-32) EXIT
        FOLD=F
      END DO
      IF(F>3.0) THEN ; F=3.0 ; ELSEIF(F<0) THEN ; F=0 ; ENDIF
      Z=M+F ; S=EXP(-Z)
      OVERAGE=C-F*N*(1.0-S)/Z
      GETF=F
      RETURN
      END

!-----------------------------------------------------------------------------
      SUBROUTINE GETWAA(y,loop)
! note: G=0 denotes ssb
!-----------------------------------------------------------------------------
      USE CNTRL ; USE STATISTICS ; IMPLICIT NONE
      INTEGER :: S,A,Y,G,loop
      REAL (KIND=8) :: Ta,TL

      DO S=1,N_BOX ; DO G=0,N_GEARS(S) ; DO A=age(s)%FIRST,age(s)%LAST
        IF((WEIGHT_TYPE%FIRST==0) .OR. (WEIGHT_TYPE%FIRST==1 .AND. Y>YEAR%LAST .AND. A==age(s)%LAST)) THEN
!          use growth curve
           IF(A==age(s)%LAST) THEN
               IF(Y>YEAR%LAST) THEN
                   Ta=PLUSAGE(S)+GROWTH(S,G)%DATE/12.0
                 ELSE
!                  this is the equilibrium approximation for historical plus-age, which is crude if F or R changes thru time
                   Ta=DBLE(age(s)%LAST)-1.0+1.0/(1.0-EXP(-m(s,a,y)-F(s,a,y)))+GROWTH(S,G)%DATE/12.0
               ENDIF
             ELSE
               Ta=DBLE(A)+GROWTH(S,G)%DATE/12.0
           ENDIF
           IF(GROWTH(S,G)%CURVE==1) THEN  ! von Bertalanffy growth
               TL=GROWTH(S,G)%LINF*(1.0-GROWTH(S,G)%M*EXP(-GROWTH(S,G)%K*(Ta-GROWTH(S,G)%T0)))**(1.0/GROWTH(S,G)%M)
               !IF(LOOP>1 .AND. Y>YEAR%LAST) TL=TL*(1.0/VARY_GROWTH+(VARY_GROWTH-1.0/VARY_GROWTH)*RANDG(LOOP,S,Y))
               W(S,g,A,Y)=GROWTH(S,G)%WA*TL**GROWTH(S,G)%WB
             ELSEIF(GROWTH(S,G)%CURVE==2) THEN  ! Gompertz growth
               W(S,g,A,Y)=GROWTH(S,G)%LINF*EXP(-1.0*EXP(-GROWTH(S,G)%K*(Ta-GROWTH(S,G)%T0)))
               !IF(LOOP>1 .AND. Y>YEAR%LAST) W(S,g,A,Y)=(1/VARY_GROWTH+(VARY_GROWTH-1.0/VARY_GROWTH)*RANDG(LOOP,S,Y))*W(S,g,A,Y)
           ENDIF
         ELSEIF(Y>year%LAST) THEN ! use weights from last year in data
           W(s,g,a,y)=W(s,g,a,year%last)
        ENDIF

      END DO ; END DO ; END DO

      RETURN ; END

!-----------------------------------------------------------------------------
      SUBROUTINE GETN(y,nyears,LOOP)
!-----------------------------------------------------------------------------
      USE CNTRL ; USE STATISTICS ; IMPLICIT NONE
      INTEGER :: A,Y,NYEARS,LOOP,O,B,DRAW,ibox
      REAL (KIND=8) :: temp(2),xmu,rec,ETA(2),RHO,STDDEV

!     project abundance of each existing cohort
      DO B=1,N_BOX
        IF(B==1) THEN ; O=2 ; ELSE ; O=1 ; ENDIF
        DO A=AGE(b)%FIRST,AGE(b)%LAST
          IF(MODEL_TYPE==1) THEN ! diffusion
              N(B,A+1,Y+1)=(N(b,A,Y)*(1.0-T(b,a))+N(o,A,Y)*T(o,a))*EXP(-M(b,A,y)-F(b,a,y))
            ELSE ! overlap
              N(B,A+1,Y+1)=N(b,A,Y)*((1.0-T(b,a))*EXP(-M(b,A,y)-F(b,a,y))+T(b,a)*EXP(-M(o,A,y)-F(o,a,y)))
          ENDIF
        END DO
      END DO

!     update age of plus group
      DO B=1,N_BOX
        IF(B==1) THEN ; O=2 ; ELSE ; O=1 ; ENDIF
        IF(MODEL_TYPE==1) THEN ! diffusion
            A=AGE(b)%LAST-1
            TEMP(b)=DBLE(age(b)%LAST)*(N(b,A,Y)*(1.0-T(b,a))+N(o,A,Y)*T(o,a))*EXP(-M(b,A,y)-F(b,a,y))+ &
                    ((PLUSAGE(b)+1.0)*N(b,A+1,Y)*(1.0-T(b,a+1))+(PLUSAGE(o)+1.0)*N(o,A+1,Y)*T(o,a+1))*EXP(-M(b,A+1,y)-F(b,a+1,y))
          ELSE ! overlap
            TEMP(b)=DBLE(age(b)%LAST)*N(b,age(b)%last,Y+1)+(PLUSAGE(b)+1.0)*N(b,age(b)%LAST+1,Y+1)
        ENDIF
      END DO
      DO B=1,N_BOX
        N(b,age(b)%last,Y+1)=N(b,age(b)%last,Y+1)+N(b,AGE(b)%last+1,Y+1)
        PLUSAGE(b)=MAX(REAL(AGE(b)%LAST),TEMP(b)/N(b,AGE(b)%last,Y+1))
      END DO

!     begin a new cohort (recruits)
      ! Deterministic part
      IF(AGE(1)%FIRST==0) THEN ! Get SSB at beginning of next forecasted year (age 0 produced at beginning of year)
        DO B=1,N_BOX ; N(b,AGE(1)%FIRST,Y+1)=0; END DO
        CALL GETWAA(Y+1,LOOP); CALL GETSSB(y+1)
      ENDIF
      XMU=0
      DO B=1,N_BOX
        IF(Y==YEAR%LAST) THEN
          IF(R_FILE%FILETYPE==0 .AND. ESTIMATE_AUTOCORRELATION<=0) THEN ;  ETA(b)=sr(6,b)  ! Starting residual from single entry in ASC-ii file
            ELSE; ETA(b)=sr(5+Y,b)                                                          ! Starting residual from binary file entry for year before patch or estimated internally
          ENDIF
        ENDIF
        CALL RECRUITS(rec,b,Y+1)
        N(b,AGE(b)%first,Y+1)=MAX(0.0d0,REC)
        XMU=XMU+MAX(0.0d0,REC)
        IF(LOOP==1) N(b,AGE(b)%first,Y+1)=N(b,AGE(b)%first,Y+1)*Rec_Mod(b,Y+1)
      END DO ! b
      ! Stochastic part
      DO B=1,N_BOX
        IF(SR_BOX>0 .AND. SR_TYPE(b)/=6) THEN ! sex specific
            REC=XMU ;
            IF(SR_BOX>2) THEN ; iBOX=1 ; ELSE ; iBOX=SR_BOX ; ENDIF  ! s/r relationship assuming sex S is the limiting factor
            IF(B>1) EXIT      ! one s/r relationship governs both sexes, so only need one pass through loop
          ELSE ! area/stock specific or SR_TYPE=6
            REC=N(B,age(b)%first,Y+1) ; iBOX=B
        ENDIF
        RHO=sr(5,ibox) ; STDDEV=ABS(SR(4,ibox))
        IF(REC>0 .AND. LOOP>1) THEN
            IF(SR_TYPE(b)==6) THEN
                DRAW=REC_YEAR(1)+INT(RAND(LOOP,b,Y)*N_REC)   ! draw at random from historical recruitments
                REC=N(b,age(ibox)%first,DRAW)
              ELSE
                ETA(ibox)=RHO*ETA(ibox)+STDDEV*U(LOOP,iBOX,y+1) ! STDDEV is std dev. of random component of rec. variation
                IF(ERRORTYPE==1) THEN
                    ! mu = E[ln(R)] = ln(E[R])-mse/2 where E is expectation
                    REC=REC*EXP(ETA(iBOX))    ! stochastic recruitment value with mean = deterministic value
                  ELSE
                    REC=REC+ETA(iBOX)
                ENDIF
            ENDIF
            N(B,age(b)%first,Y+1)=REC*Rec_Mod(B,Y+1)
          ELSEIF(REC>0 .AND. LOOP==1 .AND. ERRORTYPE==1) THEN
            N(B,age(b)%first,Y+1)=REC*BIAS_ADJ(b)  ! bias-correct MLE from lognormal distribution
        ENDIF
      END DO
      IF(SR_BOX>0.AND.SR_TYPE(BOX%STORE)/=6) THEN
        DO B=1,N_BOX ; N(B,age(b)%first,Y+1)=REC*SEX_FRACTION(B)*Rec_Mod(B,Y+1) ; END DO
      ENDIF !sex

      RETURN ; END

!-----------------------------------------------------------------------------
      SUBROUTINE RECRUITS(rec,b,y)
!     deterministic recruitment for next year. Note that recruitement is tied to the stock/recruit parameters
!     for SR_BOX (which should be set to zero if the analysis is not sex specific, but area/stock specific)
!-----------------------------------------------------------------------------
      USE CNTRL ; USE STATISTICS ; IMPLICIT NONE
      INTEGER :: Y,B,S
      REAL (KIND=8) :: ssbtemp,rec

      SSBTEMP=0

      IF(SR_BOX<=0) THEN ! this analysis is area/stock specific, not sex specific
          SSBTEMP=SSB(b,Y-AGE(b)%FIRST) ; BOX%STORE=B
        ELSE             ! sex-specific
          DO S=1,N_BOX ; IF(SR_BOX==S .OR. SR_BOX>2) SSBTEMP=SSBTEMP+SSB(S,Y-AGE(S)%FIRST) ;  END DO
          IF(SR_BOX<=2) THEN ; BOX%STORE=SR_BOX ; ELSE ; BOX%STORE=1 ; ENDIF ! Use parameters specified for sex 1 if SR_BOX>2
      ENDIF
      IF(SR_TYPE(BOX%STORE)==1) THEN ! {Beverton and Holt curve}
          REC=sr(2,BOX%STORE)*ssbtemp/(ssbtemp+sr(3,BOX%STORE))
        ELSEIF(SR_TYPE(BOX%STORE)==2) THEN ! {Ricker curve}
          REC=sr(2,BOX%STORE)*ssbtemp*EXP(-ssbtemp*sr(3,BOX%STORE))
        ELSEIF(SR_TYPE(BOX%STORE)>2 .and. SR_TYPE(BOX%STORE)<6) THEN ! {Two-line curves}
          IF(SSBTEMP<SR(3,BOX%STORE)) THEN ; REC=SR(2,BOX%STORE)*SSBTEMP/SR(3,BOX%STORE) ; ELSE ; REC=SR(2,BOX%STORE) ; ENDIF
        ELSE ! {constant recruitment at value of SR(2)}
          REC=SR(2,BOX%STORE)
      ENDIF
      IF(SR_BOX>0) THEN ! Split recruits by sex
        REC=REC*SEX_FRACTION(b)
      ENDIF

      RETURN ; END

!-----------------------------------------------------------------------------
      SUBROUTINE SR_FIT(loop)
      ! Fits spawner recruit relationships
!-----------------------------------------------------------------------------
      USE CNTRL ; USE STATISTICS ; IMPLICIT NONE
      INTEGER :: B,Y,COUNTS,NPARM,MAX_ITER=100,N_CYCLES,loop
      REAL (KIND=8) :: maxssb,minssb,aveSSB,maxR,sumR,logsumr,PARM(100),p(100,99),x(100),O(100),PB(100),FUNC,REC,SSBTEMP,REGRESS,ss
      EXTERNAL OBJECTIVE_SR,OBJECTIVE_AUTO
      DO Y=1,YEAR%LAST
        CALL GETWAA(Y,1) ; CALL GETSSB(Y)
      END DO
      IF(SR_BOX>0) THEN ; N_CYCLES=1 ; ELSE ; N_CYCLES=N_BOX ; ENDIF
      DO B=1,N_CYCLES
        maxR=-1 ; maxSSB=-1 ; minSSB=1.0D+32 ; aveSSB=0.0 ; sumR=0.0 ; logsumr=0.0 ; var_R(b)=0.0 ; COUNTS=0 ; regress=0
        DO Y=REC_YEAR(1),REC_YEAR(2)
          IF (Y>=YEAR%FIRST+age(b)%FIRST) THEN
            IF(SR_BOX>0) THEN ; REC=N(1,age(1)%FIRST,Y)+N(2,age(2)%FIRST,Y) ; ELSE ; REC=N(B,age(b)%FIRST,Y) ; ENDIF
            maxR=MAX(maxR,REC) ; sumR=sumR+REC ; logsumr=logsumr+LOG(REC) ; COUNTS=COUNTS+1
            SELECT CASE (SR_BOX)
              CASE (1) ; SSBTEMP=SSB(1,Y-age(b)%FIRST)
              CASE (2) ; SSBTEMP=SSB(2,Y-age(b)%FIRST)
              CASE (3) ; SSBTEMP=SSB(1,Y-age(b)%FIRST)+SSB(2,Y-age(b)%FIRST)
              CASE DEFAULT ; SSBTEMP=SSB(B,Y-age(b)%FIRST)
            END SELECT
            IF(ERRORTYPE==1) THEN
                var_R(b)=var_R(b)+LOG(REC)**2
                regress=regress+LOG(REC/SSBTEMP)
              ELSE
                var_R(b)=var_R(b)+REC**2
                regress=regress+REC*SSBTEMP ; ss=ss+SSBTEMP*SSBTEMP
            ENDIF
          END IF
        END DO
        ! compute variance of observed recruitments using formula that avoids rounding error
        IF(ERRORTYPE==1) THEN
            var_R(b)=var_R(b)/DBLE(COUNTS)-(logsumR/DBLE(COUNTS))**2
          ELSE
            var_R(b)=var_R(b)/DBLE(COUNTS)-(sumR/DBLE(COUNTS))**2
        ENDIF
        DO Y=SSB_YEAR(1),SSB_YEAR(2)
          SELECT CASE (SR_BOX)
            CASE (1) ; SSBTEMP=SSB(1,Y)
            CASE (2) ; SSBTEMP=SSB(2,Y)
            CASE (3) ; SSBTEMP=SSB(1,Y)+SSB(2,Y)
            CASE DEFAULT ; SSBTEMP=SSB(B,Y)
          END SELECT
          maxSSB=MAX(maxSSB,SSBTEMP)
          minSSB = MIN(minSSB,SSBTEMP)
          aveSSB = aveSSB + SSBTEMP/dble(ssb_year(2)-ssb_year(1)+1)
        END DO
        IF(SR_TYPE(b)==1) THEN ! {Beverton and Holt curve}
            IF(loop==1) THEN; SR_GUESS(2,b)=0.6*MAXR; SR_GUESS(3,b)=maxSSB/10.0; ENDIF
            NPARM=2 ; PARM(1)=SR_GUESS(2,b) ; PARM(2)=SR_GUESS(3,b) ; BOX%STORE=B
            IF(ESTIMATE_AUTOCORRELATION>0) THEN; NPARM=3; PARM(3)=SR(5,b); IF(ABS(PARM(3))<0.01) PARM(3)=0.01; ENDIF
            CALL ESTIMATE(OBJECTIVE_SR,PARM,P,X,O,PB,NPARM,SEED,MAX_ITER)
            CALL OBJECTIVE_SR(PARM,func,NPARM)
            SR(2,b)=PARM(1) ; SR(3,b)=PARM(2) ; IF(ESTIMATE_AUTOCORRELATION>0) THEN; SR(5,b)=PARM(3); ENDIF
            IF(LOOP==1) THEN; SR_GUESS(2,b)=PARM(1); SR_GUESS(3,b)=PARM(2); ENDIF
          ELSEIF(SR_TYPE(b)==2) THEN ! {Ricker curve}
            NPARM=2 ; PARM(1)=0.9*maxR*exp(1.)/aveSSB ; PARM(2)=1.0/aveSSB ; BOX%STORE=B
            IF(ESTIMATE_AUTOCORRELATION>0) THEN; NPARM=3; PARM(3)=SR(5,b); IF(ABS(PARM(3))<0.01) PARM(3)=0.01; ENDIF
            CALL ESTIMATE(OBJECTIVE_SR,PARM,P,X,O,PB,NPARM,SEED,MAX_ITER)
            CALL OBJECTIVE_SR(PARM,func,NPARM)
            SR(2,b)=PARM(1) ; SR(3,b)=PARM(2); IF(ESTIMATE_AUTOCORRELATION>0) SR(5,b)=PARM(3)
          ELSEIF(SR_TYPE(b)==3) THEN ! {Two-line curve based on linear regression with inflection at maximum SSB observed}
            IF(ERRORTYPE==1) THEN ! maximum recruitment = a*maxssb
                SR(2,b)=EXP(regress/DBLE(counts))*maxssb
              ELSE
                SR(2,b)=REGRESS/SS*maxssb
            ENDIF
            SR(3,b)=MAXSSB    ! maximum ssb observed
          ELSEIF(SR_TYPE(b)==4.OR.SR_TYPE(b)==5) THEN ! {Two-line curve based on geometric mean recruitment}
            IF(ERRORTYPE==1) THEN ! maximum recruitment = a*maxssb
                SR(2,b)=EXP(logsumr/DBLE(counts)) ! geometric mean recruitment
              ELSE
                SR(2,b)=sumR/dble(counts)         ! arithmetic mean recruitment
            ENDIF
            IF(SR_TYPE(b)==4) THEN
                SR(3,b)=minSSB ! Inflection at minimum observed SSB
              ELSE
                SR(3,b)=aveSSB ! inflection at average SSB
            ENDIF
          ELSE ! {constant average recruitment for option 6: note that this ignores ERROR_TYPE}
            SR(2,b)=sumR/dble(counts)
            SR(4,b)=MAX(SQRT(var_R(b)*DBLE(COUNTS)/(DBLE(COUNTS)-1)),SR_mse_minimum(B))
        ENDIF
        IF(SR_TYPE(b)>2 .AND. SR_TYPE(b)<6) THEN ! Estimate standard error and autocorrelation for two-line models
            SR(4,b)=0.0
            DO Y=REC_YEAR(1),REC_YEAR(2)
              IF (Y>=YEAR%FIRST+age(b)%FIRST) THEN
                IF(SR_BOX>0) THEN ; REC=N(1,age(1)%FIRST,Y)+N(2,age(2)%FIRST,Y) ; ELSE ; REC=N(B,age(b)%FIRST,Y) ; ENDIF
                SELECT CASE (SR_BOX)
                  CASE (1) ; SSBTEMP=SSB(1,Y-age(b)%FIRST)
                  CASE (2) ; SSBTEMP=SSB(2,Y-age(b)%FIRST)
                  CASE (3) ; SSBTEMP=SSB(1,Y-age(b)%FIRST)+SSB(2,Y-age(b)%FIRST)
                  CASE DEFAULT ; SSBTEMP=SSB(B,Y-age(b)%FIRST)
                END SELECT
                IF(SSBTEMP<SR(3,b)) THEN ; O(1)=SR(2,B)*SSBTEMP/SR(3,b) ; ELSE ; O(1)=SR(2,b) ; ENDIF
                IF(ERRORTYPE==1) THEN
                    SR(4,B)=SR(4,B)+(LOG(REC/O(1)))**2
                  ELSE
                    SR(4,B)=SR(4,B)+(REC-O(1))**2
                ENDIF
              END IF
            END DO
            SR(4,B)=SR(4,B)/(COUNTS-1) ! note: Objective_Auto recomputes the variance when called
            IF(ESTIMATE_AUTOCORRELATION>0) THEN;
              NPARM=1; PARM(1)=SR(5,b); IF(ABS(PARM(1))<0.01) PARM(1)=0.1; BOX%STORE=B
              CALL ESTIMATE(OBJECTIVE_AUTO,PARM,P,X,O,PB,NPARM,SEED,MAX_ITER)
              CALL OBJECTIVE_AUTO(PARM,func,NPARM)
              SR(5,b)=PARM(1)
            ENDIF
            SR(4,B)=MAX( SQRT(SR(4,B)),SR_mse_minimum(B) )  ! Impose a minimum variance
        ENDIF
      END DO
      IF(SR_BOX>0) THEN ; SR(2,2)=SR(2,1) ; SR(3,2)=SR(3,1) ; SR(4,2)=SR(4,1) ; SR(5,2)=SR(5,1) ; ENDIF

      RETURN ; END

!-----------------------------------------------------------------------------
      SUBROUTINE GETYIELD(y)
! Computes the yield by area
!-----------------------------------------------------------------------------
      USE CNTRL ; USE STATISTICS ; IMPLICIT NONE
      INTEGER :: B,A,Y,O,G
      REAL (KIND=8) :: z,temp
      DO B=1,N_BOX
        IF(B==1) THEN ; O=2 ; ELSE ; O=1 ; ENDIF ; SUM_YLD(B,Y)=0.0
        DO G=1,N_GEARS(B)
          YLD(B,G,Y)=0.0
          DO A=age(b)%FIRST, age(b)%LAST
            Z=M(b,A,y)+F(b,A,Y)
            temp=(N(b,A,Y)*(1.0-T(b,a))+N(o,a,y)*T(o,a))*Fp(b,g,A,Y)*W(b,g,A,Y)*(1.-EXP(-z))/z
            YLD(b,g,Y)=YLD(b,g,Y)+temp
          END DO
          SUM_YLD(B,Y)=SUM_YLD(B,Y)+YLD(B,G,Y)
        END DO
        IF(SUM_YLD(B,Y)<0.1) SUM_YLD(B,Y)=0.0
      END DO
     
      RETURN ; END

!-----------------------------------------------------------------------------
      SUBROUTINE GETYR(LOOP)
!-----------------------------------------------------------------------------
      USE CNTRL ; USE STATISTICS ; IMPLICIT NONE
      INTEGER :: finish(4),nparm,LOOP,MAX_ITER=100,B,S,oth
      REAL (KIND=8) :: PARM(100),p(100,99),x(100),O(100),PB(100),rmax(2),func,oldF,slope(2),RATIO,VIRGIN(2),TEMP,EQUILIBRIUM_REC
      EXTERNAL OBJECTIVE_YR
      FINISH=0

      IF(LOOP==1) WRITE(*,*) 'COMPUTING YIELD PER RECRUIT STATS'
      IF(SR_BOX>0) THEN ; NPARM=1 ; ELSE ; NPARM=N_BOX ; ENDIF ; PARM=0.4/DBLE(nparm)
      IF(ABS(F_MULT)>-999) NPARM=1

!     Fmax: maximum yield per recruit
      PR_TYPE=1
      CALL OBJECTIVE_YR(PARM,func,NPARM)
      CALL ESTIMATE(OBJECTIVE_YR,PARM,P,X,O,PB,NPARM,SEED,MAX_ITER)
      CALL OBJECTIVE_YR(PARM,func,NPARM)
      DO B=1,NPARM
        BENCHMARKS(LOOP,B)%YPRMAX=YRTEMP(b) ; BENCHMARKS(LOOP,B)%SPRMAX=SRTEMP(b) ; BENCHMARKS(LOOP,B)%FMAX=PARM(b)
      END DO
      IF(SR_BOX>0 .OR. F_MULT>-999) THEN ! Sex-specific or F-linked, so only one pass made
        BENCHMARKS(LOOP,2)%FMAX=BENCHMARKS(LOOP,1)%FMAX ; BENCHMARKS(LOOP,2)%YPRMAX=YRTEMP(2) ; BENCHMARKS(LOOP,2)%SPRMAX=SRTEMP(2)
        IF(F_MULT<0) THEN; BENCHMARKS(LOOP,2)%FMAX=ABS(F_MULT)*BENCHMARKS(LOOP,1)%FMAX
          ELSE; BENCHMARKS(LOOP,2)%FMAX=F_MULT*F_CURRENT_RATIO*BENCHMARKS(LOOP,1)%FMAX
        ENDIF
      ENDIF

!     F0.1: (F where slope = 0.1 slope at F=0) -- for two-areas with independent (unlinked) F's this is obtained assuming the other area is fished at Fmax.
      ! Forward difference to find slope at origin = (y/r{x}-0)/x
      IF(LOOP==1) WRITE(*,*) '     F0.1'
      DO B=1,NPARM; X(b)=0.000001; END DO
      CALL OBJECTIVE_YR(X,func,NPARM)
      SLOPE(1)=0.1*ABS(FUNC)/0.000001 !numerical approx. to slope of y/r curve at origin

      ! brute search for F that gives 10% slope at origin
      DO S=1,NPARM ; X(S)=BENCHMARKS(LOOP,S)%FMAX ; END DO ; CALL OBJECTIVE_YR(X,FUNC,NPARM)
      OLDF=ABS(FUNC)
      DO WHILE (FINISH(1)==0)
        IF(X(1)<0.001.or.X(NPARM)<0.001) THEN
            EXIT
          ELSE
            DO S=1,NPARM; X(S)=X(S)-0.001 ; END DO
            CALL OBJECTIVE_YR(X,FUNC,NPARM)
        ENDIF
        IF((OLDF-ABS(FUNC))/0.001>SLOPE(1)) EXIT
        OLDF=ABS(FUNC)
      END DO
      DO B=1,NPARM
        BENCHMARKS(LOOP,B)%F01=X(B) ; BENCHMARKS(LOOP,B)%YPR01=YRTEMP(B) ; BENCHMARKS(LOOP,B)%SPR01=SRTEMP(B)
      END DO
      IF(SR_BOX>0 .OR. F_MULT>-999) THEN ! Sex-specific or F-linked, so only one pass made
        BENCHMARKS(LOOP,2)%F01=BENCHMARKS(LOOP,1)%F01 ; BENCHMARKS(LOOP,2)%YPR01=YRTEMP(2) ; BENCHMARKS(LOOP,2)%SPR01=SRTEMP(2)
        IF(F_MULT<0) THEN; BENCHMARKS(LOOP,2)%F01=ABS(F_MULT)*BENCHMARKS(LOOP,1)%F01
          ELSE; BENCHMARKS(LOOP,2)%F01=F_MULT*F_CURRENT_RATIO*BENCHMARKS(LOOP,1)%F01
        ENDIF
      ENDIF

!     F at 90% max Y/R (w/mixing, many combinations of F(1) and F(2) will give 90% maxY/R, decrementing each F equally is one way)
      IF(LOOP==1) WRITE(*,*) '     F AT 90% max Y/R'
      DO B=1,NPARM ; X(B)=BENCHMARKS(LOOP,B)%FMAX ; END DO
      SELECT CASE (YR_BOX)
        CASE (1) ;     TEMP=BENCHMARKS(LOOP,1)%YPRMAX
        CASE (2) ;     TEMP=BENCHMARKS(LOOP,2)%YPRMAX
        CASE DEFAULT ; TEMP=BENCHMARKS(LOOP,1)%YPRMAX+BENCHMARKS(LOOP,2)%YPRMAX
      END SELECT
      DO WHILE (FINISH(1)==0)
        DO B=1,NPARM ; X(b)=X(b)-0.001 ; END DO
        CALL OBJECTIVE_YR(X,FUNC,NPARM)
        IF( (X(1)<0.001.OR.X(NPARM)<0.001) .OR. ABS(FUNC)<=0.9*TEMP) EXIT
      END DO
      DO B=1,NPARM
        BENCHMARKS(LOOP,b)%F90=X(b) ; BENCHMARKS(LOOP,b)%YPRf90=YRTEMP(b) ; BENCHMARKS(LOOP,b)%SPRf90=SRTEMP(b)
      END DO
      IF(SR_BOX>0 .OR. F_MULT>-999) THEN ! Sex-specific or F-linked, so only one pass made
        BENCHMARKS(LOOP,2)%F90=BENCHMARKS(LOOP,1)%F90; BENCHMARKS(LOOP,2)%YPRf90=YRTEMP(2) ; BENCHMARKS(LOOP,2)%SPRf90=SRTEMP(2)
        IF(F_MULT<0) THEN; BENCHMARKS(LOOP,2)%F90=ABS(F_MULT)*BENCHMARKS(LOOP,1)%F90
          ELSE; BENCHMARKS(LOOP,2)%F90=F_MULT*F_CURRENT_RATIO*BENCHMARKS(LOOP,1)%F90
        ENDIF
      ENDIF

!     F at 75% Fmax
      IF(LOOP==1) WRITE(*,*) '     F AT 75% F of max Y/R'
      DO B=1,NPARM ; X(B)=0.75*BENCHMARKS(LOOP,B)%FMAX ; END DO
      CALL OBJECTIVE_YR(X,FUNC,NPARM)
      DO B=1,NPARM
        BENCHMARKS(LOOP,B)%F75Fmax=X(B) ; BENCHMARKS(LOOP,B)%YPR75Fmax=YRTEMP(B) ; BENCHMARKS(LOOP,B)%SPR75Fmax=SRTEMP(B)
      END DO
      IF(SR_BOX>0 .OR. F_MULT>-999) THEN ! Sex-specific or F-linked, so only one pass made
        BENCHMARKS(LOOP,2)%YPR75Fmax=YRTEMP(2) ; BENCHMARKS(LOOP,2)%SPR75Fmax=SRTEMP(2)
        BENCHMARKS(LOOP,2)%F75Fmax= BENCHMARKS(LOOP,1)%F75Fmax 
        IF(F_MULT<0) THEN; BENCHMARKS(LOOP,2)%F75Fmax=ABS(F_MULT)*BENCHMARKS(LOOP,1)%F75Fmax
          ELSE; BENCHMARKS(LOOP,2)%F75Fmax=F_MULT*F_CURRENT_RATIO*BENCHMARKS(LOOP,1)%F75Fmax
        ENDIF
      ENDIF

!     F where spawning potential ratio is one of the three input values
      IF(LOOP==1) WRITE(*,*) '     F AT X% SPR'
      X=0  ; CALL OBJECTIVE_YR(X,func,NPARM)
      DO B=1,N_BOX ; VIRGIN(B)= SRTEMP(B) ; END DO ! spawner per recruit at F=0

      DO B=1,NPARM
         DO S=1,NPARM ; X(S)=3.0 ; END DO
         DO WHILE (FINISH(2)==0)
           DO S=1,NPARM ; X(S)=X(S)-0.001 ; END DO ; CALL OBJECTIVE_YR(X,FUNC,NPARM)
           SELECT CASE (YR_BOX)
               CASE (1) ;     RATIO=SRTEMP(1)/VIRGIN(1)                           ! only stock 1
               CASE (2) ;     RATIO=SRTEMP(2)/VIRGIN(2)                           ! only stock 2
               CASE DEFAULT ; RATIO=(SRTEMP(1)+SRTEMP(2))/(VIRGIN(1)+VIRGIN(2))   ! both stocks
           END SELECT
           SELECT CASE (SR_BOX)
               CASE (1) ;     RATIO=SRTEMP(1)/VIRGIN(1)                           ! only sex 1
               CASE (2) ;     RATIO=SRTEMP(2)/VIRGIN(2)                           ! only sex 2
               CASE (3) ;     RATIO=(SRTEMP(1)+SRTEMP(2))/(VIRGIN(1)+VIRGIN(2))   ! both sexes
               CASE DEFAULT                                                       ! not sex-specific
           END SELECT
           IF(RATIO>=SPR_TARGET(1) .AND. FINISH(4)==0) THEN
               BENCHMARKS(LOOP,b)%F20=X(b) ; BENCHMARKS(LOOP,b)%SPR20=SRTEMP(b) ; BENCHMARKS(LOOP,b)%YPR20=yrtemp(b) ; FINISH(4)=1
               IF(SR_BOX>0 .OR. F_MULT>-999) THEN ; BENCHMARKS(LOOP,2)%SPR20=SRTEMP(2) ; BENCHMARKS(LOOP,2)%YPR20=yrtemp(2) ; ENDIF
             ELSEIF(RATIO>=SPR_TARGET(2) .AND. FINISH(3)==0) THEN
               BENCHMARKS(LOOP,b)%F30=X(b) ; BENCHMARKS(LOOP,b)%SPR30=SRTEMP(b) ; BENCHMARKS(LOOP,b)%YPR30=yrtemp(b) ; FINISH(3)=1
               IF(SR_BOX>0 .OR. F_MULT>-999) THEN ; BENCHMARKS(LOOP,2)%SPR30=SRTEMP(2) ; BENCHMARKS(LOOP,2)%YPR30=yrtemp(2) ; ENDIF
             ELSEIF(RATIO>=SPR_TARGET(3) .AND. FINISH(2)==0) THEN
               BENCHMARKS(LOOP,b)%F40=X(b) ; BENCHMARKS(LOOP,b)%SPR40=SRTEMP(b) ; BENCHMARKS(LOOP,b)%YPR40=yrtemp(b) ; FINISH(2)=1
               IF(SR_BOX>0 .OR. F_MULT>-999) THEN ; BENCHMARKS(LOOP,2)%SPR40=SRTEMP(2) ; BENCHMARKS(LOOP,2)%YPR40=yrtemp(2) ; ENDIF
           ENDIF
         END DO
         FINISH=0
      END DO
      IF(SR_BOX>0 .OR. F_MULT>-999) THEN ! Sex-specific or F-linked
        BENCHMARKS(LOOP,2)%F20=BENCHMARKS(LOOP,1)%F20
        BENCHMARKS(LOOP,2)%F30=BENCHMARKS(LOOP,1)%F30
        BENCHMARKS(LOOP,2)%F40=BENCHMARKS(LOOP,1)%F40
        IF(F_MULT<0) THEN;
            BENCHMARKS(LOOP,2)%F20=ABS(F_MULT)*BENCHMARKS(LOOP,1)%F20
            BENCHMARKS(LOOP,2)%F30=ABS(F_MULT)*BENCHMARKS(LOOP,1)%F30
            BENCHMARKS(LOOP,2)%F40=ABS(F_MULT)*BENCHMARKS(LOOP,1)%F40
          ELSE;
            BENCHMARKS(LOOP,2)%F20=F_MULT*F_CURRENT_RATIO*BENCHMARKS(LOOP,1)%F20
            BENCHMARKS(LOOP,2)%F30=F_MULT*F_CURRENT_RATIO*BENCHMARKS(LOOP,1)%F30
            BENCHMARKS(LOOP,2)%F40=F_MULT*F_CURRENT_RATIO*BENCHMARKS(LOOP,1)%F40
        ENDIF
      ENDIF

!     msy calculations directly from spawner-recruit and yield-per-recruit statistics
      PR_TYPE=2
      IF(LOOP==1) WRITE(*,*) '     F at MSY'
      IF(SR_TYPE(1)==1 .OR. SR_TYPE(1)==2) THEN
          ! (1) Beverton and Holt sr curve or (2) Ricker curve
          CALL ESTIMATE(OBJECTIVE_YR,PARM,P,X,O,PB,NPARM,SEED,MAX_ITER)
          CALL OBJECTIVE_YR(PARM,func,NPARM)
          DO B=1,NPARM
            BENCHMARKS(LOOP,B)%FMSY=PARM(B)
            BENCHMARKS(LOOP,B)%MSY=MSYTEMP(B)
            BENCHMARKS(LOOP,B)%SPRMSY=SRTEMP(B) ; BENCHMARKS(LOOP,B)%SSBMSY=SSBMSYTEMP(B)
          END DO
          IF(SR_BOX>0 .OR. F_MULT>-999) THEN ! Sex-specific or F-linked, so only one pass made
            BENCHMARKS(LOOP,2)%MSY=MSYTEMP(2) ; BENCHMARKS(LOOP,2)%SPRMSY=SRTEMP(2);
            BENCHMARKS(LOOP,2)%FMSY=BENCHMARKS(LOOP,1)%FMSY ;  BENCHMARKS(LOOP,2)%SSBMSY=SSBMSYTEMP(2)
            IF(F_MULT<0) THEN; BENCHMARKS(LOOP,2)%FMSY=ABS(F_MULT)*BENCHMARKS(LOOP,1)%FMSY
              ELSE; BENCHMARKS(LOOP,2)%FMSY=F_MULT*F_CURRENT_RATIO*BENCHMARKS(LOOP,1)%FMSY
            ENDIF
          ENDIF

        ELSEIF(SR_TYPE(1)>2 .OR. SR_TYPE(1)<6) THEN
          ! Two-line linear regression sr curve
          ! by definition the minimum feasible spr is 1/(slope of two-line method).  Note
          ! that the 'MSY' is not actually sustainable unless the optimal S/R ratio is at or above
          ! S(inflection)/Rmax. Therefore, if Fmax produces SPR below S(inflection)/Rmax, the MSY level
          ! is set equivalent to the SPR level of {S(inflection)/Rmax}/{S(virgin)/Rmax}. This can be tested
          ! by making long-term projections with F values progressively lower than Fmax. Note that longterm
          ! stochastic projections at Fmsy defined as above actually tend to fall below the MSY level
          ! owing to the discontinuity at the hinge, therefore the target SPR should be slightly higher than
          ! {S(inflection)/Rmax}/{S(virgin)/Rmax}. A 1% increase in target SPR results in nearly identical
          ! longterm yields and allows stochastic projections at Fmsy to average nearer to the deterministic value.
          FINISH=0
          DO B=1,N_BOX
            RMAX(b)=SR(2,B)*BIAS_ADJ(B)
            IF(BENCHMARKS(LOOP,B)%sprmax>(1010.0*SR(3,b)/SR(2,b)) ) THEN  ! Fmsy equivalent to Fmax
                BENCHMARKS(LOOP,B)%FMSY=BENCHMARKS(LOOP,B)%Fmax; BENCHMARKS(LOOP,B)%MSY=BENCHMARKS(LOOP,B)%Yprmax*rmax(b)/1000.0
                BENCHMARKS(LOOP,B)%sprmSY=BENCHMARKS(LOOP,B)%sprmax;
                BENCHMARKS(LOOP,B)%SSBMSY=BENCHMARKS(LOOP,B)%sprmax*rmax(b)/1000.0
                FINISH(B)=1;
            ENDIF
          END DO
          IF(SUM(FINISH)<N_BOX) THEN ! Fmax not sustainable (so back off a little on F)
              FINISH=0
              DO S=1,NPARM ; X(S)=BENCHMARKS(LOOP,S)%FMAX ; END DO
              DO WHILE (SUM(FINISH)<N_BOX)
                DO S=1,NPARM ; X(S)=X(S)-0.001 ; END DO ; CALL OBJECTIVE_YR(X,FUNC,NPARM)
                DO B=1,N_BOX
                  SELECT CASE (SR_BOX)
                      CASE (1) ;     RATIO=SRTEMP(1)               ! only sex 1
                      CASE (2) ;     RATIO=SRTEMP(2)               ! only sex 2
                      CASE (3) ;     RATIO=(SRTEMP(1)+SRTEMP(2))   ! both sexes
                      CASE DEFAULT ; RATIO=SRTEMP(B)               ! not sex-specific
                  END SELECT
                  IF(RATIO>1010.0*SR(3,b)/SR(2,b) .AND. SUM(FINISH)<N_BOX) THEN
                    ! find F that gives SPR value 10% higher than slope of the descending limb (multiply by 1100 instead of 1000)
                    ! to make sure ratio sufficiently above threshold
                    BENCHMARKS(LOOP,b)%FMSY=X(b) ; BENCHMARKS(LOOP,b)%SPRMSY=SRTEMP(b) ;
                    BENCHMARKS(LOOP,B)%MSY=yrtemp(b)*rmax(b)/1000.0
                    BENCHMARKS(LOOP,B)%SSBMSY=BENCHMARKS(LOOP,B)%SPRMSY*rmax(b)/1000.0
                    FINISH(b)=1
                    IF(SR_BOX>0 .OR. F_MULT>-999) THEN
                      BENCHMARKS(LOOP,2)%SPRMSY=SRTEMP(2) ; BENCHMARKS(LOOP,2)%MSY=BENCHMARKS(LOOP,1)%MSY
                      IF(F_MULT<0) THEN; BENCHMARKS(LOOP,2)%FMSY=ABS(F_MULT)*BENCHMARKS(LOOP,1)%FMSY
                        ELSE; BENCHMARKS(LOOP,2)%FMSY=F_MULT*F_CURRENT_RATIO*BENCHMARKS(LOOP,1)%FMSY
                      ENDIF
                    ENDIF
                  ENDIF
                END DO
              END DO
          ENDIF

        ELSE ! constant recruitment (equivalent to resampling, so no bias correction)
          DO B=1,N_BOX
            rmax(b)=SR(2,B)
            IF(SR_BOX>0) RMAX(B)=RMAX(B)*SEX_FRACTION(b)
            BENCHMARKS(LOOP,B)%FMSY=BENCHMARKS(LOOP,B)%Fmax ; BENCHMARKS(LOOP,B)%MSY=BENCHMARKS(LOOP,B)%Yprmax*rmax(b)/1000.0
            BENCHMARKS(LOOP,B)%sprmSY=BENCHMARKS(LOOP,B)%sprmax ; BENCHMARKS(LOOP,B)%SSBMSY=BENCHMARKS(LOOP,B)%sprmax*rmax(b)/1000.0
          END DO
      ENDIF


      DO B=1,N_BOX
        IF(N_BOX>1) THEN ; Oth=2 ; ELSE ; Oth=1 ; ENDIF
        BENCHMARKS(LOOP,B)%SSBMAX=BENCHMARKS(LOOP,B)%sprmax*EQUILIBRIUM_REC(B,BENCHMARKS(LOOP,B)%sprmax,BENCHMARKS(LOOP,OTH)%sprmax)
        BENCHMARKS(LOOP,B)%YMAX  =BENCHMARKS(LOOP,B)%YPRMAX*EQUILIBRIUM_REC(B,BENCHMARKS(LOOP,B)%sprmax,BENCHMARKS(LOOP,OTH)%sprmax)
        BENCHMARKS(LOOP,B)%SSBf01=BENCHMARKS(LOOP,B)%spr01* EQUILIBRIUM_REC(B,BENCHMARKS(LOOP,B)%spr01,BENCHMARKS(LOOP,OTH)%spr01)
        BENCHMARKS(LOOP,B)%Y01   =BENCHMARKS(LOOP,B)%YPR01* EQUILIBRIUM_REC(B,BENCHMARKS(LOOP,B)%spr01,BENCHMARKS(LOOP,OTH)%spr01)
        BENCHMARKS(LOOP,B)%SSB20 =BENCHMARKS(LOOP,B)%spr20* EQUILIBRIUM_REC(B,BENCHMARKS(LOOP,B)%spr20,BENCHMARKS(LOOP,OTH)%spr20)
        BENCHMARKS(LOOP,B)%Y20   =BENCHMARKS(LOOP,B)%YPR20* EQUILIBRIUM_REC(B,BENCHMARKS(LOOP,B)%spr20,BENCHMARKS(LOOP,OTH)%spr20)
        BENCHMARKS(LOOP,B)%SSB30 =BENCHMARKS(LOOP,B)%spr30* EQUILIBRIUM_REC(B,BENCHMARKS(LOOP,B)%spr30,BENCHMARKS(LOOP,OTH)%spr30)
        BENCHMARKS(LOOP,B)%Y30   =BENCHMARKS(LOOP,B)%YPR30* EQUILIBRIUM_REC(B,BENCHMARKS(LOOP,B)%spr30,BENCHMARKS(LOOP,OTH)%spr30)
        BENCHMARKS(LOOP,B)%SSB40 =BENCHMARKS(LOOP,B)%spr40* EQUILIBRIUM_REC(B,BENCHMARKS(LOOP,B)%spr40,BENCHMARKS(LOOP,OTH)%spr40)
        BENCHMARKS(LOOP,B)%Y40   =BENCHMARKS(LOOP,B)%YPR40* EQUILIBRIUM_REC(B,BENCHMARKS(LOOP,B)%spr40,BENCHMARKS(LOOP,OTH)%spr40)
        BENCHMARKS(LOOP,B)%SSBf90=BENCHMARKS(LOOP,B)%sprF90*EQUILIBRIUM_REC(B,BENCHMARKS(LOOP,B)%sprF90,BENCHMARKS(LOOP,OTH)%sprF90)
        BENCHMARKS(LOOP,B)%YF90  =BENCHMARKS(LOOP,B)%YPRF90*EQUILIBRIUM_REC(B,BENCHMARKS(LOOP,B)%sprF90,BENCHMARKS(LOOP,OTH)%sprF90)
        BENCHMARKS(LOOP,B)%SSB75Fmax=BENCHMARKS(LOOP,B)%spr75Fmax* &
                                     EQUILIBRIUM_REC(B,BENCHMARKS(LOOP,B)%spr75Fmax,BENCHMARKS(LOOP,OTH)%spr75Fmax)
        BENCHMARKS(LOOP,B)%Y75Fmax  =BENCHMARKS(LOOP,B)%YPR75Fmax* &
                                     EQUILIBRIUM_REC(B,BENCHMARKS(LOOP,B)%spr75Fmax,BENCHMARKS(LOOP,OTH)%spr75Fmax)
        BENCHMARKS(LOOP,B)%SPRATIOMSY=BENCHMARKS(LOOP,B)%SPRMSY/virgin(b)
        BENCHMARKS(LOOP,B)%SPRATIOMAX=BENCHMARKS(LOOP,B)%SPRMAX/virgin(b)
        BENCHMARKS(LOOP,B)%SPRATIO01=BENCHMARKS(LOOP,B)%SPR01/virgin(b)
        X(B)=BENCHMARKS(LOOP,B)%FMSY

      END DO
      PR_TYPE=1 ; CALL OBJECTIVE_YR(X,func,NPARM) ; DO B=1,N_BOX ; BENCHMARKS(LOOP,B)%YPRMSY=YRTEMP(B) ; END DO

      IF(LOOP==1) WRITE(*,*) 'END OF SUBROUTINE GETYR'

      RETURN ; END


!-----------------------------------------------------------------------------
      FUNCTION EQUILIBRIUM_REC(I,SPR1,SPR2)
! Computes the expected recruitment at equilibrium corresponding to given S/R values
!-----------------------------------------------------------------------------
      USE STATISTICS ; USE CNTRL
      IMPLICIT NONE
      REAL (KIND=8) :: SPR(2),EQUILIBRIUM_REC,TEMP,SPR1,SPR2
      INTEGER S,I
      SPR(1)=SPR1 ; SPR(2)=SPR2
      IF(SR_BOX>0) THEN
          TEMP=0 ; DO S=1,N_BOX ; IF(SR_BOX==S .OR. SR_BOX>2) TEMP=TEMP+SPR(S)/1000.0 ; END DO
        ELSE
          TEMP=SPR(I)/1000.0
      ENDIF
      SELECT CASE (SR_TYPE(1))
        CASE (1)   ! Beverton & Holt recruitment as function of S/R
                   EQUILIBRIUM_REC=MAX(0.0d0,(sr(2,I)-sr(3,I)/TEMP)/1000.0) ! spr values < sr(3,i)/sr(2,i) are not sustainable
        CASE (2)   ! Ricker
                   IF(sr(2,I)*TEMP>1.0) THEN ; EQUILIBRIUM_REC=(LOG(sr(2,I)*TEMP)/sr(3,I)/TEMP)/1000.0
                     ELSE ; EQUILIBRIUM_REC=0
                   ENDIF
        CASE (3:5) ! two-line curves
                   IF(TEMP>=(SR(3,I)/SR(2,I))) THEN ! Assume maximum recruitment
                       EQUILIBRIUM_REC=SR(2,I)/1000.0   ! divided by 1000. because it is later multiplied by SPR
                     ELSE ! Such small s/r's are not sustainable according to 2line model, s/r ranges from 1/slope to infinity
                       EQUILIBRIUM_REC=0
                   ENDIF
        CASE DEFAULT ; EQUILIBRIUM_REC=SR(2,I)/1000.0   ! constant rec.
      END SELECT
      EQUILIBRIUM_REC=EQUILIBRIUM_REC*BIAS_ADJ(I)
      RETURN
      END

!-----------------------------------------------------------------------------
      SUBROUTINE OBJECTIVE_YR(parm,func,NPARM)
! Note here that the growth curve for Gear 1 is used to represent the average growth
! for all fisheries combined, therefore the fishery represented as 'gear 1' in each area
! should be the one with the growth curve most similar to that of the aggregate
! fishery. Alternatively, one could define gear 1 as a dummy fishery with no historical or
! projected catch. This is not an issue of course when there is only one fishery in each area.
!-----------------------------------------------------------------------------
      USE CNTRL ; USE STATISTICS ; IMPLICIT NONE
      INTEGER :: S,A,NPARM,G,O,MAXAGE,ag,I
      REAL (KIND=8) :: PARM(NPARM),z(2),ftemp(0:100),ntemp(2),d(2),e(2),SURV(2),ndiffusion(2,2), &
                       WEIGHT(2,0:1),func,Ta,tl,PENALTY,temp
      func=0 ; yrtemp=0. ; srtemp=0. ; PENALTY=0 ; d=0; e=0; ndiffusion=0; z=0; WEIGHT=0; SURV=0;
      NTEMP=SEX_FRACTION ! NOTE: SEX_FRACTION==1 if SR_BOX<1
      DO S=1,NPARM
        IF(PARM(S)<0.0) PENALTY=PENALTY+1.0D10*(1+PARM(S)**2)
        IF(PARM(S)>3.0) PENALTY=PENALTY+(PARM(S)-3.0)**2
      END DO
      IF(SR_BOX>0) THEN  ! Sex specific
          MAXAGE=MAX(AGE(1)%LAST,AGE(2)%LAST)
        ELSEIF(N_BOX==1) THEN ! Single stock
          MAXAGE=AGE(1)%LAST
        ELSE ! Two-stock
          MAXAGE=100   ! hardwire to large number because it is hard to deal with a plus-group in mixing models
      ENDIF
      DO A=AGE(1)%FIRST,MAXAGE
       DO S=1,N_BOX
         AG=MIN(A,AGE(s)%LAST) ! fish older than the plus-group have identical parameters
         IF(SR_BOX>0.AND.A>AG) THEN
             EXIT ! Exit allowing plus group of one sex to be younger than the other
           ELSEIF(SR_BOX>0) THEN
             I=1  ! Only one parameter for sex-specific analyses
           ELSEIF(F_MULT>-999) THEN
             I=1  ! F on stock 2 is multiple of F on stock 2
           ELSE
             I=S  ! One parameter for each stock
         ENDIF
         Ftemp(s)=MAX(PARM(I),0.0)*SELECTIVITY(S,0,ag) ;
         IF(S>1 .and. F_MULT>-999) THEN;
           IF(F_MULT<0) THEN; Ftemp(s)=ABS(F_MULT)*Ftemp(s); ELSE; Ftemp(s)=F_MULT*F_CURRENT_RATIO*Ftemp(s); ENDIF
         ENDIF
         Z(S)=M(S,ag,year%last)+Ftemp(S) ; Ftemp(s)=Ftemp(s)*DISCARD_MOD(s,0,ag)
         IF(A==AGE(S)%LAST.AND.(N_BOX==1.OR.SR_BOX>0)) NTEMP(S)=NTEMP(S)/(1-EXP(-Z(S))) ! adj. for + group (not for mixing models)
         ! get weights
         DO G=0,1 ! G=0 refers to ssb in regards to growth, G=1 is the first fishery (with growth equal to overall average)
           IF(WEIGHT_TYPE%FIRST==0 .OR. (WEIGHT_TYPE%FIRST==1 .AND. A>=age(s)%LAST)) THEN
               Ta=DBLE(A)+GROWTH(S,g)%DATE/12.0
               IF(A==AGE(S)%LAST.AND.(N_BOX==1.OR.SR_BOX>0)) Ta=Ta+EXP(-Z(S))/(1-EXP(-Z(S))) ! age of + group (not used for mixing)
               IF(GROWTH(S,G)%CURVE==1) THEN  ! von Bertalanffy/Chapman-Richards growth
                   TL=GROWTH(S,G)%LINF*(1.0-GROWTH(S,G)%M*EXP(-GROWTH(S,G)%K*(Ta-GROWTH(S,G)%T0)))**(1.0/GROWTH(S,G)%M)
                   WEIGHT(S,G)=GROWTH(S,G)%WA*TL**GROWTH(S,G)%WB
                 ELSEIF(GROWTH(S,G)%CURVE==2) THEN  ! Gompertz growth
                   WEIGHT(S,G)=GROWTH(S,g)%LINF*EXP(-1.0*EXP(-GROWTH(S,g)%K*(Ta-GROWTH(S,g)%T0)))
               ENDIF
             ELSE
               WEIGHT(S,G)=w(S,G,A,YEAR%LAST)
           ENDIF
         END DO
         d(S)=ftemp(s)*(1.0-EXP(-z(S)))/z(S)
         e(S)=maturity(S,AG)*Weight(s,0)
         surv(s)=EXP(-z(S)*GROWTH(S,0)%DATE/12.0)
       END DO
       ! get per-recruit statistics
       DO S=1,N_BOX
         AG=MIN(A,AGE(s)%LAST) ; IF(SR_BOX>0.AND.A>AG) EXIT
         IF(s==1) THEN ; O=2 ; ELSE ; O=1 ; ENDIF
         IF(MODEL_TYPE==1) THEN
             IF(N_BOX==1.OR.SR_BOX>0) THEN  ! one-box or sex-specific model
                 TEMP=NTEMP(s)*(1.0-T(s,ag))+NTEMP(o)*T(o,ag)
                 YRTEMP(S)=YRTEMP(S)+TEMP*Weight(S,1)*d(S)
                 SRTEMP(S)=SRTEMP(S)+TEMP*e(S)*surv(s)     ! assumes all fish spawn where they are
               ELSE ! diffusion two-box model
                 IF(A==AGE(1)%FIRST) THEN ;  ndiffusion(s,s)=1.0-T(s,ag);  ndiffusion(s,o)=T(s,ag)  ; ENDIF
                 YRTEMP(S)=YRTEMP(S)+ndiffusion(s,s)*d(s)*Weight(s,1)+ndiffusion(s,o)*d(o)*Weight(o,1)
                 SRTEMP(S)=SRTEMP(S)+ndiffusion(s,s)*e(s)*surv(s)+ndiffusion(s,o)*e(o)*surv(o)     ! assumes all fish spawn where they are
             ENDIF
           ELSE ! overlap-- ntemp(s) represents numbers belonging to stock s (regardless of area they are in)
             TEMP=1.0-T(s,ag)
             YRTEMP(S)=YRTEMP(S)+Weight(s,1)*NTEMP(s)*(TEMP*d(s)   +T(s,ag)*d(o))
             SRTEMP(S)=SRTEMP(S)+e(s)       *NTEMP(s)*(TEMP*surv(s)+T(s,ag)*surv(o)) ! assumes all fish go back home to spawn
         ENDIF
       END DO
       DO S=1,N_BOX
         AG=MIN(A,AGE(s)%LAST) ; IF(SR_BOX>0.AND.A>AG) EXIT
         IF(s==1) THEN ; O=2 ; ELSE ; O=1 ; ENDIF
         IF(MODEL_TYPE==1) THEN
             IF(N_BOX==1.OR.SR_BOX>0) THEN ! one-box or sex-specific model
                 NTEMP(S)=(NTEMP(s)*(1.0-T(s,ag))+NTEMP(o)*T(o,ag))*EXP(-z(S))
               ELSE ! diffusion two-box model
                 ndiffusion(s,s)=ndiffusion(s,s)*EXP(-z(s))*(1.0-T(s,ag)) + ndiffusion(s,o)*EXP(-z(o))*T(o,ag)
                 ndiffusion(s,o)=ndiffusion(s,s)*EXP(-z(s))*T(s,ag)       + ndiffusion(s,o)*EXP(-z(o))*(1-T(o,ag))
             ENDIF
           ELSE ! overlap-- ntemp(s) represents numbers belonging to stock s (regardless of area they are in)
             NTEMP(S)=NTEMP(s)*((1.0-T(s,ag))*EXP(-z(S))+T(s,ag)*EXP(-z(o)))
         ENDIF
       END DO
      END DO ! a
      IF(SR_BOX>=1) THEN ; TEMP=0 ; DO S=1,N_BOX ; IF(SR_BOX==S .OR. SR_BOX>2) TEMP=TEMP+SRTEMP(S) ; END DO ; ENDIF ! sex-specific
      DO S=1,N_BOX
!       yr_box tells whether to maximize y/r or msy with respect to one or both stocks
        IF(SR_BOX<1) TEMP=SRTEMP(S) ! not sex-specific
        IF(TEMP<0 .OR. YRTEMP(s)<0) THEN
            func=1E+32
          ELSEIF(PR_TYPE==1) THEN
!           {Compute y/r for combined stocks or focusing on one stock in particular}
            IF(YR_BOX==S .OR. YR_BOX<1 .OR. YR_BOX>2) FUNC=FUNC-YRTEMP(S)
          ELSEIF(PR_TYPE==2) THEN
!           {Compute msy for combined stocks or focusing on one stock in particular}
            IF(SR_TYPE(s)==1) THEN ! Beverton and Holt
                SSBMSYtemp(s)=bias_adj(s)*sr(2,s)*TEMP/1000.0-sr(3,s)
              ELSEIF(SR_TYPE(s)==2) THEN ! Ricker
                SSBMSYtemp(s)=LOG(bias_adj(s)*sr(2,s)*TEMP/1000.0)/sr(3,s)
            ENDIF
            IF(TEMP>0) THEN ; MSYtemp(s)=SSBMSYtemp(s)*YRTEMP(s)/TEMP ; ELSE ; MSYTEMP(s)=0 ; ENDIF
            IF(YR_BOX==S .OR. YR_BOX<1 .OR. YR_BOX>2) FUNC=FUNC-MSYtemp(s)
        END IF
      END DO ! s
      FUNC=FUNC+PENALTY
      RETURN ; END

!-----------------------------------------------------------------------------
      SUBROUTINE OBJECTIVE_SR(parm,func,NPARM)
      ! Compute likelihood components when fitting Beverton-Holt or ricker functions
!-----------------------------------------------------------------------------
      USE CNTRL ; USE STATISTICS ; IMPLICIT NONE
      INTEGER :: I,COUNTS,Y,NPARM,YS
      REAL (KIND=8) :: PARM(NPARM),RSS,FUNC,PRED,maxR,PENALTY,NEGLOGLIKELIHOOD,MSE,TEMP,REC,RESIDUAL,RESIDUAL_LAG1,MULT,RHO
      RSS=0 ; FUNC=0 ; COUNTS=0 ; MAXR=0; PENALTY=0; RESIDUAL_LAG1=0;
      IF(ESTIMATE_AUTOCORRELATION>0) THEN ; RHO=PARM(NPARM); ELSE ; RHO=SR(5,BOX%STORE); ENDIF
      MULT=1.0-RHO**2.;

      DO I=1,2 ; IF(PARM(I)<=0) THEN ; FUNC=1.D32 ; RETURN ; ENDIF ; END DO
      IF(ABS(RHO)>1.0) THEN ; FUNC=1.D32 ; RETURN ; ENDIF

      DO Y=1,YEAR%LAST

        IF(SR_BOX>0) THEN ; REC=N(1,age(1)%FIRST,Y)+N(2,age(2)%FIRST,Y) ; ELSE ; REC=N(BOX%STORE,age(BOX%STORE)%FIRST,Y) ; ENDIF
        IF(SR_PENALTY(BOX%STORE)>0 .AND. REC_PENALTY(1)>0 .AND. Y>=REC_PENALTY(1) .AND. Y<=REC_PENALTY(2)) &
          maxR=maxR+REC/(REC_PENALTY(2)-REC_PENALTY(1)+1) ! arithmetic average recruitment
        IF(ERRORTYPE==1) REC=LOG(REC)

        YS=Y-age(1)%FIRST
        IF (YS>=YEAR%FIRST .AND. Y>=REC_YEAR(1) .AND. Y<=REC_YEAR(2)) THEN
          SELECT CASE (SR_BOX)
              CASE (1) ;     TEMP=SSB(1,YS)
              CASE (2) ;     TEMP=SSB(2,YS)
              CASE (3) ;     TEMP=SSB(1,YS)+SSB(2,YS)
              CASE DEFAULT ; TEMP=SSB(BOX%STORE,YS)
          END SELECT
          IF(SR_TYPE(1)==1) THEN ! {Beverton and Holt curve}
              PRED=PARM(1)*TEMP/(TEMP+PARM(2))
              IF(Y==REC_YEAR(2)) PENALTY=PARM(1) ! predicted maximum recruitment
            ELSEIF(SR_TYPE(1)==2) THEN ! {Ricker curve}
              PRED=PARM(1)*TEMP*EXP(-PARM(2)*TEMP)
              IF(Y==REC_YEAR(2)) PENALTY=PARM(1)*EXP(-1.0)/PARM(2) ! predicted maximum recruitment
          END IF
          IF(ERRORTYPE==1) THEN ; PRED=LOG(PRED) ; IF(PENALTY>0) PENALTY=LOG(PENALTY) ; ENDIF
          RESIDUAL=PRED-REC; sr(5+Y,BOX%STORE)=RESIDUAL
          IF(SR_PENALTY(BOX%STORE)>0 .AND. REC_PENALTY(1)<=0) maxR=MAX(REC,maxR)
          RSS=RSS+MULT*(RESIDUAL-RHO*RESIDUAL_LAG1)**2. ! first instance RSS=(1-rho^2)*(residual)^2
          COUNTS=COUNTS+1
          RESIDUAL_LAG1=RESIDUAL ! Previous year's residual
          MULT=1.0
        ENDIF

      END DO
      MSE=MAX(RSS/DBLE(COUNTS),1.0D-32)
      SR(4,BOX%STORE)=MAX(SQRT(MSE),SR_mse_minimum(BOX%STORE))
      NEGLOGLIKELIHOOD=0.5*DBLE(COUNTS)*(1.0+LOG(MSE))-0.5*LOG(1-RHO**2.)
      FUNC=NEGLOGLIKELIHOOD
      IF(SR_PENALTY(BOX%STORE)>0) THEN ! add penalty
        IF(REC_PENALTY(1)>0 .AND. ERRORTYPE==1) maxR=LOG(maxR)
        IF(PENALTY>maxR) FUNC=FUNC+0.5*((PENALTY-maxR)**2)/VAR_R(BOX%STORE)
      ENDIF
      RETURN ; END

!-----------------------------------------------------------------------------
      SUBROUTINE OBJECTIVE_AUTO(parm,func,NPARM)
      ! Compute autocorrelation for two-line models
!-----------------------------------------------------------------------------
      USE CNTRL ; USE STATISTICS ; IMPLICIT NONE
      INTEGER :: I,COUNTS,Y,NPARM,YS
      REAL (KIND=8) :: PARM(NPARM),RSS,FUNC,PRED,maxR,PENALTY,NEGLOGLIKELIHOOD,MSE,TEMP,REC,RESIDUAL,RESIDUAL_LAG1,MULT,RHO
      RSS=0 ; FUNC=0 ; COUNTS=0 ; MAXR=0; PENALTY=0; RESIDUAL_LAG1=0;
      RHO=PARM(NPARM) ; IF(ABS(RHO)>1.0) THEN ; FUNC=1.D32 ; RETURN ; ENDIF
      MULT=1.0-RHO**2.;

      DO Y=1,YEAR%LAST

        IF(SR_BOX>0) THEN ; REC=N(1,age(1)%FIRST,Y)+N(2,age(2)%FIRST,Y) ; ELSE ; REC=N(BOX%STORE,age(BOX%STORE)%FIRST,Y) ; ENDIF
        IF(ERRORTYPE==1) REC=LOG(REC)

        YS=Y-age(1)%FIRST
        IF (YS>=YEAR%FIRST .AND. Y>=REC_YEAR(1) .AND. Y<=REC_YEAR(2)) THEN
          SELECT CASE (SR_BOX)
              CASE (1) ;     TEMP=SSB(1,YS)
              CASE (2) ;     TEMP=SSB(2,YS)
              CASE (3) ;     TEMP=SSB(1,YS)+SSB(2,YS)
              CASE DEFAULT ; TEMP=SSB(BOX%STORE,YS)
          END SELECT
          IF(TEMP<SR(3,BOX%STORE)) THEN ; PRED=SR(2,BOX%STORE)*TEMP/SR(3,BOX%STORE) ; ELSE ; PRED=SR(2,BOX%STORE) ; ENDIF
          IF(ERRORTYPE==1) THEN ; PRED=LOG(PRED) ; IF(PENALTY>0) PENALTY=LOG(PENALTY) ; ENDIF
          RESIDUAL=REC-PRED; sr(5+Y,BOX%STORE)=RESIDUAL
          RSS=RSS+MULT*(RESIDUAL-RHO*RESIDUAL_LAG1)**2. ! first instance RSS=(1-rho^2)*(residual)^2
          COUNTS=COUNTS+1
          RESIDUAL_LAG1=RESIDUAL ! Previous year's residual
          MULT=1.0
        ENDIF

      END DO
      MSE=MAX(RSS/DBLE(COUNTS),1.0D-32)
      SR(4,BOX%STORE)=MSE
      NEGLOGLIKELIHOOD=0.5*DBLE(COUNTS)*(1.0+LOG(MSE))-LOG(1-RHO**2.)
      FUNC=NEGLOGLIKELIHOOD
      RETURN ; END

!-----------------------------------------------------------------------------
      SUBROUTINE OBJECTIVE_F(parm,func,NPARM)
!-----------------------------------------------------------------------------
      USE CNTRL ; USE STATISTICS ; IMPLICIT NONE
      INTEGER :: I,B,A,Y,O,G,NPARM,YS,its,j,jj,bb,gg
      REAL (KIND=8) :: PARM(NPARM),FUNC,PENALTY,TEMP,z(2,0:50),N0(2,0:50),e(2,0:50),fvec(nparm),S(100,0:50),FF(200)

      FUNC=0 ; DO I=1,NPARM ; IF(PARM(I)<=0) THEN ; FUNC=1.D32 ; RETURN ; ENDIF ; END DO
      B=BOX%STORE ; Y=YEAR%STORE ; z=0
      ! get abundance in each area after mixing
      IF(B==1) THEN ; O=2 ; ELSE ; O=1 ; ENDIF
      DO A=age(b)%FIRST,age(b)%Last
        N0(b,A)=N(b,A,Y)*(1.0-T(b,a))+N(o,a,y)*T(o,a)
      END DO
      ! assign parameter values
      J=0
      DO G=1,N_GEARS(B)
        IF(SIGNAL(B,G)>0) THEN ; j=j+1 ; ff(g)=parm(j) ; ELSE ; ff(g)=F_APICAL(b,g,y) ; ENDIF
        DO A=age(b)%FIRST,age(b)%Last ; z(b,a)=z(b,a)+ff(g)*SELECTIVITY(b,g,A) ; END DO
        IF(SR_BOX>0) THEN ; DO A=age(o)%FIRST,age(o)%Last ; z(o,a)=z(o,a)+ff(g)*SELECTIVITY(o,g,A) ; END DO ; ENDIF
      END DO
      DO A=age(b)%FIRST,age(b)%Last ; z(b,a)=z(b,a)+m(b,a,y) ; e(b,a)=EXP(-z(b,A)) ; END DO
      IF(SR_BOX>0) THEN ; DO A=age(o)%FIRST,age(o)%Last ; z(o,a)=z(o,a)+m(o,a,y) ; e(o,a)=EXP(-z(o,A)) ; END DO ; ENDIF

      FVEC=0 ; J=0
      DO G=1,N_GEARS(B)
        ! function evaluation (fvec(j) = catch expected with given value of parm(j))
        IF(SIGNAL(B,G)>0) THEN
          j=j+1
          DO A=age(b)%FIRST, age(b)%LAST
            fvec(j)=fvec(j)+N0(b,a)*ff(g)*SELECTIVITY(b,g,A)*DISCARD_MOD(b,g,a)*(1.-e(b,a))/z(b,a)*W(b,g,A,Y)
          END DO
          IF(SR_BOX>0) THEN
            DO A=age(O)%FIRST, age(O)%LAST
              fvec(j)=fvec(j)+N(o,a,y)*ff(g)*SELECTIVITY(o,g,A)*DISCARD_MOD(o,g,a)*(1.-e(o,a))/z(o,a)*W(o,g,A,Y)
            END DO
          ENDIF
          func=func+(fvec(j)-TAC(b,g,y))**2
        END IF
      END DO
      RETURN ; END

!-----------------------------------------------------------------------------
      SUBROUTINE GETBIOMASS(y)
!-----------------------------------------------------------------------------
! determines the relative exploitable biomass (relative selectivity
! multiplied by biomass at age) for all years in metric tons
! This routine uses growth coefficients from fishery 1
      USE CNTRL ; USE STATISTICS ; IMPLICIT NONE
      INTEGER :: S,A,Y,O
      REAL (KIND=8) :: Z,E
      DO s=1,N_BOX
        IF(S==1) THEN ; O=2 ; ELSE ; O=1 ; ENDIF
        BIOMASS(S,Y)=0 ; FBIOMASS(S,Y)=0
        DO A=age(s)%FIRST,age(s)%LAST
          Z=M(S,A,y)+F(s,a,y) ; E=(1-EXP(-Z))/Z/1000.0
          BIOMASS(S,Y) = BIOMASS(S,Y)+(N(S,A,Y)*(1.0-T(s,a))+N(o,A,Y)*T(o,a))*W(S,1,A,Y)*E
          IF(F100(S,Y)>0) FBIOMASS(S,Y)=FBIOMASS(S,Y)+(N(S,A,Y)*(1.0-T(s,a))+N(o,A,Y)*T(o,a))*W(S,1,A,Y)*E*F(s,a,y)/F100(s,y)
        END DO
      END DO

      RETURN ; END

!-----------------------------------------------------------------------------
      SUBROUTINE GETSSB(Y)
!-----------------------------------------------------------------------------
      USE CNTRL ; USE STATISTICS ; IMPLICIT NONE
      INTEGER :: S,A,Y,I,o
      REAL (KIND=8) :: E,EO
      DO s=1,N_BOX
        IF(S==1) THEN ; O=2 ; ELSE ; O=1 ; ENDIF
        SSB(S,Y)=0 ; SSN(S,Y)=0 ; I=0 ; IF(age(s)%FIRST==0) I=1
        DO A=age(s)%FIRST+I,age(s)%LAST
          E=EXP((-M(S,A,y)-F(s,a,y))*GROWTH(S,0)%DATE/12.0)
          IF(MODEL_TYPE==1) THEN
              SSBa(S,A,Y)=(N(S,A,Y)*(1.0-T(s,a))+N(o,A,Y)*T(o,a))*MATURITY(S,A)*E*W(S,0,A,Y)/1000.0
              SSB(S,Y)=SSB(S,Y)+SSBa(S,A,Y)
              SSN(S,Y)=SSN(S,Y)+(N(S,A,Y)*(1.0-T(s,a))+N(o,A,Y)*T(o,a))*MATURITY(S,A)*E
              !WRITE(*,*)  N(S,A,Y),W(S,0,A,Y),M(S,A,y),F(s,a,y)
            ELSE
              EO=EXP((-M(O,A,y)-F(O,a,y))*GROWTH(O,0)%DATE/12.0)
              SSBa(S,A,Y)=N(S,A,Y)*((1.0-T(s,a))*E+T(S,a)*EO)*MATURITY(S,A)*W(S,0,A,Y)/1000.0
              SSB(S,Y)=SSB(S,Y)+SSBa(S,A,Y)
              SSN(S,Y)=SSN(S,Y)+N(S,A,Y)*((1.0-T(s,a))*E+T(S,a)*EO)*MATURITY(S,A)
          ENDIF
        END DO
        IF(SSB(S,Y)<0.000001) SSB(S,Y)=0.000001
      END DO
      RETURN ; END

!-----------------------------------------------------------------------------
      SUBROUTINE GETFAA(Y)
!-----------------------------------------------------------------------------
      USE CNTRL ; USE STATISTICS ; IMPLICIT NONE
      INTEGER :: B,A,Y,G

      DO B=1,N_BOX
        F100(b,Y)=0
        DO A=age(b)%FIRST,age(b)%LAST
          F(b,a,y)=0
          DO G=1,N_GEARS(B)
            Fp(b,g,A,Y)=F_APICAL(b,g,y)*SELECTIVITY(b,g,A)
            F(b,a,y)=F(b,a,y)+Fp(b,g,a,y)
            Fp(b,g,A,Y)=Fp(b,g,A,Y)*DISCARD_MOD(b,g,a)
          END DO
          F100(b,Y)=MAX(F100(b,Y),F(b,a,y))
        END DO
      END DO

      RETURN ; END

!-----------------------------------------------------------------------------
      SUBROUTINE RESULTS
! summarize output
!-----------------------------------------------------------------------------
      USE CNTRL ; USE STATISTICS, ONLY : BENCHMARKS
      IMPLICIT NONE
      INTEGER :: NCOL,S,IC,I,NS,J
      CHARACTER (LEN=12) :: FILNAME
      CHARACTER (LEN=6) :: IDBOX(2)
      CHARACTER (LEN=4) :: NOBOX
      DATA IDBOX/'-1.STA','-2.STA'/
      DATA NOBOX/'.STA'/

      write(*,*) 'COMPUTING SUMMARY STATISITICS'

      DO S=1,N_BOX

        ic=20*S
        DO J=1,8 ; REWIND(ic+J) ; END DO

        I=0
        WRITE(*,*) 'SUMMARIZING BENCHMARKS'
        NCOL=43; NS=1
        IF(N_BOX==1) THEN; FILNAME='BENCH'//NOBOX; ELSE; FILNAME='BENCH'//IDBOX(s); ENDIF
        CALL SUMMARIZE(IC,I,NCOL,ci,FILNAME,N_boots,NS,YEAR%display)

        NCOL=YEAR%PRJCT ; NS=N_SCENARIOS
        DO J=1,8
          IF(j==8) THEN ; NCOL=n_gears(s) ; NS=1 ; ENDIF
          WRITE(*,'(1x,a12,a8)') 'SUMMARIZING ',CLASS(J)
          IF(N_BOX==1) THEN; FILNAME=CLASS(J)//NOBOX; ELSE; FILNAME=CLASS(J)//IDBOX(S); ENDIF
          CALL SUMMARIZE(J+IC,I,NCOL,ci,FILNAME,N_boots,NS,YEAR%display)
        END DO

        I=20*s ; NCOL=YEAR%PRJCT
        WRITE(*,*) 'SUMMARIZING SPAWNING STOCK RELATIVE TO REFERENCE POINTS'
        REWIND(1+ic)
        IF(N_BOX==1) THEN; FILNAME='SSmsy'//NOBOX; ELSE; FILNAME='SSmsy'//IDBOX(S); ENDIF
        CALL SUMMARIZE(1+IC,I,NCOL,ci,FILNAME,N_boots,N_SCENARIOS,YEAR%display)
        REWIND(1+ic)
        i=81+10*s;
        IF(N_BOX==1) THEN; FILNAME='SSf01'//NOBOX; ELSE; FILNAME='SSf01'//IDBOX(S); ENDIF
        CALL SUMMARIZE(1+IC,I,NCOL,ci,FILNAME,N_boots,N_SCENARIOS,YEAR%display)
        REWIND(1+ic)
        i=82+10*s;
        IF(N_BOX==1) THEN; FILNAME='SSspr1'//NOBOX; ELSE; FILNAME='SSspr1'//IDBOX(S); ENDIF
        CALL SUMMARIZE(1+IC,I,NCOL,ci,FILNAME,N_boots,N_SCENARIOS,YEAR%display)
        REWIND(1+ic)
        i=83+10*s;
        IF(N_BOX==1) THEN; FILNAME='SSspr2'//NOBOX; ELSE; FILNAME='SSspr2'//IDBOX(S); ENDIF
        CALL SUMMARIZE(1+IC,I,NCOL,ci,FILNAME,N_boots,N_SCENARIOS,YEAR%display)
        REWIND(1+ic)
        i=84+10*s;
        IF(N_BOX==1) THEN; FILNAME='SSspr3'//NOBOX; ELSE; FILNAME='SSspr3'//IDBOX(S); ENDIF
        CALL SUMMARIZE(1+IC,I,NCOL,ci,FILNAME,N_boots,N_SCENARIOS,YEAR%display)
        REWIND(1+ic)
        i=85+10*s;
        IF(N_BOX==1) THEN; FILNAME='SSfmax'//NOBOX; ELSE; FILNAME='SSfmax'//IDBOX(S); ENDIF
        CALL SUMMARIZE(1+IC,I,NCOL,ci,FILNAME,N_boots,N_SCENARIOS,YEAR%display)

        DO J=1,8 ; CLOSE(ic+J) ; END DO
        CLOSE(ic+31)

      END DO ! STOCK
      CLOSE(98); CLOSE(99); CLOSE(70)

      RETURN ; END

!----------------------------------------------------------------------------
      SUBROUTINE SUMMARIZE(IFILE,IFILE2,NCOL,CI,FILNAME,Nboots,nscenarios,display)
!----------------------------------------------------------------------------
!  Program to summarize variables : col=var, row=observation
      USE STATISTICS , ONLY : BENCHMARKS
      USE CNTRL , ONLY : SR_BOX,SPR_TARGET
      IMPLICIT NONE
      REAL (KIND=8) :: x(500,2000),t(2000),median(500),hi(500),low(500),ci,alpha,divisor,mean(1000),std(1000)
      INTEGER IFILE,IFILE2,Nboots,LOOP,ncol,SCENARIO,s,display,NSCENARIOS,ILOW,IHI,IMED,IVAR,I,J,b
      CHARACTER (LEN=12) FILNAME,dummy
      OPEN(7,FILE=FILNAME,STATUS='UNKNOWN')
      if(Nboots>0) then
        alpha=(1.-ci/100.)/2.
        if(real(Nboots)*alpha.lt.1) THEN
          alpha=1/real(Nboots)
          if(Nboots.gt.0) then
            write(*,*) 'Warning: there were not enough data points to construct the requested confidence intervals.'
            write(*,*) '  The confidence intervals have been changed to ',(1-2.*alpha)*100,' percent'
          endif
        endif
        ilow=int(dble(Nboots)*alpha+.51)
        imed=int(dble(Nboots)/2+.51)
        ihi=int(dble(Nboots)*(1-alpha)+1.51)
      endif
      DO S=1,nscenarios
       DO  I=1,Nboots+1
        IF(ifile==20 .OR. ifile==40) THEN
          b=ifile/20
          X(1,I)=BENCHMARKS(I,B)%FMSY ; X(2,I)=BENCHMARKS(I,B)%MSY ; X(3,I)=BENCHMARKS(I,B)%YPRMSY ; X(4,I)=BENCHMARKS(I,B)%SPRMSY
          X(5,I)=BENCHMARKS(I,B)%SPRATIOMSY ; X(6,I)=BENCHMARKS(I,B)%SSBMSY ; X(7,I)=BENCHMARKS(I,B)%FMAX
          X(8,I)=BENCHMARKS(I,B)%YPRMAX ; X(9,I)=BENCHMARKS(I,B)%SPRMAX ; X(10,I)=BENCHMARKS(I,B)%SPRATIOMAX
          X(11,I)=BENCHMARKS(I,B)%SSBMAX ; X(12,I)=BENCHMARKS(I,B)%F01 ; X(13,I)=BENCHMARKS(I,B)%YPR01
          X(14,I)=BENCHMARKS(I,B)%SPR01 ; X(15,I)=BENCHMARKS(I,B)%SPRATIO01 ; X(16,I)=BENCHMARKS(I,B)%SSBF01
          X(17,I)=BENCHMARKS(I,B)%F20; X(18,I)=BENCHMARKS(I,B)%YPR20 ; X(19,I)=BENCHMARKS(I,B)%SPR20 ; X(20,I)=BENCHMARKS(I,B)%SSB20
          X(21,I)=BENCHMARKS(I,B)%F30; X(22,I)=BENCHMARKS(I,B)%YPR30 ; X(23,I)=BENCHMARKS(I,B)%SPR30 ; X(24,I)=BENCHMARKS(I,B)%SSB30
          X(25,I)=BENCHMARKS(I,B)%F40; X(26,I)=BENCHMARKS(I,B)%YPR40 ; X(27,I)=BENCHMARKS(I,B)%SPR40 ; X(28,I)=BENCHMARKS(I,B)%SSB40
          X(29,I)=BENCHMARKS(I,B)%F90; X(30,I)=BENCHMARKS(I,B)%Yf90 ; X(31,I)=BENCHMARKS(I,B)%yPRf90; X(32,I)=BENCHMARKS(I,B)%Sprf90
          X(33,I)=BENCHMARKS(I,B)%ssbf90; X(34,I)=BENCHMARKS(I,B)%F75fmax; X(35,I)=BENCHMARKS(I,B)%Y75fmax ;
          X(36,I)=BENCHMARKS(I,B)%yPR75fmax ; X(37,I)=BENCHMARKS(I,B)%Spr75fmax; X(38,I)=BENCHMARKS(I,B)%ssb75fmax
          X(39,I)=BENCHMARKS(I,B)%ymax ; X(40,I)=BENCHMARKS(I,B)%Y01; X(41,I)=BENCHMARKS(I,B)%Y20
          X(42,I)=BENCHMARKS(I,B)%y30 ; X(43,I)=BENCHMARKS(I,B)%Y40;
         ELSE
          READ(IFILE,*) scenario,loop,(X(j,I),j=1,ncol) ; loop=loop+1
        ENDIF
        IF(IFILE2.EQ.20 .OR. IFILE2.EQ.40) THEN
!           {divisor is ssb at MSY}
            b=ifile2/20      ; DIVISOR=BENCHMARKS(loop,b)%SSBMSY
          ELSEIF(IFILE2.eq.91 .or. IFILE2.eq.101) THEN
            b=(ifile2-81)/10 ; DIVISOR=BENCHMARKS(loop,b)%SSBf01
          ELSEIF(IFILE2.eq.92  .or. IFILE2.eq.102) THEN
            b=(ifile2-82)/10 ; DIVISOR=BENCHMARKS(loop,b)%SSB20
          ELSEIF(IFILE2.eq.93  .or. IFILE2.eq.103) THEN
            b=(ifile2-83)/10 ; DIVISOR=BENCHMARKS(loop,b)%SSB30
          ELSEIF(IFILE2.eq.94 .or. IFILE2.eq.104) THEN
            b=(ifile2-84)/10 ; DIVISOR=BENCHMARKS(loop,b)%SSB40
          ELSEIF(IFILE2.eq.95 .or. IFILE2.eq.105) THEN
            b=(ifile2-85)/10 ; DIVISOR=BENCHMARKS(loop,b)%SSBMAX
          ELSEIF(IFILE2.LT.0) THEN
!           {divisor is the value during a particular year}
            DIVISOR=X(IABS(IFILE2),I)
          ELSE
            DIVISOR=1.
        ENDIF
        DO J=1,NCOL
          IF(DIVISOR.GT.0) THEN
             X(J,I)=X(J,I)/DIVISOR
           ELSE
             X(J,I)=-1
          ENDIF
        END DO
       END DO
       DO ivar=1,ncol
         std(ivar)=0.d0 ; mean(ivar)=0.d0
         IF(Nboots>0) THEN
              do 11 i=2,Nboots+1 ! exclude loop 1, which contains deterministic result
                  t(i-1)=x(ivar,i)
                  mean(ivar)=mean(ivar)+x(ivar,i)/dble(Nboots)
11            continue
              do 12 i=2,Nboots+1 ! exclude loop 1, which contains deterministic result
                if(Nboots.gt.1) std(ivar)=std(ivar)+(x(ivar,i)-mean(ivar))**2.d0
12            continue
              if(Nboots.gt.1) then
                std(ivar)=dsqrt(std(ivar)/dble(Nboots-1))
                call sort(Nboots,t)
              endif
              low(ivar)=t(ilow)
              hi(ivar)=t(ihi)
              median(ivar)=t(imed)
              if(mod(Nboots,2).eq.0) median(ivar)=(t(imed)+t(imed+1))/2
           ELSE
              mean(ivar)=x(ivar,1) ; median(ivar)=x(ivar,1) ; low(ivar)=x(ivar,1) ; hi(ivar)=x(ivar,1)
         ENDIF
       END DO ! ivar
       IF(IFILE.EQ.20.or.IFILE.EQ.40) THEN
        WRITE(7,'(13X,A15,1X,99(A12,2X))') 'MEASURE        ','LOWER_CL ','MEDIAN  ','UPPER_CL ','AVERAGE ','  RUN_0  ','STD_DEV. '
        write(7,101) 'F at MSY       ',low(1),median(1),hi(1),mean(1),X(1,1),std(1)
        write(7,101) 'MSY            ',low(2),median(2),hi(2),mean(2),X(2,1),std(2)
        write(7,101) 'Y/R at MSY     ',low(3),median(3),hi(3),mean(3),X(3,1),std(3)
        write(7,101) 'S/R at MSY     ',low(4),median(4),hi(4),mean(4),X(4,1),std(4)
        write(7,101) 'SPR AT MSY     ',low(5),median(5),hi(5),mean(5),X(5,1),std(5)
        write(7,101) 'SS AT MSY      ',low(6),median(6),hi(6),mean(6),X(6,1),std(6)
        write(7,101) 'F at max. Y/R  ',low(7),median(7),hi(7),mean(7),X(7,1),std(7)
        write(7,101) 'Y at max. Y/R  ',low(39),median(39),hi(39),mean(39),X(39,1),std(39)
        write(7,101) 'Y/R maximum    ',low(8),median(8),hi(8),mean(8),X(8,1),std(8)
        write(7,101) 'S/R at Fmax    ',low(9),median(9),hi(9),mean(9),X(9,1),std(9)
        write(7,101) 'SPR at Fmax    ',low(10),median(10),hi(10),mean(10),X(10,1),std(10)
        write(7,101) 'SS at Fmax     ',low(11),median(11),hi(11),mean(11),X(11,1),std(11)
        write(7,101) 'F 0.1          ',low(12),median(12),hi(12),mean(12),X(12,1),std(12)
        write(7,101) 'Y at F0.1      ',low(40),median(40),hi(40),mean(40),X(40,1),std(40)
        write(7,101) 'Y/R at F0.1    ',low(13),median(13),hi(13),mean(13),X(13,1),std(13)
        write(7,101) 'S/R at F0.1    ',low(14),median(14),hi(14),mean(14),X(14,1),std(14)
        write(7,101) 'SPR at F0.1    ',low(15),median(15),hi(15),mean(15),X(15,1),std(15)
        write(7,101) 'SS at F0.1     ',low(16),median(16),hi(16),mean(16),X(16,1),std(16)
        J=INT(100*SPR_TARGET(1))
        write(7,102) 'F',J,'% SPR      ',low(17),median(17),hi(17),mean(17),X(17,1),std(17)
        write(7,104) 'Y at  F',J,'%SPR ',low(41),median(41),hi(41),mean(41),X(41,1),std(41)
        write(7,103) 'Y/R at F',J,'%SPR',low(18),median(18),hi(18),mean(18),X(18,1),std(18)
        write(7,103) 'S/R at F',J,'%SPR',low(19),median(19),hi(19),mean(19),X(19,1),std(19)
        write(7,104) 'SS at F',J,'%SPR ',low(20),median(20),hi(20),mean(20),X(20,1),std(20)
        J=INT(100*SPR_TARGET(2))
        write(7,102) 'F',J,'% SPR      ',low(21),median(21),hi(21),mean(21),X(21,1),std(21)
        write(7,104) 'Y at  F',J,'%SPR ',low(42),median(42),hi(42),mean(42),X(42,1),std(42)
        write(7,103) 'Y/R at F',J,'%SPR',low(22),median(22),hi(22),mean(22),X(22,1),std(22)
        write(7,103) 'S/R at F',J,'%SPR',low(23),median(23),hi(23),mean(23),X(23,1),std(23)
        write(7,104) 'SS at F',J,'%SPR ',low(24),median(24),hi(24),mean(24),X(24,1),std(24)
        J=INT(100*SPR_TARGET(3))
        write(7,102) 'F',J,'% SPR      ',low(25),median(25),hi(25),mean(25),X(25,1),std(25)
        write(7,104) 'Y at  F',J,'%SPR ',low(43),median(43),hi(43),mean(43),X(43,1),std(43)
        write(7,103) 'Y/R at F',J,'%SPR',low(26),median(26),hi(26),mean(26),X(26,1),std(26)
        write(7,103) 'S/R at F',J,'%SPR',low(27),median(27),hi(27),mean(27),X(27,1),std(27)
        write(7,104) 'SS at F',J,'%SPR ',low(28),median(28),hi(28),mean(28),X(28,1),std(28)
        write(7,101) 'F 90% max Y/R  ',low(29),median(29),hi(29),mean(29),X(29,1),std(29)
        write(7,101) 'Y 90% max Y/R  ',low(30),median(30),hi(30),mean(30),X(30,1),std(30)
        write(7,101) 'Y/R 90% max Y/R',low(31),median(31),hi(31),mean(31),X(31,1),std(31)
        write(7,101) 'S/R 90% max Y/R',low(32),median(32),hi(32),mean(32),X(32,1),std(32)
        write(7,101) 'SS 90% max Y/R ',low(33),median(33),hi(33),mean(33),X(33,1),std(33)
        write(7,101) 'F 75% of Fmax  ',low(34),median(34),hi(34),mean(34),X(34,1),std(34)
        write(7,101) 'Y 75% of Fmax  ',low(35),median(35),hi(35),mean(35),X(35,1),std(35)
        write(7,101) 'Y/R at 75% Fmax',low(36),median(36),hi(36),mean(36),X(36,1),std(36)
        write(7,101) 'S/R at 75% Fmax',low(37),median(37),hi(37),mean(37),X(37,1),std(37)
        write(7,101) 'SS at 75% Fmax ',low(38),median(38),hi(38),mean(38),X(38,1),std(38)
        IF(SR_BOX>0) WRITE(7,*) '     Note: All statistics are sex-specific except F and SSB values.'
       ELSE
         IF(IFILE==28 .OR. IFILE==48) THEN
            IF(SCENARIO.EQ.1) WRITE(7,'(1X,99(A10,4X))') 'SCENARIO','FLEET   ','LOWER_CL ','MEDIAN  ','UPPER_CL ','AVERAGE ',  &
                                                     ' RUN_0  ','STD_DEV. '
            DO j=1,ncol ; write(7,100) scenario,j,low(j),median(j),hi(j),mean(j),X(j,1),std(j) ; END DO
           ELSE
            IF(SCENARIO.EQ.1) WRITE(7,'(1X,99(A10,4X))') 'SCENARIO','YEAR    ','LOWER_CL ','MEDIAN  ','UPPER_CL ','AVERAGE ',  &
                                                     ' RUN_0  ','STD_DEV. '
            DO j=1,ncol ; write(7,100) scenario,j+display-1,low(j),median(j),hi(j),mean(j),X(j,1),std(j) ; END DO
         ENDIF
       ENDIF

      END DO ! Scenarios

      close(7)
100   format(4X,I4,9X,I4,6X,100(e13.6,1X))
101   format(13x,A15,1X,100(e13.6,1X))
102   format(13x,A1,I3,A11,1X,100(e13.6,1X))
103   format(13x,A8,I3,A4,1X,100(e13.6,1X))
104   format(13x,A7,I3,A4,2X,100(e13.6,1X))
      return
      end


!----------------------------------------------------------------------------
      SUBROUTINE SORT(N,RA)
!----------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER N,L,IR,I,J
      REAL (KIND=8) :: RA(N),RRA
      L=N/2+1 ; IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1 ; RRA=RA(L)
        ELSE
          RRA=RA(IR) ; RA(IR)=RA(1) ; IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RETURN
          ENDIF
        ENDIF
        I=L ; J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J) ; I=J ; J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
      GO TO 10
      END


!----------------------------------------------------------------------------
      SUBROUTINE ESTIMATE(OBJECTIVE,PARM,P,X,O,PB,NPARM,SEED,MAX_ITER)
!----------------------------------------------------------------------------
      
      IMPLICIT NONE
      INTEGER :: MFLAG
      PARAMETER(MFLAG=3)
      INTEGER :: I,ITER,J,MP,NDIM,NPARM,NTRY,GOOD,SEED,MAX_ITER,LOOPKILL
      REAL (KIND=8) :: PARM(NPARM+1),p(NPARM+1,NPARM),x(NPARM),O(NPARM+1),ftol,RSS,RNORM,ran,PDEV,PB(1)
      EXTERNAL OBJECTIVE
      MP=NPARM+1
      ndim=NPARM
      ntry=0
      FTOL=1.d-12
      PDEV=0.5d0
      LOOPKILL=MAX_ITER

! {The downhill simplex method must be started with NPARM+1 (MP) points--
!  defining an initial simplex (a total of NP*NP+NP values).  If you
!  think of one of these points as the initial point (P[1,j], j=1,np), the
!  second point (P[2,j], j=1,np) would be obtained adding to it the quantity
!  lambda[1] (which is the approximate length scale associate with parameter 1.
!  The third point would be obtained by adding lambda[2] to P[1,2] and so on.

210   DO  I=1,NPARM+1
        DO J=1,NPARM
          IF(j.NE.i-1) THEN
              P(i,j)=PARM(j)
            ELSE
!             {do random search in vicinity of current parameter estimates}
              P(i,j)=PARM(j)*EXP(PDEV*RNORM(SEED))
          ENDIF
      END DO ; END DO

      DO I=1,MP
          DO  j=1,NPARM ; x(j)=p(i,j) ; END DO
          CALL OBJECTIVE(X,O(i),NPARM)
      END DO

! {let amoeba do its thing}
211   CALL AMOEBA(p,O,mp,NPARM,ndim,ftol,OBJECTIVE,iter)

! {find the vertex with the lowest (best) objective function value}
      rss=1.d+32
      do 4 i=1,mp
        do 5 j=1,NPARM
          x(j)=p(i,j)
5       continue
        CALL OBJECTIVE(X,O(i),NPARM) 
        if(O(i).lt.rss) then
          rss=O(i)
          DO J=1,NPARM ; parm(J)=x(J) ; END DO
          PARM(NPARM+1)=RSS
        endif
4     continue

! {check for convergence}
      NTRY=NTRY+1
!      WRITE(*,*) '         Amoeba restart # ',NTRY,', RSS = ',RSS
!      WRITE(*,*) '         ',(PARM(I),I=1,NPARM)
      CALL CHECK(parm,ntry,good,NPARM,MFLAG)
      IF((GOOD.EQ.1).OR.(NTRY.GT.LOOPKILL)) then
        RETURN
      endif

! {restart the minimization routine at the point where it claims to have
!  found a minimum by setting P(1,J)=P(minimum RSS) and then reinitializing
!  the remaining NPARM number of simplex vertices}
      GO TO 210
      END
      
! ----------------------------------------------------
      SUBROUTINE CHECK(parmest,ntry,good,NP,NFLAG)
! check for convergence by comparing the results of the last 5 iterations
!----------------------------------------------------------------------------
      IMPLICIT NONE      
      INTEGER :: I,J,K,NP,NTRY,GOOD,NFLAG
      REAL (KIND=8) :: parmest(NP+1),parmold(10,401),r

      GOOD=0
      
      if(ntry.le.NFLAG) then
        DO J=1,NP+1 ; parmold(ntry,J)=parmest(J) ; END DO
        return
      endif

      DO J=1,NP+1
        DO I=1,NFLAG-1 ; parmold(I,J)=parmOLD(I+1,J) ; END DO
        PARMOLD(NFLAG,J)=PARMEST(J)
      END DO

!  Determine whether the value of the objective function varies by less
!  than 1% among the last NFLAG iterations.  If not, then return.  If so, then
!  apply the same criteria to the parameter estimates.
      DO K=NP+1,1,-1
        DO i=1,NFLAG
          DO J=i,NFLAG
              IF(PARMOLD(J,K).NE.0.d0) THEN
                  r=parmold(i,K)/parmold(J,K)
                  IF(R.LT.0.99d0.or.R.GT.1.01d0) return
              ENDIF
      END DO ; END DO ; END DO
      GOOD=1
      return
      end

