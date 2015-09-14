!-----------------------------------------------------------------
      SUBROUTINE READ_DATA  
! reads the data to be used in the objective function
!-----------------------------------------------------------------
      USE STATISTICS ; USE DATUM
      USE PARAMETERS, ONLY : INFILE,OPTION,SCALEs,MODEL_TYPE,PDF_TAG,NBOX,CV_OVERIDE,LINK_OLDEST,LINK_YOUNGEST
      IMPLICIT NONE
      REAL (KIND=8) :: TEMP(0:100)
      INTEGER :: S,G,Y,A,ITEMP(6),J,RECORD=0,STAT,COUNT,COUNTNEGS=0
      CHARACTER (LEN=10) :: CH
      CHARACTER (LEN=1) :: CH1
!     initialize all data (negative values mean no data)
      CATCH_DATA = 1.0 ; EFFORT_DATA = -9.0 ; SEL_EFFORT = 1.0 ; PSEL = -1.0 ; INDEX_MEAN=0 ; PDF_EFFORT=0 ; EFFORT_DATA_STORE=-1.0
      SIGMA_Q=0 ; SIGMA_F=0 ; SIGMA_M=0 ; SIGMA_T=0 ; SIGMA_TERMINAL=0 ; WEIGHT_CATCH=0 ; WEIGHT_SSB=0
      ETA_Q=0   ; ETA_F=0   ; ETA_M=0   ; ETA_T=0  ; ETA_TERMINAL=0
      SIGMA_TAG_REPORT=0 ; SIGMA_TAG_SURV=0 ; SIGMA_TAG_LOSS=0 ; ETA_TAG_REPORT=0 ; ETA_TAG_SURV=0 ; ETA_TAG_LOSS=0
      SIGMA_TAG_NOMIX2=0 ; SIGMA_TAG_NOMIX=0 ; ETA_TAG_NOMIX=0 ;  N_Qs=0
      OPEN(10,FILE=INFILE(1),STATUS='OLD',ERR=999)
      IF(MODEL_TYPE==0) THEN
          IF(OPTION(1)/=2) THEN
            OPTION(1)=2 ! Combine data from two areas for 1 vpa
            WRITE(*,*) 'WARNING: The model type was specified as 0, which is the option to  '
            WRITE(*,*) '         combine 2 stocks/areas, but something other than 2 stocks/areas'
            WRITE(*,*) '         was indicated on the preceding line.'
            WRITE(*,*) '         The program will now assume there are data entered for 2 stocks.'
            WRITE(*,*) '         Do you wish to continue? Type Yes or No: '; READ(*,*) CH1 ; IF(CH1=='N' .OR. CH1=='n') STOP
          END IF
        ELSE IF(MODEL_TYPE==2 .AND. OPTION(1)==1) THEN
          MODEL_TYPE=1
        ELSE IF(MODEL_TYPE<0 .or. MODEL_TYPE>2) THEN
            WRITE(*,*) 'WARNING:   The model type was specified as = ',MODEL_TYPE
            WRITE(*,*) '           The only options available are 0, 1, or 2.'
            WRITE(*,*) '           The program will now assume you meant to do a single stock analysis.'
            WRITE(*,*) '         Do you wish to continue? Type Yes or No: '; READ(*,*) CH1 ; IF(CH1=='N' .OR. CH1=='n') STOP
      END IF
      NBOX=OPTION(1);
      IF(NBOX<1) THEN
          WRITE(*,*) 'WARNING: The number of stocks/areas was specified as ',NBOX
          WRITE(*,*) '         The program will now assume you meant 1 stock/area.'
          WRITE(*,*) '         Do you wish to continue? Type Yes or No: '; READ(*,*) CH1 ; IF(CH1=='N' .OR. CH1=='n') STOP
          NBOX=1
      ENDIF
!     read general specifications
      WRITE(*,*) 'READING GENERAL SPECIFICATIONS IN DATA FILE'
      CALL BLANK(10,RECORD)
      READ(10,*,IOSTAT=STAT,ERR=1000,END=1000) FIRSTYEAR,LASTYEAR ; IF(STAT.NE.0) GOTO 1000   
      DISPLAYYEAR=FIRSTYEAR ; LASTYEAR=LASTYEAR-FIRSTYEAR+1 ; FIRSTYEAR=1 ; NYEARS=LASTYEAR
      IGNORE_RECRUIT=IGNORE_RECRUIT-DISPLAYYEAR+1
      CALL BLANK(10,RECORD)
      READ(10,*,IOSTAT=STAT,ERR=1000,END=1000) FIRSTAGE,LASTAGE,PLUSGROUP,EXPANDED_PLUSGROUP
      IF(STAT.NE.0) GOTO 1000   
      IF(PLUSGROUP.GE.0) PLUSGROUP=MIN(LASTAGE,PLUSGROUP); IF(EXPANDED_PLUSGROUP<PLUSGROUP) EXPANDED_PLUSGROUP=PLUSGROUP
      IF(LINK_OLDEST>LASTAGE) LINK_OLDEST=LASTAGE ; IF(LINK_YOUNGEST<FIRSTAGE) LINK_YOUNGEST=FIRSTAGE

!     read area-specific data
      DO 1 BOX=1,NBOX
        CALL BLANK(10,RECORD)
        READ(10,*,IOSTAT=STAT,ERR=1000,END=1000) NGEARS(box) 
        IF(STAT.NE.0) GOTO 1000   
        CALL BLANK(10,RECORD)
        READ(10,*) SEASON_SSB(box)
        CALL BLANK(10,RECORD)
        READ(10,*) (FECUNDITY(box,a),a=FIRSTAGE,EXPANDED_PLUSGROUP)
!       read catch data
        WRITE(*,*) 'READING CATCH DATA' ; COUNT=0
        CALL BLANK(10,RECORD)
        READ(10,*) TITLE_FISHERY(box),PDF_CATCH(box),SIGMA_CATCH(box)
10      CALL BLANK(10,RECORD)
        READ(10,*) Y  
        IF(Y.LE.0) THEN ; GOTO 19; ELSE ; Y=Y-DISPLAYYEAR+1 ; CALL CHECKDIMENSION(Y,Ys,LASTYEAR,'YEAR      ',RECORD) ; ENDIF
        BACKSPACE(10)
        READ(10,*) S,(CATCH_DATA(box,a,y), a=FIRSTAGE,LASTAGE)  ; COUNT=COUNT+1
        DO A=FIRSTAGE,LASTAGE ; IF(CATCH_DATA(box,a,y)<=0.0) CATCH_DATA(box,a,y)=1.0 ; END DO
        GOTO 10

19      IF(COUNT/=LASTYEAR) THEN
          WRITE(*,*) 'WARNING: Number of rows with catch data does not equal the number of years.'
          WRITE(*,*) '         The program automatically assigned the missing catch data value of 1.0'
          WRITE(*,*) '         Do you wish to continue? Type Yes or No: '; READ(*,*) CH1 ; IF(CH1=='N' .OR. CH1=='n') STOP
        ENDIF
!       read effort (cpue) specifications
        WRITE(*,*) 'READING CPUE DATA' ; CH='FISHERY  '
30      CALL BLANK(10,RECORD)
        READ(10,*) G
        IF(G.LE.0) THEN ; COUNTNEGS=COUNTNEGS+1 ; GOTO 40 ; ELSE ; CALL CHECKDIMENSION(G,Gs,NGEARS(BOX),CH,RECORD) ; ENDIF
        BACKSPACE(10)
        READ(10,*) G,ITEMP,TITLE_EFFORT(g,box)
        WRITE(*,*) g,ITEMP,TITLE_EFFORT(g,box)
        PDF_EFFORT(g,box)=ITEMP(1) ; BIO_EFFORT(g,box)=abs(ITEMP(2)) ; SEL_TYPE(g,box)=ITEMP(3)
        SEASON_EFFORT(g,box)=iTEMP(4) ; AGE_EFFORT(1,g,box)=ITEMP(5) ; AGE_EFFORT(2,g,box)=ITEMP(6) 
        IF(PDF_EFFORT(g,box)/=0.AND.OPTION(4)<=0.AND.PDF_EFFORT(g,box)/=12.AND.BIO_EFFORT(g,box)<7) N_Qs=N_Qs+1  ! number of Q's estimated in conc. like.
        CALL CHECKDIMENSION(AGE_EFFORT(1,g,box),As,EXPANDED_PLUSGROUP,'AGE       ',RECORD)
        CALL CHECKDIMENSION(AGE_EFFORT(2,g,box),As,EXPANDED_PLUSGROUP,'AGE       ',RECORD)
!       hardwire
        IF(ITEMP(2)<0 .AND. NBOX>1) THEN ; BOTH_STOCKS(g,box)=1 ; ELSE ; BOTH_STOCKS(g,box)=-1 ; ENDIF
        IF(PDF_EFFORT(g,box) <= 0) THEN
            WRITE(*,*) 'YOU HAVE CHOSEN NOT TO USE THE INDEX FOR SURVEY/FISHERY ',G
          ELSE IF( (BIO_EFFORT(g,box)==7 .OR. BIO_EFFORT(g,box)==8) .AND. MODEL_TYPE==1 ) THEN
            WRITE(*,*) 'ERROR: PROPORTIONS CAN ONLY BE ESTIMATED WITH THE OVERLAP MODEL (2). CHECK SPECIFICATIONS FOR INDEX ',G  ;  STOP
        ENDIF
        IF((SEL_TYPE(g,box)>1).AND.(AGE_EFFORT(1,g,box)>PLUSGROUP)) THEN
          WRITE(*,*) 'ERROR: SELECTIVITY MUST BE FIXED TO 1.0 WHEN THE INDEX APPLIES ONLY TO AGES OLDER THAN THE PLUS GROUP.'
          write(*,*) '       CHECK SPECIFICATIONS FOR INDEX ',G  ;  STOP
        ENDIF
        GOTO 30

!       read effort abundance index data
40      CALL BLANK(10,RECORD)
        READ(10,*) G
        IF(G.LE.0) THEN ; COUNTNEGS=COUNTNEGS+1 ; GOTO 41; ELSE ; CALL CHECKDIMENSION(G,Gs,NGEARS(BOX),'FISHERY   ',RECORD) ; ENDIF
        BACKSPACE(10)
        READ(10,*) G,Y,TEMP(1),TEMP(2)
        Y=Y-DISPLAYYEAR+1 ; CALL CHECKDIMENSION(Y,Ys,LASTYEAR,'YEAR      ',RECORD)
        EFFORT_DATA(g,box,y) = TEMP(1) ; SIGMA_EFFORT(g,box,y) = TEMP(2) ; EFFORT_DATA_STORE(g,box,y)=TEMP(1)
        IF(PDF_EFFORT(g,box)==1 .and. BIO_EFFORT(g,box)>2 .AND. BIO_EFFORT(g,box)<7) EFFORT_DATA(g,box,y)=EFFORT_DATA(g,box,y)+0.00001
        GOTO 40
     
41      IF((CV_OVERIDE>0.00001 .AND. INT(CV_OVERIDE)/=999) .or. CV_OVERIDE<-0.00001) SIGMA_EFFORT=CV_OVERIDE
!       read effort selectivity data
        CALL BLANK(10,RECORD)
        READ(10,*) G
        IF(G.LE.0) THEN ; COUNTNEGS=COUNTNEGS+1 ; GOTO 49; ELSE ; CALL CHECKDIMENSION(G,Gs,NGEARS(BOX),'FISHERY   ',RECORD) ; ENDIF
        BACKSPACE(10)
        READ(10,*) G,Y
        Y=Y-DISPLAYYEAR+1 ; CALL CHECKDIMENSION(Y,Ys,LASTYEAR,'YEAR      ',RECORD)
        BACKSPACE(10)
        READ(10,*) G,j,(PSEL(g,box,a,y),a=FIRSTAGE,LASTAGE)
        PSEL_MEAN(g,box,y)=0.0
        DO 42 A=AGE_EFFORT(1,g,box),LASTAGE
          IF(SEL_TYPE(g,box)==1) THEN ! fixed selectivities
              SEL_EFFORT(g,box,a,y)=PSEL(g,box,a,y)
            ELSEIF(SEL_TYPE(g,box)>=3) THEN ! partial catches
              PSEL(g,box,a,y)=PSEL(g,box,a,y)/CATCH_DATA(box,a,y)
          ENDIF
          PSEL_MEAN(g,box,y)=PSEL_MEAN(g,box,y)+PSEL(g,box,a,y)/DBLE(LASTAGE-AGE_EFFORT(1,g,box)+1)
42      CONTINUE
        DO A=PLUSGROUP,LASTAGE; IF(PSEL(g,box,a,y)<=0.0) PSEL(g,box,a,y)=PSEL_MEAN(g,box,y)*1.0D-3; END DO
        GOTO 41

!       read weight at age data 
49      WRITE(*,*) 'READING WEIGHT AT AGE DATA FOR FISHERIES AND SURVEYS'
50      CALL BLANK(10,RECORD)
        READ(10,*) G
        IF(G.LE.0) THEN ; COUNTNEGS=COUNTNEGS+1 ; GOTO 59; ELSE ; CALL CHECKDIMENSION(G,Gs,NGEARS(box),'INDEX     ',RECORD) ; ENDIF
        BACKSPACE(10)
        READ(10,*) G,Y,(TEMP(a), a=FIRSTAGE,LASTAGE) ; Y=Y-DISPLAYYEAR+1 ; CALL CHECKDIMENSION(Y,Ys,LASTYEAR,'YEAR      ',RECORD)
        DO A=FIRSTAGE,LASTAGE
          IF(A>=PLUSGROUP .and. PLUSGROUP<LASTAGE) THEN
              IF(SEL_TYPE(g,box)>=2 .AND. PSEL(g,box,a,y)>0) THEN
                  WEIGHT_CATCH(g,box,PLUSGROUP,y)=WEIGHT_CATCH(g,box,PLUSGROUP,y)+TEMP(a)*PSEL(g,box,a,y)
                ELSE
                  WEIGHT_CATCH(g,box,PLUSGROUP,y)=WEIGHT_CATCH(g,box,PLUSGROUP,y)+TEMP(a)*CATCH_DATA(box,a,y)
              ENDIF
            ELSE
              WEIGHT_CATCH(g,box,a,y)=TEMP(a)
          ENDIF
        END DO
        GOTO 50

59      WRITE(*,*) 'READING WEIGHT AT AGE DATA FOR SPAWNING STOCK BIOMASS'
60      CALL BLANK(10,RECORD)
        READ(10,*) Y ; Y=Y-DISPLAYYEAR+1
        IF(Y.LE.0) THEN ; COUNTNEGS=COUNTNEGS+1 ; GOTO 69; ELSE ; CALL CHECKDIMENSION(Y,Ys,LASTYEAR,'YEAR      ',RECORD) ; ENDIF
        BACKSPACE(10)
        READ(10,*) j,(TEMP(a), a=FIRSTAGE,LASTAGE)
        DO A=FIRSTAGE,LASTAGE
          IF(A>=PLUSGROUP .and. PLUSGROUP<LASTAGE) THEN
              WEIGHT_SSB(box,PLUSGROUP,y)=WEIGHT_SSB(box,PLUSGROUP,y)+TEMP(a)*CATCH_DATA(box,a,y)
            ELSE
              WEIGHT_SSB(box,a,y)=TEMP(a)
          ENDIF
        END DO
        GOTO 60

!       combine catches and weights into plus group
69      DO Y=1,LASTYEAR
          DO A=PLUSGROUP+1,LASTAGE
            CATCH_DATA(box,PLUSGROUP,y)=CATCH_DATA(box,PLUSGROUP,y)+CATCH_DATA(box,a,y)
            DO G=1,NGEARS(box)
              IF(SEL_TYPE(g,box)>=2) PSEL(g,box,PLUSGROUP,y)=PSEL(g,box,PLUSGROUP,y)+PSEL(g,box,a,y)
            END DO
          END DO
          IF(PLUSGROUP<LASTAGE) THEN
            DO G=1,NGEARS(box)
              IF(SEL_TYPE(g,box)>=2 .and. PSEL(g,box,PLUSGROUP,y)>0) THEN ;
                  WEIGHT_CATCH(g,box,PLUSGROUP,y)=WEIGHT_CATCH(g,box,PLUSGROUP,y)/PSEL(g,box,PLUSGROUP,y)
                ELSE
                  WEIGHT_CATCH(g,box,PLUSGROUP,y)=WEIGHT_CATCH(g,box,PLUSGROUP,y)/CATCH_DATA(box,PLUSGROUP,y)
              ENDIF
            END DO
            WEIGHT_SSB(box,PLUSGROUP,y)=WEIGHT_SSB(box,PLUSGROUP,y)/CATCH_DATA(box,PLUSGROUP,y)
          ENDIF
        END DO

!       scale indices by their mean
        IF(SCALEs.GT.0) THEN ! Scale indices by their mean (geometric in case of lognormal model)
          write(*,*) 'Scaling indices of abundance by their mean'
          DO G=1,NGEARS(box)
            IF(PDF_EFFORT(g,box)==12) exit
            COUNT=0
            DO 83 Y=1,LASTYEAR
              IF(PDF_EFFORT(g,box)>1 .and. EFFORT_DATA(g,box,y)>=0) THEN
                  INDEX_MEAN(g,box)=INDEX_MEAN(g,box)+EFFORT_DATA(g,box,y) ; COUNT=COUNT+1
                ELSEIF(PDF_EFFORT(g,box)==1 .and. EFFORT_DATA(g,box,y)>0) THEN
                  INDEX_MEAN(g,box)=INDEX_MEAN(g,box)+dLOG(EFFORT_DATA(g,box,y)) ; COUNT=COUNT+1
              ENDIF
83          CONTINUE
            INDEX_MEAN(g,box)=INDEX_MEAN(g,box)/COUNT
            IF(PDF_EFFORT(g,box)==1) INDEX_MEAN(g,box)=DEXP(INDEX_MEAN(g,box))
            IF(BIO_EFFORT(g,box)>6) INDEX_MEAN(g,box)=1.0
            DO 84 Y=1,LASTYEAR
              IF(EFFORT_DATA(g,box,y)>0 .OR. (PDF_EFFORT(g,box)>1 .AND. EFFORT_DATA(g,box,y)>=0)) THEN
                EFFORT_DATA(g,box,y)=EFFORT_DATA(g,box,y)/INDEX_MEAN(g,box)  
                EFFORT_DATA_STORE(g,box,y)=EFFORT_DATA(g,box,y)/INDEX_MEAN(g,box)
              ENDIF
84          CONTINUE
          END DO
         ELSE ! Do not scale by the mean
          INDEX_MEAN=1
        ENDIF

      IF(BOX==2 .AND. MODEL_TYPE==0) THEN ! Combine statistics for both areas  
        DO 100 GEAR=1,NGEARS(2)
          G=GEAR+NGEARS(1) ; IF(G>Gs) THEN ; WRITE(*,*) 'Dimensions for indices too small',GS ; PAUSE; STOP ; ENDIF
          PDF_EFFORT(g,1) = PDF_EFFORT(gear,2) ; BIO_EFFORT(g,1)=BIO_EFFORT(gear,2) ; SEL_TYPE(g,1)=SEL_TYPE(gear,2)
          SEASON_EFFORT(g,1)=SEASON_EFFORT(gear,2) ; AGE_EFFORT(1,g,1)=AGE_EFFORT(1,gear,2) ; AGE_EFFORT(2,g,1)=AGE_EFFORT(2,gear,2)
          INDEX_MEAN(g,1)=INDEX_MEAN(gear,2) ; TITLE_EFFORT(g,1)=TITLE_EFFORT(gear,2)
          DO 110 Y=FIRSTYEAR,LASTYEAR
            EFFORT_DATA(g,1,y)=EFFORT_DATA(gear,2,y) ; SIGMA_EFFORT(g,1,y)=SIGMA_EFFORT(g,2,y)
            EFFORT_DATA_STORE(g,1,y)=EFFORT_DATA_STORE(gear,2,Y)
            DO 120 A=FIRSTAGE,LASTAGE
              SEL_EFFORT(g,1,a,y)=SEL_EFFORT(gear,2,a,y) ; PSEL(g,1,a,y)=PSEL(gear,2,a,y)
              WEIGHT_CATCH(g,1,a,y)=WEIGHT_CATCH(gear,2,a,y)              
              IF(GEAR==1) CATCH_DATA(1,A,Y)=CATCH_DATA(1,A,Y)+CATCH_DATA(2,A,Y)
120         CONTINUE
110       CONTINUE
100     CONTINUE
        NGEARS(1)=NGEARS(1)+NGEARS(2) ; MODEL_TYPE=1 ; NBOX=1
      ENDIF

1     CONTINUE
      LASTAGE=PLUSGROUP
      WRITE(*,*) 'FINISHED READING DATA FILE' ; WRITE(*,'(i5,a21)') RECORD,' RECORDS READ IN FILE'
      IF(COUNTNEGS<5) THEN ; WRITE(*,*) 'INPUT ERROR: Too few categories in data file (you are missing some -1)';  ENDIF
      IF(COUNTNEGS<5) THEN ; WRITE(12,*) 'INPUT ERROR: Too few categories in data file (you are missing some -1)'; PAUSE;STOP; ENDIF
      IF(PDF_TAG>0) CALL READ_TAG ; CLOSE(10) ; RETURN

999   WRITE(*,'(1X,A21,A120)') 'ERROR: The data file ',INFILE(1)
      WRITE(*,'(8X,A55)') 'does not exist or is being used by another application.'
      WRITE(12,'(1X,A21,A120)') 'ERROR: The data file ',INFILE(1)
      WRITE(12,'(8X,A55)') 'does not exist or is being used by another application.'  ; pause ; STOP
1000  CALL MISTAKE(STAT,RECORD)
      END

!-----------------------------------------------------------------
      SUBROUTINE READ_TAG
!-----------------------------------------------------------------
      USE STATISTICS ; USE DATUM
      USE PARAMETERS, ONLY : INFILE,OPTION,SCALEs,MODEL_TYPE,NBOX,SIGMA_TAG_MOD
      IMPLICIT NONE
      INTEGER :: RECORD=0,J,YREC,BREC,Y_end
      REAL (KIND=8) :: SUM_TAG_FISH,SUM_FISH=0
      CHARACTER (LEN=3) :: tag_format_type='new'
      OPEN(10,FILE=INFILE(6),STATUS='OLD',ERR=999)
      WRITE(*,*) 'READING TAG COHORT DEFINITIONS'
      SUMRDATA=0 ; N_Tag_cohorts=0; RECAPTURE_DATA=-1; RECOVERY_DATA=-1;
1     CALL BLANK(10,RECORD)
      READ(10,*,END=2) J ; IF(J<=0) THEN ; GOTO 2 ; ELSE ; BACKSPACE(10) ; ENDIF
      READ(10,*,END=2) J,TC(j) ;
      TC(j)%Yr=TC(j)%Yr-DISPLAYYEAR+1; TC(j)%Ye=TC(j)%Ye-DISPLAYYEAR+1
      T_mod_month(j)=int(TC(j)%Mr)+1; if(T_mod_month(j)>12) T_mod_month(j)=12;
      TC(j)%Mr=TC(j)%Mr/12.0; TC(j)%Me=TC(j)%Me/12.0
      TC(j)%Sigma=TC(j)%Sigma/sqrt(SIGMA_TAG_MOD)
      IF(J>N_TAG_COHORTS) N_TAG_COHORTS=J;
      IF(TC(j)%Surv>N_SURV_PARS) N_SURV_PARS=TC(j)%Surv;
      IF(TC(j)%Shed>N_SHED_PARS) N_SHED_PARS=TC(j)%Shed;
      IF(TC(j)%Rpt0>N_REPT_PARS) N_REPT_PARS=TC(j)%Rpt0;
      IF(TC(j)%Rpt1>N_REPT_PARS) N_REPT_PARS=TC(j)%Rpt1;
      IF(TC(j)%Rpt2>N_REPT_PARS) N_REPT_PARS=TC(j)%Rpt2;
      IF(TC(j)%Mix1>N_MIX_PARS)  N_MIX_PARS=TC(j)%Mix1;
      IF(TC(j)%Mix2>N_MIX_PARS)  N_MIX_PARS=TC(j)%Mix2;

      GOTO 1

2     WRITE(*,*) 'READING TAG RECOVERY DATA'
      IF(J<=-2) THEN; WRITE(*,*) ' (using the old format with entries for all years)' ; tag_format_type='old' ; ENDIF

3     CALL BLANK(10,RECORD)
      READ(10,*,END=100) J ; IF(J<=0) THEN ; GOTO 100 ; ELSE ; BACKSPACE(10) ; ENDIF
      IF(J>N_tag_cohorts) THEN; WRITE(*,*) 'ERROR: Tag cohort ',j,' was not defined above.'; STOP; ENDIF
      IF(TC(j)%N<=0.0) THEN
          ! Ignore recovery data
          READ(10,*,END=100) J
        ELSE
          IF(tag_format_type=='new') THEN;
              Y_END = TC(j)%Ye; IF(Y_END>LASTYEAR) Y_END=LASTYEAR;
              READ(10,*,END=100) J,Brec,RECOVERY_DATA(j,brec),(RECAPTURE_DATA(j,brec,yrec),YREC=TC(j)%Yr,Y_end)
            ELSE
              READ(10,*,END=100) J,Brec,(RECAPTURE_DATA(j,brec,yrec),YREC=FIRSTYEAR,LASTYEAR)
              RECOVERY_DATA(j,brec)=-1 ! old format cannot be used for electronic data
          ENDIF
      ENDIF
      GOTO 3

100   WRITE(*,*) 'END OF TAGGING FILE REACHED'
      DO J=1,N_TAG_COHORTS
        IF(TC(j)%N<=0.0) CYCLE
        DO BREC=1,NBOX
          DO YREC=TC(j)%Yr,TC(j)%Ye
            IF(YREC<=LASTYEAR .AND. RECAPTURE_DATA(j,brec,yrec)>0) SUMRDATA(j)=SUMRDATA(j)+RECAPTURE_DATA(j,brec,yrec);
          END DO
          IF(TC(j)%Ye>LASTYEAR) RECOVERY_DATA(j,brec)=-1
          IF(RECOVERY_DATA(j,brec)>0) SUMRDATA(j)=SUMRDATA(j)+RECOVERY_DATA(j,brec);
        END DO
        IF(TC(j)%Ye>LASTYEAR) THEN ; TC(j)%Ye=LASTYEAR; TC(j)%Me=1.0; ENDIF  
        IF(SUMRDATA(j)>TC(j)%N) THEN
          WRITE(*,*) 'ERROR: There are more recoveries than releases for tag cohort',j
          PAUSE; STOP
        ENDIF
      END DO
      CLOSE(10) ; RETURN

999   WRITE(*,'(1X,A25,A120)') 'ERROR: The tag data file ',INFILE(6)
      WRITE(*,'(8X,A55)') 'does not exist or is being used by another application.'
      WRITE(12,'(1X,A25,A120)') 'ERROR: The tag data file ',INFILE(6)
      WRITE(12,'(8X,A55)') 'does not exist or is being used by another application.'  ; PAUSE ; STOP
      END

!------------------------------------------------------------------------------
      SUBROUTINE OBJECTIVE(X,FUNC,NPARM)
!------------------------------------------------------------------------------
      USE PARAMETERS, ONLY : ICALL,FIRSTCALL
      USE LOGLIKELIHOODS, ONLY : PENALTY
      IMPLICIT NONE
      INTEGER NPARM
      REAL (KIND=8) :: X(NPARM),LIKE,FUNC
      PENALTY=0
      IF(NPARM==0) GOTO 1
      CALL DEFINE_PARAMETERS(X,NPARM)
      CALL MODEL
1     CALL LIKELIHOOD(LIKE,X,NPARM)
      FUNC = -1.d0*LIKE+PENALTY
      ICALL=ICALL+1 ; IF(FIRSTCALL==0) FIRSTCALL=1
      RETURN ; END

!------------------------------------------------------------------------------
      SUBROUTINE DEFINE_PARAMETERS(X,NPARM)
! Defines the parameters read by READ_PARAMETERS and tells how they should be estimated
!------------------------------------------------------------------------------
      USE STATISTICS ; USE PARAMETERS
      USE DATUM
      IMPLICIT NONE
      INTEGER I,J,K,A,Y,NPARM,G
      REAL (KIND=8) :: X(NPARM), TEMP, DUMMY, SURV_PAR(0:1000),LOSS_PAR(0:1000),REPT_PAR(0:1000),NMIX_PAR(0:1000)
      K=0 ; I=1 ! (I gives the next variable, K gives the current estimated parameter)
      NPI=0

      DO 1 BOX=1,NBOX
! {terminal-year age structure of population}
       DO 10 A=FIRSTAGE,LASTAGE-1
        IF(ICALL==0) THEN
            CALL GETPAR(I,K,X,NPARM,TEMP,ETA_Terminal(box,a),SIGMA_Terminal(box,a))
            IF(PARM_SPECS(I-1,4)>0 .AND. PARM_SPECS(I-1,4)<0.3) GOTO 991 ! no random deviations models
            IF(PARM_SPECS(I-1,4)==0 .AND. TEMP<10.0) THEN ! treat as a relative selectivity
                SEL_TERMINAL(box,a,1)=TEMP ; SEL_TERMINAL(box,a,2) = PARM_SPECS(I-1,5)
              ELSE 
                IF(OPTION(2)==2) THEN ! parameter measures abundance
                    IF(MODEL_TYPE==1) THEN ; N_AREA(box,a+1,lastyear+1)=TEMP ; ELSE ; N_STOCK(box,a+1,lastyear+1)=TEMP ; ENDIF
                  ELSE ! parameter measures fishing mortality
                    F(box,a,lastyear)=TEMP 
                ENDIF
                SEL_TERMINAL(box,a,2) = -1 
            ENDIF
          ELSEIF(SEL_TERMINAL(box,a,2)<0) THEN
            IF(OPTION(2)==2) THEN
                IF(MODEL_TYPE==1) THEN
                    CALL GETPAR(I,K,X,NPARM,N_AREA(box,a+1,lastyear+1),ETA_Terminal(box,a),SIGMA_Terminal(box,a))
                  ELSE
                    CALL GETPAR(I,K,X,NPARM,N_STOCK(box,a+1,lastyear+1),ETA_Terminal(box,a),SIGMA_Terminal(box,a))
                ENDIF
              ELSE
                  CALL GETPAR(I,K,X,NPARM,F(box,a,lastyear),ETA_Terminal(box,a),SIGMA_Terminal(box,a))
            ENDIF
          ELSE
            CALL GETPAR(I,K,X,NPARM,TEMP,ETA_Terminal(box,a),SIGMA_Terminal(box,a))
        ENDIF
10     CONTINUE
       SEL_TERMINAL(box,lastage,2) = -1

! {terminal-age f-ratios)
       DO 20 Y=FIRSTYEAR,LASTYEAR
         IF(FIRSTCALL==0.and.Y==FIRSTYEAR.AND.((PARM_SPECS(I,4)<0.and.PARM_SPECS(I,4)>-1.0) .or. &
           (PARM_SPECS(I,4)>0.0001 .AND. PARM_SPECS(I,4)<0.299))) GOTO 992
         CALL GETPAR(I,K,X,NPARM,F_RATIO(box,y),ETA_F(box,y),SIGMA_F(box,y))
20     CONTINUE

! {natural mortality}
       DO A=FIRSTAGE,LASTAGE
         IF(FIRSTCALL==0.and.a==FIRSTAGE.AND.((PARM_SPECS(I,4)<0.and.PARM_SPECS(I,4)>-1.0) .or. &
           (PARM_SPECS(I,4)>0.0001 .AND. PARM_SPECS(I,4)<0.299))) GOTO 993
         CALL GETPAR(I,K,X,NPARM,M(box,a,firstyear),ETA_M(box,a,firstyear),SIGMA_M(box,a,firstyear))
         DO Y=FIRSTYEAR+1,LASTYEAR ; M(box,a,y)=M(box,a,firstyear) ; END DO
       END DO

! {transfer coefficients)
       DO A=FIRSTAGE,LASTAGE
         IF(FIRSTCALL==0.and.a==FIRSTAGE.AND.((PARM_SPECS(I,4)<0.and.PARM_SPECS(I,4)>-1.0) .or. &
           (PARM_SPECS(I,4)>0.0001 .AND. PARM_SPECS(I,4)<0.299))) GOTO 994
         CALL GETPAR(I,K,X,NPARM,T(box,a,firstyear),ETA_T(box,a,firstyear),SIGMA_T(box,a,firstyear))
         DO Y=FIRSTYEAR+1,LASTYEAR ; T(box,a,y)=T(box,a,firstyear) ; END DO
       END DO

! {Stock recruitment relationship}
       DO A=1,5
         IF(FIRSTCALL==0) THEN
           IF(PDF_STOCKRECRUIT==0 .AND. (PARM_SPECS(I,4)>0.0001 .OR. PARM_SPECS(I,4)<0)) THEN
             WRITE(*,*) 'ERROR: The stock-recruit penalty in the control file: ',CONFILE
             WRITE(*,*) '       has been turned off (set to 0), but the stock-recruit parameters are'
             WRITE(*,*) '       still being estimated (or linked). Either set the method indicators'
             WRITE(*,*) '       for these parameters to 0 or activate the stock-recruit penalty'
             WRITE(12,*) 'ERROR: The stock-recruit penalty in the control file: ',CONFILE
             WRITE(12,*) '       has been turned off (set to 0), but the stock-recruit parameters are'
             WRITE(12,*) '       still being estimated (or linked). Either set the method indicators'
             WRITE(12,*) '       for these parameters to 0 or activate the stock-recruit penalty' ; PAUSE; STOP
            ELSEIF(PDF_STOCKRECRUIT /=0 .AND. A<3 .AND. PARM_SPECS(I,4)<0.001 .and. PARM_SPECS(I,4)>-.0001) THEN
             WRITE(*,*) 'WARNING: The stock-recruit penalty in the control file: ',CONFILE
             WRITE(*,*) '       has been turned on, but the stock-recruit parameters are not'
             WRITE(*,*) '       being estimated. Either set the method indicators for these'
             WRITE(*,*) '       parameters to 1 or deactivate the stock-recruit penalty'
             WRITE(*,*) '       Press enter to continue, Control-C to abort.' ; READ(*,*)
             WRITE(12,*) 'WARNING: The stock-recruit penalty in the control file: ',CONFILE
             WRITE(12,*) '       has been turned on, but the stock-recruit parameters are not'
             WRITE(12,*) '       being estimated. Either set the method indicators for these'
             WRITE(12,*) '       parameters to 1 or deactivate the stock-recruit penalty'
             WRITE(12,*) '       Press enter to continue, Control-C to abort.' ; READ(*,*)
            ELSEIF(PDF_STOCKRECRUIT <0 .AND. A==5 .AND. PARM_SPECS(I,4)>=0.001) THEN
             WRITE(*,*) 'ERROR: The stock-recruit penalty in the control file: ',CONFILE
             WRITE(*,*) '       has been turned on with a negative value, so you are trying to'
             WRITE(*,*) '       estimate the standard error term both by the concentrated method'
             WRITE(*,*) '       and via the parameter file. Either Either set the method indicator'
             WRITE(*,*) '       for the 5th stock-recruit parameter to 0 or inactivate the'
             WRITE(*,*) '       concentrated method by making the control file entry positive.'
             WRITE(12,*) 'ERROR: The stock-recruit penalty in the control file: ',CONFILE
             WRITE(12,*) '       has been turned on with a negative value, so you are trying to'
             WRITE(12,*) '       estimate the standard error term both by the concentrated method'
             WRITE(12,*) '       and via the parameter file. Either Either set the method indicator'
             WRITE(12,*) '       for the 5th stock-recruit parameter to 0 or inactivate the'
             WRITE(12,*) '       concentrated method by making the control file entry positive.' ; PAUSE; STOP
            ELSEIF((PARM_SPECS(I,4)>0.0001 .AND. PARM_SPECS(I,4)<.3) .OR. PARM_SPECS(I,4)<0) THEN
             WRITE(*,*) 'ERROR: Stock-recruit parameters may not be linked to one another or to transfer rates'
             WRITE(12,*) 'ERROR: Stock-recruit parameters may not be linked to one another or to transfer rates'  ; PAUSE; STOP
           ENDIF
         ENDIF
         CALL GETPAR(I,K,X,NPARM,STOCK_RECRUIT(box,a),DUMMY,DUMMY)
       END DO

! {Variance scaling parameter}
       TEMP=0
       DO G=1,NGEARS(box)
         IF(FIRSTCALL==0) THEN
           IF(PARM_SPECS(I,4)>0) THEN ! variance parameter being estimated
               TEMP=TEMP+1; NPI(g,box)=NPI(g,box)+1
           ENDIF
           IF(INT(CV_OVERIDE)==999 .AND. (PARM_SPECS(I,4)>0.0001 .OR. PARM_SPECS(I,4)<0)) THEN
             WRITE(*,*) 'ERROR: The concentrated-likelihood (iterative-reweighting) variance option was'
             WRITE(*,*) '       selected in the control file: ',CONFILE
             WRITE(*,*) '       The variance scaling parameters in file ',INFILE(2),' are therefore redundant'
             WRITE(*,*) '       and must be set to a constant (say 1.0) with the method indicator set to 0.'
             WRITE(12,*) 'ERROR: The concentrated-likelihood (iterative-reweighting) variance option was'
             WRITE(12,*) '       selected in the control file: ',CONFILE
             WRITE(12,*) '       The variance scaling parameters in file ',INFILE(2),' are therefore redundant'
             WRITE(12,*) '       and must be set to a constant (say 1.0) with the method indicator set to 0.' ; PAUSE; STOP
            ELSEIF(PDF_EFFORT(G,box)==0 .AND. PARM_SPECS(I,4)>0.0001) THEN
             WRITE(*,*) 'ERROR: Index',g,'is not being used, so the corresponding variance parameter'
             WRITE(*,*) '       cannot be estimated and the method indicator should be set to 0 or -0.1).'
             WRITE(12,*) 'ERROR: Index',g,'is not being used, so the corresponding variance parameter'
             WRITE(12,*) '       cannot be estimated and the method indicator should be set to 0 or -0.1).'   ; PAUSE; STOP
            ELSEIF(G==1.AND.((PARM_SPECS(I,4)>0.0001 .AND. PARM_SPECS(I,4)<.3) .OR. PARM_SPECS(I,4)<0)) THEN ; GOTO 1000
           ENDIF
         ENDIF
         CALL GETPAR(I,K,X,NPARM,SCALE_VARIANCE(g,box),ETA_SCALE_VARIANCE(g,box),SIGMA_SCALE_VARIANCE(g,box))
       END DO
       IF(INT(TEMP+0.5)==1) NPI(1,box)=NPI(1,box)-1 ! only one overall variance was estimated: not unique to that estimate.

! {catchability coefficients for indices)
       IF(OPTION(4)>0) THEN ! Estimate catchabilities in search
         DO G=1,NGEARS(box)
          DO Y=FIRSTYEAR,LASTYEAR
            IF(ICALL==0) THEN
              IF(INT(PARM_SPECS(I,4))==1) NPI(g,box)=NPI(g,box)+1 ! Bayesian types not counted because 1-npi/ndi would be >=0
              IF(Y==FIRSTYEAR.AND.(PARM_SPECS(I,4)<0 .or. (PARM_SPECS(I,4)>0.0001 .AND. PARM_SPECS(I,4)<0.3))) THEN ; GOTO 995
               ELSEIF(PDF_EFFORT(G,box)==0 .AND. PARM_SPECS(I,4)>0.0001) THEN
                WRITE(*,*) 'ERROR: Index',g,'is not being used, so the corresponding catchability parameters'
                WRITE(*,*) '       cannot be estimated and the method indicator should be set to 0 or -0.1).'
                WRITE(12,*) 'ERROR: Index',g,'is not being used, so the corresponding catchability parameters'
                WRITE(12,*) '       cannot be estimated and the method indicator should be set to 0 or -0.1).' ; PAUSE; STOP
              ENDIF
            ENDIF
            CALL GETPAR(I,K,X,NPARM,Q_Effort(g,box,y),ETA_Q(g,box,y),SIGMA_Q(g,box,y))
            IF(EFFORT_DATA_STORE(g,box,y)>0 .OR. (PDF_EFFORT(g,box).NE.1 .AND. EFFORT_DATA_STORE(g,box,y)==0)) THEN
!              OK
             ELSE
               Q_Effort(g,box,y)=-1 ; SIGMA_Q(g,box,y)=-1
            ENDIF
          END DO
         END DO
        ELSE
         NPI=NPI+1
       ENDIF

1     CONTINUE

      IF(PDF_TAG>0) THEN ! Estimate tagging parameters

!       {fraction of animals surviving the initial tagging procedure}
        
        DO A=1,N_SURV_PARS
          IF(ICALL==0 .and. a==1 .AND. (PARM_SPECS(I,4)<0 .or. (PARM_SPECS(I,4)>0.0001 .AND. PARM_SPECS(I,4)<0.3))) GOTO 996
          CALL GETPAR(I,K,X,NPARM,SURV_PAR(a),ETA_TAG_SURV(a),SIGMA_TAG_SURV(a))
        END DO ; 
        DO J=1,N_TAG_COHORTS
          TAG_SURV(j)=SURV_PAR(TC(j)%SURV) 
        END DO

!       {chronic tag shedding or mortality}
        DO A=1,N_SHED_PARS
          IF(ICALL==0 .and. a==1 .AND. ((PARM_SPECS(I,4)<0 .and. PARM_SPECS(I,4)>-1.0) .or. &
            (PARM_SPECS(I,4)>0.0001 .AND. PARM_SPECS(I,4)<0.3))) GOTO 997
          CALL GETPAR(I,K,X,NPARM,LOSS_PAR(a),ETA_TAG_LOSS(a),SIGMA_TAG_LOSS(a))
        END DO ;
        DO J=1,N_TAG_COHORTS
          TAG_LOSS(j)=LOSS_PAR(TC(j)%SHED)
        END DO

!       {fraction of recaptured tags that are reported}
        DO A=1,N_REPT_PARS
          IF(ICALL==0 .and. a==1 .AND. ((PARM_SPECS(I,4)<0 .and. PARM_SPECS(I,4)>-1.0) .or. &
            (PARM_SPECS(I,4)>0.0001 .AND. PARM_SPECS(I,4)<0.3))) GOTO 998
          CALL GETPAR(I,K,X,NPARM,REPT_PAR(a),ETA_TAG_report(a),SIGMA_TAG_REPORT(a))
        END DO
        DO J=1,N_TAG_COHORTS
          TAG_REPORT(j,0)=REPT_PAR(TC(j)%RPT0)
          TAG_REPORT(j,1)=REPT_PAR(TC(j)%RPT1)
          TAG_REPORT(j,2)=REPT_PAR(TC(j)%RPT2)
        END DO

!       {inflation factors owing to nonmixing}
        DO A=1,N_MIX_PARS
          IF(ICALL==0. and. a==1 .AND. ((PARM_SPECS(I,4)<0 .and. PARM_SPECS(I,4)>-1.0) .or. &
            (PARM_SPECS(I,4)>0.0001 .AND. PARM_SPECS(I,4)<0.3))) GOTO 999
          CALL GETPAR(I,K,X,NPARM,NMIX_PAR(a),ETA_TAG_NOMIX(a),SIGMA_TAG_NOMIX(a))
        END DO
        DO J=1,N_TAG_COHORTS
          TAG_NOMIX(j)=NMIX_PAR(TC(j)%MIX1)
          TAG_NOMIX2(j)=NMIX_PAR(TC(j)%MIX2)
        END DO

      ENDIF

      DO A=FIRSTAGE,LASTAGE
       IF(OPTION(2)==2 .AND. MODEL_TYPE>1) THEN
         IF( ( SEL_TERMINAL(1,A,2)<0 .AND. SEL_TERMINAL(2,A,2)>=0 ) .OR.            &
             ( SEL_TERMINAL(2,A,2)<0 .AND. SEL_TERMINAL(1,A,2)>=0 ) ) THEN
           WRITE(*,'(1x,a72)')    'ERROR: You are attempting to estimate the abundance at the end of the   '
           WRITE(*,'(1x,a72)')    '       last year using the overlap recursion, which solves for both     '
           WRITE(*,'(1x,a72)')    '       stock simultaneously and therefore requires the same age classes '
           WRITE(*,'(1x,a72)')    '       to be estimated for both stocks.                                 '
           WRITE(*,'(1x,a56,i2)') '       Here you tried to estimated the abundance of age ',a
           WRITE(*,'(1x,a72)')    '       at the end of the last year for only one of the two stocks.      '
           WRITE(12,'(1x,a72)')    'ERROR: You are attempting to estimate the abundance at the end of the   '
           WRITE(12,'(1x,a72)')    '       last year using the overlap recursion, which solves for both     '
           WRITE(12,'(1x,a72)')    '       stock simultaneously and therefore requires the same age classes '
           WRITE(12,'(1x,a72)')    '       to be estimated for both stocks.                                 '
           WRITE(12,'(1x,a56,i2)') '       Here you tried to estimated the abundance of age ',a
           WRITE(12,'(1x,a72)')    '       at the end of the last year for only one of the two stocks.      '
           PAUSE ; STOP
         ENDIF
       ENDIF
      END DO
      IF(FIRSTCALL==0 .AND. I<PARM_KEY(NPAR+1)+1) THEN
          WRITE(*,'(1x,a46,a50)') 'ERROR: Too many parameters in parameter file: ',infile(2)
          WRITE(12,'(1x,a46,a50)') 'ERROR: Too many parameters in parameter file: ',infile(2) ; PAUSE; STOP
        ELSEIF(FIRSTCALL==0 .AND. I>PARM_KEY(NPAR+1)+1) THEN
          WRITE(*,'(1x,a45,a50)') 'ERROR: Too few parameters in parameter file: ',infile(2)
          WRITE(12,'(1x,a45,a50)') 'ERROR: Too few parameters in parameter file: ',infile(2) ; PAUSE; STOP
      ENDIF
      RETURN
991   WRITE(*,*) 'ERROR: The terminal year parameters do not accomodate random walk or random deviation models'
      WRITE(12,*) 'ERROR: The terminal year parameters do not accomodate random walk or random deviation models'; PAUSE ; STOP
992   WRITE(*,*) 'ERROR: You are attempting to link an F-ratio parameter to a terminal year parameter',I,PARM_SPECS(I,4)
      WRITE(12,*) 'ERROR: You are attempting to link an F-ratio parameter to a terminal year parameter',I,PARM_SPECS(I,4); PAUSE ; STOP
993   WRITE(*,*) 'ERROR: You are attempting to link a natural mortality parameter to an F-ratio'
      WRITE(12,*) 'ERROR: You are attempting to link a natural mortality parameter to an F-ratio'; PAUSE ;STOP
994   WRITE(*,*) 'ERROR: You are attempting to link a transfer coefficient to a natural mortality rate'
      WRITE(12,*) 'ERROR: You are attempting to link a transfer coefficient to a natural mortality rate'; PAUSE ;STOP
995   WRITE(*,*) 'ERROR: You are attempting to link a catchability coefficient to a variance-scaling parameter'
      WRITE(12,*) 'ERROR: You are attempting to link a catchability coefficient to a variance-scaling parameter'; PAUSE ;STOP
996   WRITE(*,*) 'ERROR: You are attempting to link tag survival to a stock-recruit (or catchability) parameter'
      WRITE(12,*) 'ERROR: You are attempting to link tag survival to a stock-recruit (or catchability) parameter'; PAUSE ;STOP
997   WRITE(*,*) 'ERROR: You are attempting to link a tag shedding parameter to a tag survival fraction'
      WRITE(12,*) 'ERROR: You are attempting to link a tag shedding parameter to a tag survival fraction'; PAUSE ; STOP
998   WRITE(*,*) 'ERROR: You are attempting to link a tag reporting parameter to a tag shedding parameter'
      WRITE(12,*) 'ERROR: You are attempting to link a tag reporting parameter to a tag shedding parameter'; PAUSE ; STOP
999   WRITE(*,*) 'ERROR: You are attempting to link a tag nonmixing parameter to a tag reporting parameter'
      WRITE(12,*) 'ERROR: You are attempting to link a tag nonmixing parameter to a tag reporting parameter'; PAUSE ; STOP
1000  WRITE(*,*) 'ERROR: You are attempting to link a variance parameter to a stock-recruit parameter'
      WRITE(12,*) 'ERROR: You are attempting to link a variance parameter to a stock-recruit parameter'; PAUSE ; STOP
      END

!------------------------------------------------------------------------------
      SUBROUTINE WRITE_PARAMETERS
!------------------------------------------------------------------------------
      USE STATISTICS ; USE PARAMETERS, ONLY : INFILE,NBOX,OPTION,PDF_TAG
      IMPLICIT NONE
      INTEGER I,K,A,Y,g
      K=0 ; I=1
      OPEN(50,FILE=INFILE(4),STATUS='UNKNOWN')
      WRITE(*,*) 'WRITING TO FILE ',INFILE(4)

      DO BOX=1,NBOX

        WRITE(50,*) '# terminal age structure of population'
        DO A=FIRSTAGE,LASTAGE-1 ; CALL WRITEPAR(I,K) ; END DO

        WRITE(50,*) '# terminal-age f-ratios'
        DO Y=FIRSTYEAR,LASTYEAR ; CALL WRITEPAR(I,K) ; END DO

        WRITE(50,*) '# natural mortality'
        DO A=FIRSTAGE,LASTAGE ; CALL WRITEPAR(I,K) ; END DO

        WRITE(50,*) '# transfer coefficients'
        DO A=FIRSTAGE,LASTAGE ; CALL WRITEPAR(I,K) ; END DO

        WRITE(50,*) '# Stock recruitment relationship'
        DO A=1,5 ; CALL WRITEPAR(I,K) ; END DO

        WRITE(50,*) '# Variance scaling parameters'
        DO G=1,NGEARS(box) ; CALL WRITEPAR(I,K) ; END DO

        IF(OPTION(4)>0) THEN
          WRITE(50,*) '# Catchability coefficients'
            DO G=1,NGEARS(box)
              DO y=FIRSTYEAR,LASTYEAR ; CALL WRITEPAR(I,K) ; END DO
            END DO
        ENDIF

      END DO ! BOX

      IF(PDF_TAG>0) THEN

        WRITE(50,*) '# fraction surviving the intial tagging process'
        DO A=1,N_SURV_PARS ; CALL WRITEPAR(I,K) ; END DO

        WRITE(50,*) '# chronic tag loss (shedding) rate'
        DO A=1,N_SHED_PARS ; CALL WRITEPAR(I,K) ; END DO

        WRITE(50,*) '# tag reporting rate'
        DO Y=1,N_REPT_PARS ; CALL WRITEPAR(I,K) ; END DO

        WRITE(50,*) '# tag nonmixing factor for first year'
        DO Y=1,N_MIX_PARS ;  CALL WRITEPAR(I,K) ; END DO

      ENDIF

      WRITE(50,*) '# Number of parameters = ',i-1
      CLOSE(50) ; RETURN ; END

!------------------------------------------------------------------------------
      SUBROUTINE REDEFINE_PARAMETERS
! redefine parameter for retrospective analyses to reflect shortened time series
!------------------------------------------------------------------------------
      USE STATISTICS ; USE PARAMETERS; USE DATUM, ONLY: SUMRDATA,RECAPTURE_DATA,RECOVERY_DATA
      IMPLICIT NONE
      INTEGER :: I,K,A,Y,g,DROP,SUMDROP,IOLD,BREC,YREC,J
      K=0 ; I=0 ; SUMDROP=0 ; IOLD=0

      DO BOX=1,NBOX

! terminal age structure of population'
        DROP=0
        DO A=FIRSTAGE,LASTAGE-1 ; CALL COUNTPAR(I,K,IOLD,DROP,SUMDROP) ; END DO

! terminal-age f-ratios'
        DO Y=FIRSTYEAR,LASTYEAR-1 ; CALL COUNTPAR(I,K,IOLD,DROP,SUMDROP) ; END DO
        DROP=1 ; CALL COUNTPAR(I,K,IOLD,DROP,SUMDROP)

! natural mortality'
        DROP=0
        DO A=FIRSTAGE,LASTAGE ; CALL COUNTPAR(I,K,IOLD,DROP,SUMDROP) ; END DO

! transfer coefficients'
        DO A=FIRSTAGE,LASTAGE ; CALL COUNTPAR(I,K,IOLD,DROP,SUMDROP) ; END DO

! Stock recruitment relationship'
        DO A=1,5 ; CALL COUNTPAR(I,K,IOLD,DROP,SUMDROP) ; END DO
      
! Variance scaling parameters'
        DO G=1,NGEARS(box) ; CALL COUNTPAR(I,K,IOLD,DROP,SUMDROP) ; END DO

        IF(OPTION(4)>0) THEN
! Catchability coefficients'
           DO G=1,NGEARS(box)
             DROP=0
             DO y=FIRSTYEAR,LASTYEAR-1 ; CALL COUNTPAR(I,K,IOLD,DROP,SUMDROP) ; END DO
             DROP=1 ; CALL COUNTPAR(I,K,IOLD,DROP,SUMDROP)
           END DO
        ENDIF

      END DO ! BOX

      IF(PDF_TAG>0) THEN

! fractionS surviving the intial tagging process'
        DROP=0
        DO A=1,N_SURV_PARS ; CALL COUNTPAR(I,K,IOLD,DROP,SUMDROP) ; END DO

! tag shedding rates
        DO A=1,N_SHED_PARS ; CALL COUNTPAR(I,K,IOLD,DROP,SUMDROP) ; END DO

! tag reporting rates
        DO A=1,N_REPT_PARS ; CALL COUNTPAR(I,K,IOLD,DROP,SUMDROP) ; END DO

! tag nonmixing factors
        DO A=1,N_MIX_PARS ; CALL COUNTPAR(I,K,IOLD,DROP,SUMDROP) ; END DO

! tag recovery data
        SUMRDATA=0
        DO J=1,N_TAG_COHORTS
          DO BREC=1,NBOX
            DO YREC=TC(j)%Yr,TC(j)%Ye
              IF(YREC<=LASTYEAR .AND. RECAPTURE_DATA(j,brec,yrec)>0) SUMRDATA(j)=SUMRDATA(j)+RECAPTURE_DATA(j,brec,yrec);
            END DO
            IF(TC(j)%Ye>LASTYEAR .AND. RECOVERY_DATA(j,brec)>0)      SUMRDATA(j)=SUMRDATA(j)+RECOVERY_DATA(j,brec);
          END DO
          IF(TC(j)%Ye>LASTYEAR) THEN ; TC(j)%Ye=LASTYEAR; TC(j)%Me=1.0; ENDIF
        END DO
        CLOSE(10) ; RETURN

      ENDIF

      PARM_KEY(0)=K     ! revised total number of parameters to be estimated
      PARM_KEY(NPAR+1)=I ! revised total number of parameters whether estimated or not
      RETURN ; END

!----------------------------------------------------------------------------        
      SUBROUTINE COUNTPAR(I,K,IOLD,DROP,SUMDROP)
! Cuts out parameters associated with dropped years in retrospective analysis
!----------------------------------------------------------------------------
      USE PARAMETERS, ONLY : PARM_SPECS, PARM_KEY, PARM_EST
      IMPLICIT NONE
      INTEGER :: I,K,DROP,SUMDROP,IOLD,J
      IOLD=IOLD+1 ; SUMDROP=SUMDROP+DROP
      IF(DROP==0) THEN
        I=I+1
        DO J=1,5 ; PARM_SPECS(I,J)=PARM_SPECS(IOLD,J) ; END DO
        IF(PARM_SPECS(I,4)>0) THEN
          K=K+1 ; PARM_KEY(K)=I ; PARM_EST(K)=PARM_SPECS(IOLD,2)
        ENDIF
      ENDIF
      RETURN ; END


!-----------------------------------------------------------------------------
      SUBROUTINE MODEL
! Performs a VPA (backwards, exact solution) to start of terminal year
!-----------------------------------------------------------------------------
      USE STATISTICS ; USE DATUM ; USE PARAMETERS, ONLY : MODEL_TYPE, OPTION, PDF_TAG, SEED, NBOX
      IMPLICIT NONE
      INTEGER Y,G,A,S,LAST,YOUNGER,OTH,NYR
      REAL (KIND=8) :: Z(Bs,0:As,0:Ys),SURV(Bs,0:As,0:Ys),NTEMP(Bs,0:As,0:Ys),NTIME,TEMP(2),SUMTEMP(2), &
                       F2(2),fvec2(2),df2(2,2),RAN1,STD,AVE,MAX_EFF,SUM_NTIME
      LOGICAL CHECK
      EXTERNAL DIFFUSION_CATCH, DIFFUSION_PLUS, OVERLAP_CATCH, OVERLAP_PLUS, RAN1
      ssb=0.d0 
      
      IF(PLUSGROUP>1) THEN ; LAST=2 ; ELSE ; LAST=1 ; ENDIF

!!!!!!!!!!!!!!!!!!!!
! get F's and N's  !
!!!!!!!!!!!!!!!!!!!!
!     NOTE: A maximum value of 5.0 is allowed because it seems unlikely in most applications
!           that fishing could reduce the population by over 99% in a single year (1-exp{-5.0}=99.33)
!           The minimum value (0.00000001) is an arbitrary small value (since F>0 in VPA).
!           Without these constraints the VPA sometimes will not arrive at a solution.

      IF(MODEL_TYPE==1) THEN ! Diffusion (no-memory) and single-stock model
      
        YEAR=LASTYEAR
        IF(OPTION(2)==2) THEN
!         Get terminal-year F's for the age groups with estimates of N
          DO AGE=FIRSTAGE,LASTAGE-1
            DO BOX=1,NBOX
              IF(SEL_TERMINAL(box,age,2)<0) THEN
                CALL DIFFUSION_CATCH(1,F2,FVEC2,DF2,-1) ! obtain initial quess for F
                CALL NEWTON(F2,1,CHECK,DIFFUSION_CATCH) ! Use Newton method to iterate final solution
!                IF(.NOT. CHECK .OR. F2(1)<0.00001 .OR. F2(1)>5.0) THEN ! A good solution was not obtained, try quasi-Newton method
!                  CALL DIFFUSION_CATCH(1,F2,FVEC2,DF2,-1) ; CALL NEWT(F2,1,CHECK,DIFFUSION_CATCH) ; ENDIF
                F(box,age,year)=MIN(5.0,MAX(F2(1),0.000000001))
              ENDIF
          END DO ; END DO
        ENDIF

!       Get terminal-year F's for the age groups with assumed relative selectivities or estimated F-ratios
        DO 2 A=FIRSTAGE,LASTAGE
          DO 3 BOX=1,NBOX
            IF(SEL_TERMINAL(box,a,2)>0 .OR. A==LASTAGE) THEN
                IF(A==LASTAGE .AND. PLUSGROUP>1) THEN ! Use the F-ratio
                    F(box,a,year)=F(box,a-1,year)*F_RATIO(box,year)
                  ELSE ! Use the relative selectivities
                    S=INT(SEL_TERMINAL(box,a,2)) ; F(box,a,year)=F(box,s,year)*SEL_TERMINAL(box,a,1)
                ENDIF
            ENDIF
            Z(box,a,year)=F(box,a,year)+M(box,a,year) ; IF(Z(box,a,year)>20) Z(box,a,year)=20 ; SURV(box,a,year)=EXP(-Z(box,a,year))
            Ntemp(box,a,year)=CATCH_DATA(box,a,year)*Z(box,a,year)/F(box,a,year)/(1-SURV(box,a,year))
3         CONTINUE
!         get N's prior to mixing
          CALL CRAMERS_RULE(N_AREA(1,a,year),N_AREA(2,a,year),1-T(1,a,year),T(2,a,year),T(1,a,year),1-T(2,a,year), &
                            Ntemp(1,a,year),Ntemp(2,a,year))
!          write(*,*) a,year,ntemp(1,a,year),ntemp(1,a,year),N_AREA(1,a,year),N_AREA(2,a,year) ; pause
2       CONTINUE

!       Get F's from estimates of N in remaining years
        DO 10 Y=LASTYEAR-1,FIRSTYEAR, -1
          YEAR=Y
!         VPA for LAST ages
          IF(PLUSGROUP>1) THEN
              YOUNGER=PLUSGROUP-1
!             get F's and N's subsequent to mixing
              DO 13 BOX=1,NBOX
                CALL DIFFUSION_PLUS(1,F2,FVEC2,DF2,-1) ! obtain initial quess for F
!                write(*,*) '2',f2(1)
                CALL NEWTON(F2,1,CHECK,DIFFUSION_PLUS) ! Use Newton method to iterate final solution
!                IF(.NOT. CHECK .OR. F2(1)<0.00001 .OR. F2(1)>5.0) THEN ! A good solution not obtained, try quasi-Newton method
!                  CALL DIFFUSION_PLUS(1,F2,FVEC2,DF2,-1) ; CALL NEWT(F2,1,CHECK,DIFFUSION_PLUS) ; ENDIF
                F(box,younger,y)=MIN(5.0,MAX(F2(1),0.000000001))
                F(box,plusgroup,y)=F_RATIO(box,year)*F(box,younger,y)
                Z(box,younger,y)=F(box,younger,y)+M(box,younger,y) ; IF(Z(box,younger,y)>20) Z(box,younger,y)=20
                SURV(box,younger,y)=EXP(-Z(box,younger,y))
                Z(box,plusgroup,y)=F(box,plusgroup,y)+M(box,plusgroup,y) ; IF(Z(box,plusgroup,y)>20) Z(box,plusgroup,y)=20
                SURV(box,plusgroup,y)=EXP(-Z(box,plusgroup,y))
                Ntemp(box,younger,y)=CATCH_DATA(box,younger,y)*Z(box,younger,y)/F(box,younger,y)/(1-SURV(box,younger,y))
                Ntemp(box,plusgroup,y)=CATCH_DATA(box,plusgroup,y)*Z(box,plusgroup,y)/F(box,plusgroup,y)/(1-SURV(box,plusgroup,y))
13            CONTINUE
!             get N's prior to mixing
              CALL CRAMERS_RULE(N_AREA(1,younger,y),N_AREA(2,younger,y),1-T(1,younger,y),T(2,younger,y),T(1,younger,y), &
                                1-T(2,younger,y),Ntemp(1,younger,y),Ntemp(2,younger,y))
              CALL CRAMERS_RULE(N_AREA(1,plusgroup,y),N_AREA(2,plusgroup,y),1-T(1,plusgroup,y),T(2,plusgroup,y),T(1,plusgroup,y),&
                                1-T(2,plusgroup,y),Ntemp(1,plusgroup,y),Ntemp(2,plusgroup,y))
            ELSE
              WRITE(*,*) 'THERE ARE NOT YET ANY ALTERNATIVES TO HAVING A PLUS GROUP'
              PAUSE; STOP
          ENDIF

!         VPA for ages other than LAST oldest
          DO 11 age = FIRSTAGE,LASTAGE-LAST
!           get F's and N's subsequent to mixing
            DO 12 BOX=1,NBOX
              CALL DIFFUSION_CATCH(1,F2,FVEC2,DF2,-1) ! obtain initial quess for F
              CALL NEWTON(F2,1,CHECK,DIFFUSION_CATCH) ! Use Newton method to iterate final solution
!              IF(.NOT. CHECK .OR. F2(1)<0.00001 .OR. F2(1)>5.0) THEN ! A good solution was not obtained, try quasi-Newton method
!                CALL DIFFUSION_CATCH(1,F2,FVEC2,DF2,-1) ; CALL NEWT(F2,1,CHECK,DIFFUSION_CATCH) ; ENDIF
              F(box,age,y)=MIN(5.0,MAX(F2(1),0.000000001))
              Z(box,age,y)=F(box,age,y)+M(box,age,y) ; SURV(box,age,y)=EXP(-Z(box,age,y))
              Ntemp(box,age,y)=N_AREA(box,age+1,y+1)/SURV(box,age,y)
!             write(*,*) '3',age,f2(1)
12          CONTINUE
!           get N's prior to mixing
            CALL CRAMERS_RULE(N_AREA(1,age,y),N_AREA(2,age,y),1-T(1,age,y),T(2,age,y),T(1,age,y),1-T(2,age,y), &
                              Ntemp(1,age,y),Ntemp(2,age,y))
11        CONTINUE

10      CONTINUE

      ELSE  ! Overlap model  

!       Get terminal-year F's for the age groups with estimates of N
        YEAR=LASTYEAR
        IF(OPTION(2)==2) THEN
          DO AGE=FIRSTAGE,LASTAGE-1
            IF(SEL_TERMINAL(1,age,2)<0) THEN
                CALL OVERLAP_CATCH(2,F2,FVEC2,DF2,-1) ! obtain initial quess for F
                CALL NEWTON(F2,2,CHECK,OVERLAP_CATCH) ! Use Newton method to iterate final solution
                F(1,age,year)=MIN(5.0,MAX(F2(1),0.00000001)) ; F(2,age,year)=MIN(5.0,MAX(F2(2),0.00000001))
            ENDIF
          END DO
        ENDIF

!       Get terminal-year F's for the age groups with relative selectivities
        DO 102 AGE=FIRSTAGE,LASTAGE
          DO 101 BOX=1,2
            IF(SEL_TERMINAL(box,age,2)>0 .OR. Age==LASTAGE) THEN
              IF(Age==LASTAGE .AND. PLUSGROUP>1) THEN ! Use the F-ratio
                  F(box,age,year)=F(box,age-1,year)*F_RATIO(box,year)
                ELSE ! Use the relative selectivities
                  S=INT(SEL_TERMINAL(box,age,2))
!                  write(*,*) SEL_TERMINAL(box,age,2) ; write(*,*) box,age,year,s,plusgroup
                  F(box,age,year)=F(box,s,year)*SEL_TERMINAL(box,age,1)
              ENDIF
            ENDIF
101       CONTINUE
          Z(1,age,year)=F(1,age,year)+M(1,age,year) ; SURV(1,age,year)=EXP(-Z(1,age,year))
          Z(2,age,year)=F(2,age,year)+M(2,age,year) ; SURV(2,age,year)=EXP(-Z(2,age,year))
          TEMP(1)=CATCH_DATA(1,age,year)*Z(1,age,year)/F(1,age,year)/(1-SURV(1,age,year))
          TEMP(2)=CATCH_DATA(2,age,year)*Z(2,age,year)/F(2,age,year)/(1-SURV(2,age,year))
          CALL CRAMERS_RULE(N_STOCK(1,age,year),N_STOCK(2,age,year),1-T(1,age,year),T(2,age,year),T(1,age,year),1-T(2,age,year),  &
                            temp(1),temp(2))
102     CONTINUE

!       Get F's from estimates of N in remaining years
        DO 110 Y=LASTYEAR-1,FIRSTYEAR, -1
          YEAR=Y
!         VPA for LAST ages
          IF(PLUSGROUP>1) THEN
              YOUNGER=PLUSGROUP-1
!             get F's and N's subsequent to mixing
              CALL OVERLAP_PLUS(2,F2,FVEC2,DF2,-1) ! obtain initial quess for F
              CALL NEWTON(F2,2,CHECK,OVERLAP_PLUS) ! Use Newton method to iterate final solution
!              IF(.NOT. CHECK .OR. F2(1)<0.0001 .OR. F2(1)>5.0 .OR. F2(2)<0.0001 .OR. F2(2)>5.0) THEN ! Try quasi-Newton method
!                CALL OVERLAP_PLUS(2,F2,FVEC2,DF2,-1) ; CALL NEWT(F2,2,CHECK,OVERLAP_PLUS) ; ENDIF
                F(1,younger,y)=MIN(5.0,MAX(F2(1),0.00000001)) ; F(2,younger,y)=MIN(5.0,MAX(F2(2),0.00000001))
                F(1,plusgroup,y)=F_RATIO(1,year)*F(1,younger,y)
                F(2,plusgroup,y)=F_RATIO(2,year)*F(2,younger,y)
                Z(1,younger,y)=F(1,younger,y)+M(1,younger,y) ; SURV(1,younger,y)=EXP(-Z(1,younger,y))
                Z(2,younger,y)=F(2,younger,y)+M(2,younger,y) ; SURV(2,younger,y)=EXP(-Z(2,younger,y))
                Z(1,plusgroup,y)=F(1,plusgroup,y)+M(1,plusgroup,y) ; SURV(1,plusgroup,y)=EXP(-Z(1,plusgroup,y))
                Z(2,plusgroup,y)=F(2,plusgroup,y)+M(2,plusgroup,y) ; SURV(2,plusgroup,y)=EXP(-Z(2,plusgroup,y))
                TEMP(1)=CATCH_DATA(1,younger,year)*Z(1,younger,year)/F(1,younger,year)/(1-SURV(1,younger,year))
                TEMP(2)=CATCH_DATA(2,younger,year)*Z(2,younger,year)/F(2,younger,year)/(1-SURV(2,younger,year))
                CALL CRAMERS_RULE(N_STOCK(1,younger,y),N_STOCK(2,younger,y),1-T(1,younger,y),T(2,younger,y),T(1,younger,y), &
                                  1-T(2,younger,y),temp(1),temp(2))
                TEMP(1)=CATCH_DATA(1,plusgroup,year)*Z(1,plusgroup,year)/F(1,plusgroup,year)/(1-SURV(1,plusgroup,year))
                TEMP(2)=CATCH_DATA(2,plusgroup,year)*Z(2,plusgroup,year)/F(2,plusgroup,year)/(1-SURV(2,plusgroup,year))
                CALL CRAMERS_RULE(N_STOCK(1,plusgroup,y),N_STOCK(2,plusgroup,y),1-T(1,plusgroup,y),T(2,plusgroup,y),  &
                                  T(1,plusgroup,y),1-T(2,plusgroup,y),temp(1),temp(2))
113           CONTINUE
            ELSE
              WRITE(*,*) 'THERE ARE NOT YET ANY ALTERNATIVES TO HAVING A PLUS GROUP'
              PAUSE; STOP
          ENDIF

!         VPA for ages other than LAST oldest
          DO 111 age = FIRSTAGE,LASTAGE-LAST
!           get F's and N's subsequent to mixing
            CALL OVERLAP_CATCH(2,F2,FVEC2,DF2,-1) ! obtain initial quess for F
            CALL NEWTON(F2,2,CHECK,OVERLAP_CATCH) ! Use Newton method to iterate final solution
!            IF(.NOT. CHECK .OR. F2(1)<0.0001 .OR. F2(1)>5.0 .OR. F2(2)<0.0001 .OR. F2(2)>5.0) THEN ! Try the quasi-Newton method
!              CALL OVERLAP_CATCH(2,F2,FVEC2,DF2,-1) ; CALL NEWT(F2,2,CHECK,OVERLAP_CATCH) ; ENDIF
            F(1,age,year)=MIN(5.0,MAX(F2(1),0.000001)) ; F(2,age,year)=MIN(5.0,MAX(F2(2),0.000001))
!            F(1,age,year)=F2(1) ; F(2,age,year)=F2(2)
            Z(1,age,year)=F(1,age,year)+M(1,age,year) ; SURV(1,age,year)=EXP(-Z(1,age,year))
            Z(2,age,year)=F(2,age,year)+M(2,age,year) ; SURV(2,age,year)=EXP(-Z(2,age,year))
            TEMP(1)=CATCH_DATA(1,age,year)*Z(1,age,year)/F(1,age,year)/(1-SURV(1,age,year))
            TEMP(2)=CATCH_DATA(2,age,year)*Z(2,age,year)/F(2,age,year)/(1-SURV(2,age,year))
            CALL CRAMERS_RULE(N_STOCK(1,age,y),N_STOCK(2,age,y),1-T(1,age,y),T(2,age,y),T(1,age,y),1-T(2,age,y),temp(1),temp(2))
111       CONTINUE

110     CONTINUE

      ENDIF ! {model type} 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Expand plus-group                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(EXPANDED_PLUSGROUP>PLUSGROUP) THEN

!       equilibrium approximation in first year
        y=FIRSTYEAR
        DO BOX=1,NBOX
          sumtemp(box)=1; Ntemp(box,plusgroup,y)=1;
          DO a = PLUSGROUP+1,EXPANDED_PLUSGROUP
            age=a-1; IF(a>PLUSGROUP) age=PLUSGROUP
            Ntemp(box,a,y)=Ntemp(box,a-1,y)*SURV(box,age,y);
            IF(a==EXPANDED_PLUSGROUP) Ntemp(box,a,y)=Ntemp(box,a,y)/(1.0-SURV(box,plusgroup,y))
            sumtemp(box)=sumtemp(box)+Ntemp(box,a,y)
          END DO
          IF(MODEL_TYPE==1) THEN ! Diffusion (no-memory) and single-stock model
              temp(box)=N_AREA(box,plusgroup,y)
              DO a = PLUSGROUP,EXPANDED_PLUSGROUP; N_AREA(box,a,y)=temp(box)*Ntemp(box,a,y)/sumtemp(box); END DO
            ELSE
              temp(box)=N_STOCK(box,plusgroup,y)
              DO a = PLUSGROUP,EXPANDED_PLUSGROUP; N_STOCK(box,a,y)=temp(box)*Ntemp(box,a,y)/sumtemp(box); END DO
          ENDIF
        END DO

!       track cohorts in subsequent years
        DO Y=FIRSTYEAR,LASTYEAR
          DO BOX=1,NBOX
            IF(BOX==1) THEN ; OTH=2 ; ELSE ; OTH=1 ; ENDIF
            DO a = FIRSTAGE,EXPANDED_PLUSGROUP
              age=a; IF(age>PLUSGROUP) age=PLUSGROUP
              IF(MODEL_TYPE==1) THEN ! diffusion and no-mixing models
                  N_AREA(box,a+1,y+1)=( N_AREA(box,a,y)*(1-T(box,age,y))+N_AREA(oth,a,y)*T(oth,age,y) )*SURV(box,age,y)
                  IF(A==EXPANDED_PLUSGROUP) N_AREA(box,a,y+1)=N_AREA(box,a,y+1)+N_AREA(box,a+1,y+1)
                ELSE ! overlap model
                  N_STOCK(box,a+1,y+1)=N_STOCK(box,a,y)*( (1-T(box,age,y))*SURV(box,age,y) + T(box,age,y)*SURV(oth,age,y) )
                  IF(A==EXPANDED_PLUSGROUP) N_STOCK(box,a,y+1)=N_STOCK(box,a,y+1)+N_STOCK(box,a+1,y+1)
              ENDIF
            END DO
          END DO
        END DO

      ENDIF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute abundance indices and spawning stock biomass !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO Y=FIRSTYEAR,LASTYEAR
        DO a = FIRSTAGE,EXPANDED_PLUSGROUP
          age=a; IF(age>PLUSGROUP) age=PLUSGROUP
          DO BOX=1,NBOX
            IF(BOX==1) THEN ; OTH=2 ; ELSE ; OTH=1 ; ENDIF
            IF(SEASON_SSB(box)<0) THEN ; NTIME=(1-SURV(box,age,y))/Z(box,age,y)
              ELSE ; NTIME=SURV(box,age,y)**(SEASON_SSB(box)/12.d0) ; ENDIF
!           recall that Ntemp is the total number in each area after mixing, Nstock is the number in each stock regardless of area.
            IF(MODEL_TYPE==2) THEN ! overlap
                NTEMP(box,a,y)=N_STOCK(box,a,y)*(1-T(box,age,y))+N_STOCK(oth,a,y)*T(oth,age,y)
                ssb(box,y)=ssb(box,y)+FECUNDITY(box,a)*N_STOCK(box,a,y)*WEIGHT_SSB(box,age,y)*NTIME/1000
              ELSE ! diffusion and single-stock stock
                Ntemp(box,a,y)=N_AREA(box,a,y)*(1-T(box,age,y))+N_AREA(oth,a,y)*T(oth,age,y)
                ssb(box,y)=ssb(box,y)+FECUNDITY(box,a)*Ntemp(box,a,y)*WEIGHT_SSB(box,age,y)*NTIME/1000
            ENDIF
      END DO ; END DO ; END DO

!     indices of abundance
      EFFORT=0.D0  
      DO 30 BOX=1,NBOX
        IF(BOX==1) THEN ; OTH=2 ; ELSE ; OTH=1 ; ENDIF
        DO 31 G=1,NGEARS(box)
          IF(PDF_EFFORT(g,box)==0) GOTO 31 ; GEAR=g ; MAX_EFF=-1
          DO 32 Y = FIRSTYEAR,LASTYEAR
            IF(EFFORT_DATA_STORE(g,box,y)>0 .OR. (PDF_EFFORT(g,box)>1 .AND. EFFORT_DATA_STORE(g,box,y)>-9)) THEN
              IF(SEL_TYPE(g,box)>1) CALL GET_SELECTIVITY ! compute from fractional catches
              IF(AGE_EFFORT(1,g,box)<=PLUSGROUP .AND. AGE_EFFORT(2,g,box)>=PLUSGROUP  ) THEN ! Set selectivity values for expanded ages equal to that of the plus-group
                 DO A=PLUSGROUP+1,AGE_EFFORT(2,gear,box); SEL_EFFORT(gear,box,a,y)=SEL_EFFORT(gear,box,plusgroup,y); END DO
              ENDIF
              SUM_NTIME=0
              DO 33 A=AGE_EFFORT(1,g,box),AGE_EFFORT(2,g,box)
                age=a; IF(age>PLUSGROUP) age=PLUSGROUP
                IF(SEASON_EFFORT(g,box)<0) THEN ; NTIME=(1-SURV(box,age,y))/Z(box,age,y)
                  ELSE ; NTIME=SURV(box,age,y)**(SEASON_EFFORT(g,box)/12.d0)
                ENDIF
                IF(BIO_EFFORT(g,box)==2) THEN  ! index is biomass-based
                    NTIME=NTIME*WEIGHT_CATCH(g,box,age,y) ! note: for overlap model one could modify to allowed weights by stock
                  ELSEIF(BIO_EFFORT(g,box)==3 .OR. BIO_EFFORT(g,box)==5) THEN ! index is in terms of fishing mortality rate
                    SUM_NTIME=SUM_NTIME+NTIME*NTEMP(box,a,y)*SEL_EFFORT(g,box,a,y) ; NTIME=NTIME*F(box,age,y)
                  ELSEIF(BIO_EFFORT(g,box)==4 .OR. BIO_EFFORT(g,box)==6) THEN ! index is in terms of total mortality rate
                    SUM_NTIME=SUM_NTIME+NTIME*NTEMP(box,a,y)*SEL_EFFORT(g,box,a,y) ; NTIME=NTIME*Z(box,age,y)
                  ELSEIF(BIO_EFFORT(g,box)==7) THEN ! index is the proportion of number in area that originated in box
                    SUM_NTIME=SUM_NTIME+SEL_EFFORT(g,box,a,y)*NTIME*N_STOCK(box,a,y)*(1-T(box,age,y))
                  ELSEIF(BIO_EFFORT(g,box)==8) THEN ! index is the proportion of biomass in area that originated in box
                    NTIME=NTIME*WEIGHT_CATCH(g,box,age,y)
                    SUM_NTIME=SUM_NTIME+SEL_EFFORT(g,box,a,y)*NTIME*N_STOCK(box,a,y)*(1-T(box,age,y))
                ENDIF
                EFFORT(g,box,y)=EFFORT(g,box,y)+SEL_EFFORT(g,box,a,y)*NTIME*NTEMP(box,a,y)
                ! write(*,*) g,y,effort(g,box,y),NTIME,NTEMP(box,a,y),SEL_EFFORT(g,box,a,y)
!               hardwire
                IF(BOTH_STOCKS(g,box)>0) THEN ! add abundance in other area as well (i.e., index represents both areas)
                  IF(SEASON_EFFORT(g,box)<0) THEN ; NTIME=(1.0-SURV(oth,age,y))/Z(oth,age,y)
                    ELSE ; NTIME=SURV(oth,age,y)**(SEASON_EFFORT(g,box)/12.d0)
                  ENDIF
                  IF(BIO_EFFORT(g,box)==2) NTIME=NTIME*WEIGHT_SSB(oth,age,y) ! this is a patch, better to specify g in WEIGHT_CATCH(g,oth... )
                  EFFORT(g,box,y)=EFFORT(g,box,y)+SEL_EFFORT(g,box,a,y)*NTIME*NTEMP(oth,a,y)
                ENDIF
33            CONTINUE
              IF(BIO_EFFORT(g,box)>2 .AND. BIO_EFFORT(g,box)<7) THEN ! F or Z indices
                EFFORT(g,box,y)=EFFORT(g,box,y)/SUM_NTIME
                IF(BIO_EFFORT(g,box)<5) MAX_EFF=MAX(MAX_EFF,EFFORT(g,box,y))
              ENDIF
            ENDIF
            IF(BIO_EFFORT(g,box)==7 .or. BIO_EFFORT(g,box)==8) THEN
              IF(EFFORT_DATA_STORE(g,box,y)>0 .AND. EFFORT_DATA(g,box,y)>0) THEN
                if(SUM_NTIME/EFFORT(g,box,y)>1) WRITE(*,*) g,box,y,EFFORT(g,box,y),SUM_NTIME,SUM_NTIME/EFFORT(g,box,y)
                EFFORT(g,box,y)=SUM_NTIME/EFFORT(g,box,y)
              ENDIF
            ENDIF

32        CONTINUE ! year
          IF(BIO_EFFORT(g,box)==3 .or. BIO_EFFORT(g,box)==4) THEN
            DO Y=FIRSTYEAR,LASTYEAR
              IF(EFFORT_DATA_STORE(g,box,y)>0 .OR. (PDF_EFFORT(g,box)>1 .AND. EFFORT_DATA_STORE(g,box,y)>-9)) THEN
                EFFORT(g,box,y)=MAX_EFF-EFFORT(g,box,y)
              ENDIF
            END DO
          ENDIF
          IF(OPTION(4)<=0 .and. BIO_EFFORT(g,box)<7 .AND. PDF_EFFORT(g,box)/=12 ) CALL GET_Q_BY_MLE ! Estimate q by mle
          DO 36 Y=FIRSTYEAR,LASTYEAR
            IF(BIO_EFFORT(g,box)>6) Q_EFFORT(g,box,y)=1.0
            IF(EFFORT_DATA_STORE(g,box,y)>0 .OR. (PDF_EFFORT(g,box)>1 .AND. EFFORT_DATA_STORE(g,box,y)>-9)) THEN
              EFFORT(g,box,y)=EFFORT(g,box,y)*Q_EFFORT(g,box,y)
              IF(PDF_EFFORT(g,box)==1.AND.BIO_EFFORT(g,box)>2.AND.BIO_EFFORT(g,box)<7) EFFORT(g,box,y)=EFFORT(g,box,y)+0.00001
            ENDIF
36        CONTINUE
          IF(PDF_EFFORT(g,box)==12) THEN
            AVE=0 ; NYR=0 ; STD=0
            DO y=firstyear,lastyear
              IF(EFFORT_DATA_STORE(g,box,y)>0 .OR. (PDF_EFFORT(g,box)>1 .AND. EFFORT_DATA_STORE(g,box,y)>-9)) THEN
                AVE=AVE+EFFORT(g,box,y) ; NYR=NYR+1
              ENDIF
            END DO
            AVE=AVE/NYR
            DO y=firstyear,lastyear
              IF(EFFORT_DATA_STORE(g,box,y)>0 .OR. (PDF_EFFORT(g,box)>1 .AND. EFFORT_DATA_STORE(g,box,y)>-9)) THEN
                STD=STD+(EFFORT(g,box,y)-AVE)**2
              ENDIF
            END DO
            STD=SQRT(STD/(NYR-1))
            DO y=firstyear,lastyear
              IF(EFFORT_DATA_STORE(g,box,y)>0 .OR. (PDF_EFFORT(g,box)>1 .AND. EFFORT_DATA_STORE(g,box,y)>-9)) THEN
                EFFORT(g,box,y)=(EFFORT(g,box,y)-AVE)/STD
              ENDIF
            END DO
          ENDIF
31      CONTINUE            
30    CONTINUE      
      
      
!     Tagging statistics 
      IF(PDF_TAG>0) CALL TAG_RECAPTURES

      RETURN ; END
 
! -----------------------------------------------------------------
      SUBROUTINE GET_Q_BY_MLE
! Computes q values by maximum likelihood estimator
! Note that the formulas differ depending on whether
! the objective function is in logarithms or not.
!------------------------------------------------------------------------------
      USE STATISTICS  ; USE DATUM, ONLY : SIGMA_EFFORT, EFFORT_DATA, PDF_EFFORT,BIO_EFFORT, EFFORT_DATA_STORE
      USE PARAMETERS, ONLY : ADD_VAR
      IMPLICIT NONE
      REAL (KIND=8) :: sumn,cross,SIG2
      INTEGER Y,G
              
      CROSS=0.d0 ; SUMN=0.d0 ; g = GEAR
      DO 2 Y=FIRSTYEAR,LASTYEAR
        IF(EFFORT_DATA_STORE(g,box,y)>=0 .AND. EFFORT(g,box,y)>0) THEN
          CALL GETVARIANCE(EFFORT_DATA_STORE(g,box,y),EFFORT(g,box,y),SIGMA_EFFORT(g,box,y),SCALE_VARIANCE(g,box), &
                           PDF_EFFORT(g,box),SIG2)
          IF(PDF_EFFORT(g,box).EQ.1) THEN
!             {cpue is lognormally distributed}
              IF(EFFORT_DATA_STORE(g,box,y).GT.0 .AND. EFFORT(g,box,y)>0) THEN
                IF(BIO_EFFORT(g,box)==3 .or. BIO_EFFORT(g,box)==4) THEN ! the year at Max_eff is excluded
                    CROSS=CROSS+DLOG(EFFORT_DATA(g,box,y)/(EFFORT(g,box,y)+0.0001))/SIG2 ; SUMN=SUMN+1.d0/SIG2
                  ELSE
                    CROSS=CROSS+DLOG(EFFORT_DATA(g,box,y)/EFFORT(g,box,y))/SIG2 ; SUMN=SUMN+1.d0/SIG2
                ENDIF
              ENDIF
            ELSEIF(PDF_EFFORT(g,box)==4) THEN
!             {cpue is Poisson distributed}
              CROSS=cross+EFFORT_DATA(g,box,y) ; SUMN=sumn+EFFORT(g,box,y)
            ELSE
!             {cpue is normally distributed}
              CROSS=cross+EFFORT(g,box,y)*EFFORT_DATA(g,box,y)/SIG2 ; SUMN=sumn+(EFFORT(g,box,y)**2)/SIG2
          ENDIF
        ENDIF
2     CONTINUE
      IF(SUMN.EQ.0.) THEN
          Q_EFFORT(g,box,firstyear)=0.D0
        ELSEIF(PDF_EFFORT(g,box)==1) THEN
          Q_EFFORT(g,box,firstyear)=dEXP(CROSS/SUMN)
        ELSEIF(PDF_EFFORT(g,box)>1) THEN
          Q_EFFORT(g,box,firstyear)=CROSS/SUMN
      ENDIF
      DO 3 Y=FIRSTYEAR+1,LASTYEAR
          Q_EFFORT(g,box,y)=Q_EFFORT(g,box,firstyear)
3     CONTINUE
      RETURN ; END
     
!-------------------------------------------------------------------
      SUBROUTINE GET_SELECTIVITY
!-------------------------------------------------------------------
      USE STATISTICS  ; USE DATUM, ONLY : SEL_TYPE, PSEL, EFFORT_DATA, PDF_EFFORT,BIO_EFFORT, PSEL_MEAN
      IMPLICIT NONE
      REAL (KIND=8) :: PMAX,FBEST(0:As)
      INTEGER Y,AGEMAX
      
      AGEMAX=AGE_EFFORT(2,gear,box); IF(AGEMAX>PLUSGROUP) AGEMAX=PLUSGROUP;
      IF(SEL_TYPE(gear,box)<=3) THEN ! Year-specific selectivities-- Powers and Restrepo (1992)
          DO Y=FIRSTYEAR,LASTYEAR
            PMAX=-1.D0
            DO AGE=AGE_EFFORT(1,gear,box),AGEMAX
              FBEST(age)=F(box,age,y)*PSEL(gear,box,age,y)
              PMAX=MAX(pmax,fbest(age))
            END DO
            DO AGE=AGE_EFFORT(1,gear,box),AGEMAX
              SEL_EFFORT(gear,box,age,y)=FBEST(age)/PMAX
            END DO
          END DO
        ELSE ! constant selectivities-- Butterworth and Geromont (1999)
          PMAX=-1.D0 ; FBEST=0
          DO AGE=AGE_EFFORT(1,gear,box),AGEMAX
            DO Y=FIRSTYEAR,LASTYEAR
              IF(EFFORT_DATA(gear,box,y).GT.0 .or. (EFFORT_DATA(gear,box,y)>-9 .and. PDF_EFFORT(gear,box)>1)) THEN
                IF(PSEL_MEAN(GEAR,BOX,Y)>0.0) FBEST(age)=FBEST(age)+F(box,age,y)*PSEL(gear,box,age,y)
              ENDIF
            END DO ! Y
            PMAX=MAX(pmax,fbest(age))
          END DO ! AGE
          DO Y=FIRSTYEAR,LASTYEAR
            IF(EFFORT_DATA(gear,box,y).GT.0 .or. (EFFORT_DATA(gear,box,y)>-9 .and. PDF_EFFORT(gear,box)>1)) THEN
              DO AGE=AGE_EFFORT(1,gear,box),AGEMAX; SEL_EFFORT(gear,box,age,y)=FBEST(age)/PMAX; END DO
            ENDIF
          END DO
      ENDIF

      RETURN ; END

!------------------------------------------------------------------------------
      SUBROUTINE TAG_RECAPTURES
! Tag attrition model
!------------------------------------------------------------------------------
      USE STATISTICS ; USE PARAMETERS ; USE DATUM, ONLY : T_PROPORTION
      IMPLICIT NONE
      INTEGER OTH,j
      REAL (KIND=8) :: BOTTOM,z1,z2,F1,S1,S2,T1,T2,T_mod,DELTA
      RECAPTURES=0; RECOVERIES=0;

! \\\\\\\\\\\\\\\\\\\\\\\\\\\              
! table of variables
! ///////////////////////////

! N_TAGS(j,k,year)    =   number of tagged fish in cohort j that belong to stock k in year
!                           (note that for the diffusion model the number in stock k is
!                           the same as the number in area k)
! RECAPTURES(j,k,year)=   number of tagged fish in cohort j that were reported as
!                           recaptured in area k
! RECOVERIES(j,k)     =   number of tagged fish in cohort j that were detected by
!                           satellite or other means at the end of the experiment
! TAG_REPORT(j,0)     =   reporting rate to satellite or other means of detection
! TAG_REPORT(j,k)     =   reporting rate by fisheries operating in area k
! T_proportion(i)     =   chance of a tagged animal moving after month i relative to
!                           the chance that it would have moved had it been tagged at
!                           the beginning of the year

! \\\\\\\\\\\\\\\\\\\\\\\\\\\
! calculations
! ///////////////////////////

      DO J=1,N_TAG_COHORTS

        IF(TC(j)%N<=0) CYCLE

        BOX=TC(j)%Br; IF(BOX==1) THEN ; OTH=2 ; ELSE ; OTH=1 ; ENDIF

        DO YEAR=TC(j)%Yr,TC(j)%Ye

          AGE=TC(j)%Ar+Year-TC(j)%Yr; IF(AGE>PLUSGROUP) AGE=PLUSGROUP
          IF( TC(j)%Yr==TC(j)%Ye  ) THEN;  DELTA=TC(j)%Me-TC(j)%Mr
            ELSE IF(YEAR==TC(j)%Yr) THEN; DELTA=1.0 - TC(j)%Mr
            ELSE IF(YEAR< TC(j)%Ye) THEN; DELTA=1.0
            ELSE                        ; DELTA=TC(j)%Me
          ENDIF
          T1=T(box,age,year) ; T2=T(oth,age,year)
          IF(YEAR==TC(j)%Yr) THEN ;
!             {accomodates possibility that tags do not randomly mix throughout the zone in their first calendar year at large
!              and so may have experience a different fishing mortality rate than the general population, but
!              assumes tags that cross to the other zone do mix throughout the zone}
              F1=TAG_NOMIX(j)*F(box,age,year);
!             {accomodate possibility that tags released in certain months during the first calendar year may be less likely to move}
              T1=T(box,age,year)*T_proportion(T_mod_month(j));
            ELSE IF(YEAR==TC(j)%Yr+1) THEN ;
!             {accomodate possibility that tags do not randomly mix throughout the zone in their second calendar year at large}
              F1=TAG_NOMIX2(j)*F(box,age,year);
            ELSE ;
              F1=F(box,age,year);
          ENDIF
          Z1=( M(box,age,year)+F1+TAG_LOSS(j) )          ; S1=EXP(-Z1*DELTA)
          Z2=M(oth,age,year)+F(oth,age,year)+TAG_LOSS(j) ; S2=DEXP(-Z2*DELTA) ! assume fish that crosses is fully mixed

          IF(MODEL_TYPE==1) THEN

!             {Diffusion model}
              IF(YEAR==TC(j)%Yr) THEN ;
!               {assume fish that had moved to the area where tagging occurred effectively recruited to that local population}
                N_TAGS(j,box,year)=TC(j)%N*TAG_SURV(j); N_TAGS(j,oth,year)=0
              ENDIF
              N_TAGS(j,box,year+1)=N_TAGS(j,box,year)*(1.0-T1)+ N_TAGS(j,oth,year)*T2
              N_TAGS(j,oth,year+1)=N_TAGS(j,box,year)*T1 + N_TAGS(j,oth,year)*(1.0-T2)
              RECAPTURES(j,box,year)=TAG_REPORT(j,box)*F1*N_TAGS(j,box,year+1)*(1.0-S1)/Z1;
              RECAPTURES(j,oth,year)=TAG_REPORT(j,oth)*F(oth,age,year)*N_TAGS(j,oth,year+1)*(1.0-S2)/Z2
              N_TAGS(j,box,year+1)=N_TAGS(j,box,year+1)*S1
              N_TAGS(j,oth,year+1)=N_TAGS(j,oth,year+1)*S2
              IF(YEAR==TC(j)%Ye) THEN
                RECOVERIES(j,box)=TAG_REPORT(j,0)*N_TAGS(j,box,year+1)
                RECOVERIES(j,oth)=TAG_REPORT(j,0)*N_TAGS(j,oth,year+1)
              ENDIF

            ELSE

!             {Overlap model}
              IF(YEAR==TC(j)%Yr) THEN ;
                !assume fish from each stock are tagged in proportion to their relative abundance in the zone (box)
                BOTTOM=N_STOCK(box,age,year)*(1-T(box,age,year))+N_STOCK(oth,age,year)*T(oth,age,year)
                N_TAGS(j,box,year)=TC(j)%N*TAG_SURV(j)*N_STOCK(box,age,year)*(1-T(box,age,year))/BOTTOM
                N_TAGS(j,oth,year)=TC(j)%N*TAG_SURV(j)*N_STOCK(oth,age,year)*T(oth,age,year)/BOTTOM
                T2=1.-(1.-T2)*T_proportion(T_mod_month(j))
              ENDIF
              RECAPTURES(j,box,year)=TAG_REPORT(j,box)*F1*(1.0-S1)/Z1*( N_TAGS(j,box,year)*(1.-T1) + N_TAGS(j,oth,year)*T2 )
              RECAPTURES(j,oth,year)=TAG_REPORT(j,oth)*F(oth,age,year)*(1.0-S2)/Z2*( N_TAGS(j,box,year)*T1 + N_TAGS(j,oth,year)*(1.-T2) )
              N_TAGS(j,box,year+1)=N_TAGS(j,box,year)*( (1.-T1)*S1 + T1*S2 )
              N_TAGS(j,oth,year+1)=N_TAGS(j,oth,year)*( (1.-T2)*S2 + T2*S1 )
              IF(YEAR==TC(j)%Ye) THEN
                RECOVERIES(j,box)=TAG_REPORT(j,0)*S1*( N_TAGS(j,box,year)*(1.-T1) + N_TAGS(j,oth,year)*T2 )
                RECOVERIES(j,oth)=TAG_REPORT(j,0)*S2*( N_TAGS(j,oth,year)*(1.-T2) + N_TAGS(j,box,year)*T1 )
              ENDIF

          ENDIF

        END DO

      END DO

      RETURN ; END

!------------------------------------------------------------------------------
      SUBROUTINE LIKELIHOOD(LIKE,X,NPARM)
!------------------------------------------------------------------------------
      USE DATUM ; USE STATISTICS ; USE PARAMETERS ; USE LOGLIKELIHOODS
      IMPLICIT NONE
      INTEGER G,Y,A,NPARM,I,B_REL,Y_REL,A_REL,FIRST_SSB,J
      REAL (KIND=8) :: ANSWER,PRED,OBS,LOGFACTORIAL,X(NPARM),LIKE,AVEREC,ND,SIGMA,SIGMA2,SUMPROB,DUMMY,SEL1,SEL2,MAXF(0:10000)
      EXTERNAL LOGFACTORIAL
      LIKELIHOODS=0 ; LIKE_EFFORT=0 ; LIKE_TAG=0 ;
      POSTERIORS=0  ; LIKE_Q=0 ; LIKE_F=0 ; LIKE_TERMINAL=0 ; LIKE_T=0 ; LIKE_M=0
                      LIKE_Tag_loss=0 ; LIKE_Tag_nomix=0 ; LIKE_Tag_report=0 ; LIKE_Tag_surv=0 ; LIKE_V=0
      CONSTRAINTS=0;  LIKE_RECRUITMENT=0 ;
      DUMMY=0
      IF(ICALL==0) N_DATA=0
      IF(NPARM==0) THEN ; EFFORT=EFFORT_DATA ; RECAPTURES=RECAPTURE_DATA ; RECOVERIES=RECOVERY_DATA; GOTO 99 ; END IF ! Deviance calculations require likelihood of perfect fit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
! Compute out of bounds penalties
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO BOX=1,NBOX
       DO Y=FIRSTYEAR,LASTYEAR
        DO A=FIRSTAGE,LASTAGE
            IF(MODEL_TYPE==1 .AND. N_AREA(box,a,y).LT.0) PENALTY=PENALTY+10+dABS(N_AREA(box,a,y))**1.5
            IF(MODEL_TYPE==2 .AND. N_STOCK(box,a,y).LT.0) THEN
              PENALTY=PENALTY+10+dABS(N_STOCK(box,a,y))**1.5
            ENDIF
            IF(MODEL_TYPE==1 .AND. N_AREA(box,a,y).gT.1.0d+11)  PENALTY=PENALTY+dABS(N_AREA(box,a,y)-1d+11)
            IF(MODEL_TYPE==2 .AND. N_STOCK(box,a,y).gT.1.0d+11) PENALTY=PENALTY+dABS(N_STOCK(box,a,y)-1d+11)
      END DO ; END DO ; END DO
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
! Compute constraints for deviations in recruitment between boxes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      IF(LINK_STOCK>0 .AND. NBOX>1) THEN
        DO Y=FIRSTYEAR,LASTYEAR
          IF(MODEL_TYPE.EQ.1) THEN ; OBS=N_AREA(1,firstage,y) ; PRED=N_AREA(2,firstage,y)
            ELSE ; OBS=N_STOCK(1,firstage,y) ; PRED=N_STOCK(2,firstage,y) ; ENDIF
          LIKE_RECRUITMENT(1)=LIKE_RECRUITMENT(1)-0.5*(LOG(RATIO_STOCK*PRED/OBS)/SIGMA_STOCK)**2.0 -LOG(SIGMA_STOCK)
        END DO
      ENDIF
      
99    DO 10 BOX=1,NBOX

        IF(NPARM==0) GOTO 999

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
! Compute constraints for auto-correlated stock-recruitment relationship
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(LINK_REC>0) THEN ! impose random walk penalty
          DO 11 Y=LASTYEAR-LINK_REC+1,LASTYEAR-1
            IF(MODEL_TYPE.EQ.1) THEN ; OBS=N_AREA(box,firstage,y) ; PRED=N_AREA(box,firstage,y+1) 
              ELSE ; OBS=N_STOCK(box,firstage,y) ; PRED=N_STOCK(box,firstage,y+1) ; ENDIF
            LIKE_RECRUITMENT(box)=LIKE_RECRUITMENT(box)-0.5*(LOG(PRED/OBS)/SIGMA_REC)**2 - LOG(SIGMA_REC)
11        CONTINUE
        ENDIF
        IF(PDF_STOCKRECRUIT.NE.0) THEN ! autocorrelated SR function
          IF(STOCK_RECRUIT(box,4).GE.1.) PENALTY=1.0D32 ! correlation can't exceed 1
          ANSWER=0.d0 ; ND=0.d0 ; AVEREC=0.d0
          FIRST_SSB=MAX(FIRSTYEAR,IGNORE_RECRUIT(1)-FIRSTAGE)
          DO 12 Y=IGNORE_RECRUIT(1),IGNORE_RECRUIT(2)
!           Y is the year of recruitment, i is the corresponding year for SSB
            I=Y-FIRSTAGE
            IF(I.GE.FIRST_SSB) THEN
              ND=ND+1.d0
              PRED=STOCK_RECRUIT(box,1)*SSB(box,i)/(STOCK_RECRUIT(box,2)+SSB(box,i))
              IF(MODEL_TYPE.EQ.1) THEN ; OBS=N_AREA(box,firstage,y) ; ELSE ; OBS=N_STOCK(box,firstage,y) ; ENDIF
              RECRUITS(box,y) = pred ; AVEREC=AVEREC+OBS
              if(ABS(PDF_STOCKRECRUIT).eq.1) then
                OBS=Dlog(OBS); PRED=Dlog(PRED)
              endif
              ETA_SR(box,y) = OBS-PRED
              IF(I.EQ.FIRST_SSB) THEN
                  ANSWER=ANSWER+(1.0-STOCK_RECRUIT(box,4)**2.0)*ETA_SR(box,y)**2.0
                ELSE
                  ANSWER=ANSWER+(ETA_SR(box,y)-STOCK_RECRUIT(box,4)*ETA_SR(box,y-1))**2.0
              ENDIF
            ENDIF
12        CONTINUE
          AVEREC=AVEREC/ND
          IF(PDF_STOCKRECRUIT>0) THEN ! do not estimate variance using concentrated likelihood approach
              SIGMA2=STOCK_RECRUIT(box,5)**2
            ELSE ! go ahead and estimate variance using concentrated likelihood approach
              SIGMA2=ANSWER/ND
              STOCK_RECRUIT(box,5)=SIGMA2
          ENDIF
          LIKE_RECRUITMENT(box)= -0.5*(ND*DLOG(SIGMA2)-DLOG(1.-STOCK_RECRUIT(box,4)**2)+ANSWER/SIGMA2)
        ENDIF
        CONSTRAINTS=CONSTRAINTS+LIKE_RECRUITMENT(BOX)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute constraints for F deviations from age to age and among last LINK_FT years
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
        DO A=FIRSTAGE,LASTAGE-1   ! the last age is handled by the F-ratio
          IF(SIGMA_TERMINAL(box,a)>0) LIKE_TERMINAL(box)=LIKE_TERMINAL(box)-0.5*(ETA_TERMINAL(box,a)/SIGMA_TERMINAL(box,a))**2.0 &
                                                         -LOG(SIGMA_TERMINAL(box,a))
        END DO
        IF(LINK_FT>1) THEN
          MAXF=-1.0
          DO Y=LASTYEAR-LINK_FT+1,LASTYEAR
            DO A=FIRSTAGE,LASTAGE-1; MAXF(Y)=MAX(MAXF(Y),F(box,a,y)) ; END DO
          END DO
          DO Y=LASTYEAR-LINK_FT+2,LASTYEAR
            DO A=LINK_YOUNGEST,LINK_OLDEST
              SEL1=F(box,a,y)/MAXF(Y) ; SEL2=F(box,a,y-1)/MAXF(Y-1)
              LIKE_TERMINAL(box)=LIKE_TERMINAL(box)-0.5*(LOG(sel1/sel2)/SIGMA_FT)**2.0 -LOG(SIGMA_FT)
            END DO
          END DO
        ENDIF
        CONSTRAINTS = CONSTRAINTS + LIKE_TERMINAL(box)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute log-likelihoods for indices of abundance and catchability deviations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

999     DO 20 G=1,NGEARS(box)
          IF(INT(CV_OVERIDE)==999 .and. nparm>0) THEN
            ! use concentrated likelihood with the nuisance parameter sigma replaced by its mle
            SIGMA=0 ; ND=0
            DO 21 Y=FIRSTYEAR,LASTYEAR
              IF(PDF_EFFORT(g,box)>1 .AND. EFFORT_DATA_STORE(g,box,y)>-9.0) THEN
                  SIGMA = SIGMA+(EFFORT_DATA(g,box,y)-EFFORT(g,box,y))**2.0 ; ND=ND+1
                ELSEIF(PDF_EFFORT(g,box)==1 .AND. EFFORT_DATA_STORE(g,box,y)>0) THEN
                  SIGMA = SIGMA+(DLOG(EFFORT_DATA(g,box,y))-DLOG(EFFORT(g,box,y)))**2.0 ; ND=ND+1.0
              ENDIF
21          CONTINUE
          ENDIF
          DO 22 Y=FIRSTYEAR,LASTYEAR
            IF((PDF_EFFORT(g,box)>1.AND.EFFORT_DATA_STORE(g,box,y)>-9).OR. &
              (PDF_EFFORT(g,box)==1.AND.EFFORT_DATA_STORE(g,box,y)>0)) THEN
              IF(INT(CV_OVERIDE)==999 .and. nparm>0) THEN ! when computing deviance we cannot estimate sigma, must use MLE
                IF(PDF_EFFORT(g,box)==2) THEN ; SIGMA_EFFORT(g,box,y)=-SQRT(SIGMA/ND) ! subroutine formulae uses -std dev.
                 ELSEIF(PDF_EFFORT(g,box)==1) THEN ; SIGMA_EFFORT(g,box,y)=SQRT(EXP(SIGMA/ND)-1) ; ENDIF ! uses +cv
              ENDIF
              CALL FORMULAE(EFFORT_DATA(g,box,y),EFFORT(g,box,y),SIGMA_EFFORT(g,box,y),SCALE_VARIANCE(g,box),PDF_EFFORT(g,box), &
                            ANSWER,CONSTANTS)
              LIKE_EFFORT(g,box) = LIKE_EFFORT(g,box)+ANSWER
              IF(ICALL==0) N_DATA=N_DATA+1
              IF(SIGMA_Q(g,box,y)>0 .AND. NPARM>0) &
                LIKE_Q(g,box)=LIKE_Q(g,box)-0.5*(eta_q(g,box,y)/SIGMA_Q(g,box,y))**2.0 - LOG(SIGMA_Q(g,box,y))
            ENDIF
22        CONTINUE
20      LIKELIHOODS = LIKELIHOODS + LIKE_EFFORT(g,box) + LIKE_Q(g,box)

        IF(NPARM==0) GOTO 9999 ! deviance pertains only to likelihoods

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute posterior for F-ratio deviations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
        DO 30 Y=FIRSTYEAR,LASTYEAR
         IF(SIGMA_F(box,y)>0) LIKE_F(box)=LIKE_F(box)-0.5*(ETA_F(box,y)/SIGMA_F(box,y))**2.0 -LOG(SIGMA_F(box,y))
30      CONTINUE
        POSTERIORS = POSTERIORS + LIKE_F(box)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute posterior for M deviations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
        DO 40 a=FIRSTAGE,LASTAGE
         IF(SIGMA_M(box,a,firstyear)>0) LIKE_M(box)=LIKE_M(box)-0.5*(ETA_M(box,a,firstyear)/SIGMA_M(box,a,firstyear))**2.0  &
                                                    -LOG(SIGMA_M(box,a,firstyear))
40      CONTINUE
        POSTERIORS = POSTERIORS + LIKE_M(box)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute posterior for T deviations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
        DO 50 a=FIRSTAGE,LASTAGE
         IF(SIGMA_T(box,a,firstyear)>0) LIKE_T(box)=LIKE_T(box)-0.5*(ETA_T(box,a,firstyear)/SIGMA_T(box,a,firstyear))**2.0 &
                                                               -LOG(SIGMA_T(box,a,firstyear))
50      CONTINUE
        POSTERIORS = POSTERIORS + LIKE_T(box)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute posterior for Variance-scaling deviations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
        DO G=1,NGEARS(box)
         IF(SIGMA_SCALE_VARIANCE(g,box)>0) LIKE_V(box)=LIKE_V(box) &
           -0.5*(ETA_SCALE_VARIANCE(g,box)/SIGMA_SCALE_VARIANCE(g,box))**2.0 - LOG(SIGMA_SCALE_VARIANCE(g,box))
        END DO
        POSTERIORS = POSTERIORS + LIKE_V(box)

10    CONTINUE


      IF(PDF_TAG>0) THEN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute posterior for tag mortality deviations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
         DO a=1,N_SURV_PARS
          IF(SIGMA_TAG_SURV(a)>0) LIKE_TAG_SURV=LIKE_TAG_SURV  &
                                  -0.5*( ETA_TAG_SURV(a)/SIGMA_TAG_SURV(a))**2.0 - LOG(SIGMA_TAG_SURV(a) )
         END DO
         POSTERIORS = POSTERIORS + LIKE_TAG_SURV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute posterior for tag shedding deviations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
        DO a=1,N_SHED_PARS
         IF(SIGMA_Tag_loss(a)>0) LIKE_Tag_loss=LIKE_Tag_loss &
                                 -0.5*( ETA_Tag_loss(a)/SIGMA_Tag_loss(a))**2.0 - LOG(SIGMA_TAG_LOSS(a) )
        END DO
        POSTERIORS = POSTERIORS + LIKE_Tag_loss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute posterior for tag reporting deviations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
        DO y=1,N_REPT_PARS
         IF(SIGMA_Tag_report(y)>0) LIKE_Tag_report=LIKE_Tag_report &
                                   -0.5*(ETA_Tag_report(y)/SIGMA_Tag_report(y))**2.0 - LOG(SIGMA_TAG_REPORT(y))
        END DO
        POSTERIORS = POSTERIORS + LIKE_Tag_report

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute posterior for tag-nonmixing deviations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
        DO a=1,N_MIX_PARS
         IF(SIGMA_Tag_nomix(a)>0) LIKE_Tag_nomix=LIKE_Tag_nomix &
                                  -0.5*(ETA_Tag_nomix(a)/SIGMA_Tag_nomix(a))**2.0 - LOG(SIGMA_TAG_nomix(a))
        END DO
        POSTERIORS = POSTERIORS + LIKE_Tag_nomix

      END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute log-likelihoods for tagging data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
9999  IF(PDF_TAG > 0) THEN

        DO J=1,N_TAG_COHORTS

          IF(TC(j)%N<=0) CYCLE

          SUMPROB=0
          DO BOX=1,NBOX
            DO YEAR=TC(j)%Yr,TC(j)%Ye
              PRED=RECAPTURES(j,box,year) ; SUMPROB=SUMPROB+PRED
              IF(PDF_TAG==3) PRED=PRED/TC(j)%N
              CALL FORMULAE(RECAPTURE_DATA(j,box,year),pred,TC(j)%sigma,DUMMY,PDF_TAG,ANSWER,CONSTANTS)
              IF(PDF_TAG==3) ANSWER=ANSWER/TC(j)%sigma**2.0
              LIKE_TAG(TC(j)%Br)=LIKE_TAG(TC(j)%Br)+ANSWER
              IF(ICALL==0) N_DATA=N_DATA+1
            END DO ! year
            !write(*,*) J,Box,RECOVERY_DATA(j,box),(RECAPTURE_DATA(j,box,year),Year=TC(j)%Yr,TC(j)%Ye) 
            !write(*,*) J,Box,RECOVERies(j,box),(RECAPTUREs(j,box,year),Year=TC(j)%Yr,TC(j)%Ye) 
            IF(RECOVERY_DATA(j,box)>=0) THEN
              PRED=RECOVERIES(j,box) ; SUMPROB=SUMPROB+PRED
              IF(PDF_TAG==3) PRED=PRED/TC(j)%N
              CALL FORMULAE(RECOVERY_DATA(j,box),pred,TC(j)%sigma,DUMMY,PDF_TAG,ANSWER,CONSTANTS)
              IF(PDF_TAG==3) ANSWER=ANSWER/TC(j)%sigma**2.0
              LIKE_TAG(TC(j)%Br)=LIKE_TAG(TC(j)%Br)+ANSWER
              IF(ICALL==0) N_DATA=N_DATA+1
            END IF
          END DO ! box
          PRED=TC(j)%N-SUMPROB
          IF(PDF_TAG==3) PRED=PRED/TC(j)%N
          OBS=TC(j)%N-SUMRDATA(j)
          CALL FORMULAE(OBS,pred,TC(j)%sigma,DUMMY,PDF_TAG,ANSWER,CONSTANTS)
          IF(PDF_TAG==3) THEN
            IF(CONSTANTS=='Y') ANSWER=ANSWER+LOGFACTORIAL(TC(j)%N-1)
            ANSWER=ANSWER/TC(j)%sigma**2.0
          ENDIF
          LIKE_TAG(TC(j)%Br)=LIKE_TAG(TC(j)%Br)+ANSWER
          IF(ICALL==0) N_DATA=N_DATA+1

        END DO ! tag cohorts

        DO BOX=1,NBOX ; LIKELIHOODS=LIKELIHOODS+LIKE_TAG(box) ; END DO

      ENDIF

      LIKE=LIKELIHOODS+POSTERIORS+CONSTRAINTS

      RETURN ; END

!-----------------------------------------------------------------
      SUBROUTINE OVERLAP_CATCH(np,Xvec,Fvec,DF,LOOP)
!-----------------------------------------------------------------
      USE DATUM, ONLY : CATCH_DATA ; USE STATISTICS, ONLY : N_STOCK,M,T,YEAR,AGE
      IMPLICIT NONE
      INTEGER LOOP,NP
      REAL (KIND = 8) ::  X,Y,F,G,FX,FY,GX,GY,a,b,c,d,K,L,v,u,t1,u1,m1,n1,z1,z2,s1,s2,s3,s4,s5,tf,ti, &
                          ab1,cd1,aub,ctd,ab2,cd2,C3,C4,ATU2BU1T1,tu1ab1,tu1aub,ct2udt1u1,a5,a6
      REAL (KIND = 8) ::  xvec(np),fvec(np),df(np,np)
      
      IF(LOOP<=1) THEN
        K=CATCH_DATA(1,AGE,YEAR) ; L=CATCH_DATA(2,AGE,YEAR)
        v=T(1,AGE,YEAR) ; u=T(2,AGE,YEAR) ; m1=M(1,AGE,YEAR) ; n1=M(2,AGE,YEAR)
        a=N_STOCK(1,AGE+1,YEAR+1)*(1-v) ; b=N_STOCK(2,AGE+1,YEAR+1)*u; c=N_STOCK(2,AGE+1,YEAR+1)*(1-u); d=N_STOCK(1,AGE+1,YEAR+1)*v
        t1=v-1 ; u1=u-1 ; tf=v*(3*u-1)-u+1 ; ti=v*(3*u-2) ; ab1=a*u1-b*v  ; cd1=c*t1-d*u ; aub=a*u-b*t1
        ctd=c*v-d*u1 ; ab2=a*t1*u1**2-b*v**2*u ; cd2=c*v**2*u+d*(1-v)*u1**2 ; tu1ab1=v*u1*ab1 ; tu1aub=v*u1*aub
        atu2bu1t1=a*v*u*u-b*u1*t1*t1 ; ct2udt1u1=c*v**2*u-d*t1*u1**2 ; c3=c*t1*tf+d*u*(2*u1-ti)
        c4=c*v*(ti-2*u1)-d*u1*tf ; a5=m1*(a*u*(ti-2*u1)-b*t1*tf) ;
        IF(LOOP<1) THEN ; XVEC(1)=K*(M1+0.05)/(EXP(0.05+M1)-1)/(a+b); XVEC(2)=l*(M1+0.05)/(EXP(0.05+M1)-1)/(c+d)
          XVEC(1)= MIN(3.0,MAX(XVEC(1),1.D-7)) ; XVEC(2)= MIN(3.0,MAX(XVEC(2),1.D-7)) ; RETURN
        ENDIF
      ENDIF

      if(xvec(1).gt.3) then ; xvec(1)=3.0 ; elseif(xvec(1).lt.1.0D-7) then ; xvec(1)=1.0D-7 ; endif
      if(xvec(2).gt.3) then ; xvec(2)=3.0 ; elseif(xvec(2).lt.1.0D-7) then ; xvec(2)=1.0D-7 ; endif
      X=XVEC(1) ; Y=XVEC(2)
      a6=m1*(2*x*tu1aub-a*u1*tf+b*v*(ti-2*u1))
      z1=x+m1 ; z2 =y+n1 ; S1 = EXP(z1) ; S2 = EXP(z2) ; S3=S1*S2 ; S4=S1*S3 ; S5=S2*S3

      F = k-x*(1-EXP(-z1))*(a/(v/S2-t1/S1)+b/(u/S1-u1/S2))/z1

      FX = -S2*(m1*v*S1**4*u1*ab1-S1**3*(S2*(x**2*ab2+m1*x*ab2+m1*(a*u1*tf+b*v*(2*u1-ti)))-tu1ab1*(x**2+m1*x-m1))+S4*(S2*(2*u* &  
           x**2*t1*ab1+2*m1*u*x*t1*ab1+a5)-2*v*x**2*u1*aub-a6)-S5*(u*S2*t1*(x**2+m1*x+m1)* &
           (a*u+b*(1-v))-x**2*atu2bu1t1-m1*x*atu2bu1t1+a5)+m1* &
           u*S2**3*t1*aub)/(z1**2*(v*S1+S2*(1-v))**2*(S1*u1-u*S2)**2)

      FY = x*S3*(1-S1)*(S1*S1*tu1ab1-2*S3*tu1aub+    &
           S2*S2*atu2bu1t1)/(z1*(v*S1-S2*t1)**2*(S1*u1-u*S2)**2)

      G = l-y*(1-EXP(-z2))*(c/(u/S1-u1/S2)+d/(v/S2-t1/S1))/z2

      GX = y*S3*(1-S2)*(S1*S1*ct2udt1u1-2*u*S3*t1*ctd+  &
           u*S2*S2*t1*cd1)/(z2*(v*S1-S2*t1)**2*(S1*u1-u*S2)**2)

      GY = S1*(v*S1**3*u1*ctd*(S2*(y**2+n1*y+n1)-n1)-S4*(S2*(2*v*y**2*u1*cd1     &
           +2*n1*v*y*u1*cd1+n1*(c*v*(ti-2*u1)-d*u1*tf))+y**2*cd2+n1*y*ct2udt1u1-n1*c4)+S5*(S2*(y**2*(c*t1**2*u1-d*v*u**2)+n1 &
           *y*(c*t1**2*u1-d*v*u**2)+n1*(c*t1*tf+d*u*(2*(u-1)-ti)))+2*u*y**2*t1*ctd+2*n1*u*y*t1*(c*v-   &
           d*u1)-n1*c3)+u*S2**3*(1-v)*cd1*(n1*S2+y**2+n1*y-n1))/(z2**                     &
           2*(v*S1-S2*t1)**2*(S1*u1-u*S2)**2)

      FVEC(1)=F ; FVEC(2)=G ; DF(1,1)=FX ; DF(1,2)=FY; DF(2,1)=GX ; DF(2,2)=GY

      RETURN ; END      

!-----------------------------------------------------------------
      SUBROUTINE OVERLAP_PLUS(np,Xvec,Fvec,DF,LOOP)
!-----------------------------------------------------------------
      USE DATUM, ONLY : CATCH_DATA ; USE STATISTICS, ONLY : N_STOCK,M,T,YEAR,AGE,PLUSGROUP,F_RATIO
      IMPLICIT NONE
      INTEGER LOOP,NP
      REAL (KIND = 8) ::  X,Y,F,G,FX,FY,GX,GY,a,b,c,d,K,J,O,P,R,S,N1,N2,XVEC(np),fvec(np),df(np,np),MFEXP

      EXTERNAL MFEXP
      
      IF(LOOP<=1) THEN
        N1=N_STOCK(1,PLUSGROUP,YEAR+1) ; N2=N_STOCK(2,PLUSGROUP,YEAR+1)
        IF(LOOP<1) THEN  
          K=CATCH_DATA(1,PLUSGROUP-1,YEAR) ; o=M(1,PLUSGROUP-1,YEAR)+0.01 ; C=T(1,PLUSGROUP-1,YEAR) ; A=1-C
          J=CATCH_DATA(2,PLUSGROUP-1,YEAR) ; p=M(2,PLUSGROUP-1,YEAR)+0.01 ; B=T(2,PLUSGROUP-1,YEAR) ; D=1-B
          XVEC(1)=MIN(3.0,MAX(k*o/(1-exp(-o))/(n1*a+n2*b),1.D-7)) ; XVEC(2)=MIN(3.0,MAX(j*p/(1-exp(-p))/(n2*d+n1*c),1.D-7)) ;  RETURN
        ENDIF
      ENDIF

      X=XVEC(1) ; Y=XVEC(2)
      K=CATCH_DATA(1,PLUSGROUP-1,YEAR) ; o=M(1,PLUSGROUP-1,YEAR) ; C=T(1,PLUSGROUP-1,YEAR) ; A=1-C
      J=CATCH_DATA(2,PLUSGROUP-1,YEAR) ; p=M(2,PLUSGROUP-1,YEAR) ; B=T(2,PLUSGROUP-1,YEAR) ; D=1-B
      
      F = (c*MFEXP(-y-p)+a*MFEXP(-x-o))*((d*k*(x+o)/x/(1-MFEXP(-x-o))-b*j*(y+p)/y/(1-MFEXP(-y-p)))/(a*d-b*c))
      
      G = (b*MFEXP(-x-o)+d*MFEXP(-y-p))*((a*j*(y+p)/y/(1-MFEXP(-y-p))-c*k*(x+o)/x/(1-MFEXP(-x-o)))/(a*d-b*c))

      FX = MFEXP(-x-y-o-p)*(c*d*k*o*y*MFEXP(3*x+3*o)*(MFEXP(y+p)-1)-MFEXP(2*x+2*o &
     )*(a*MFEXP(2*y+2*p)*(x**2*(y*(b*j-d*k)+b*j*p)-d*k*o*x*y-d*k*o*y)+d*      &
     k*y*MFEXP(y+p)*(x**2*(a-c)+o*x*(a-c)+o*(a+c))+c*d*k*y*(x**2+o*x-o))      &
     +a*MFEXP(x+y+o+p)*(MFEXP(y+p)*(2*b*j*x**2*(y+p)-d*k*o*y)+d*k*o*y)-a*      &
     b*j*x**2*MFEXP(2*y+2*p)*(y+p))/(x**2*y*(MFEXP(x+o)-1)**2*(MFEXP(y+p)-      &
     1)*(b*c-a*d))

      FY = MFEXP(-x-y-o-p)*(c*MFEXP(2*x+2*o)*(MFEXP(2*y+2*p)*(x*(y**2*(b*j-d*k)+  &
     b*j*p*y+b*j*p)-d*k*o*y**2)+MFEXP(y+p)*(x*(2*d*k*y**2-b*j*p)+2*d*k*o       &
     *y**2)-d*k*y**2*(x+o))+b*j*x*MFEXP(x+y+o+p)*(a*p*MFEXP(2*y+2*p)+MFEXP       &
     (y+p)*(y**2*(a-c)+p*y*(a-c)-p*(a+c))+c*p)-a*b*j*x*MFEXP(2*y+2*p)*(p       &
     *MFEXP(y+p)+y**2+p*y-p))/(x*y**2*(MFEXP(x+o)-1)*(MFEXP(y+p)-1)**2*(a*       &
     d-b*c))

      GX = MFEXP(-x-y-o-p)*(c*d*k*o*y*MFEXP(3*x+3*o)*(MFEXP(y+p)-1)-MFEXP(2*x+2*o &
     )*(b*MFEXP(2*y+2*p)*(x**2*(y*(a*j-c*k)+a*j*p)-c*k*o*x*y-c*k*o*y)+c*       &
     k*y*MFEXP(y+p)*(x**2*(b-d)+o*x*(b-d)+o*(b+d))+c*d*k*y*(x**2+o*x-o))       &
     +b*MFEXP(x+y+o+p)*(MFEXP(y+p)*(2*a*j*x**2*(y+p)-c*k*o*y)+c*k*o*y)-a*       &
     b*j*x**2*MFEXP(2*y+2*p)*(y+p))/(x**2*y*(MFEXP(x+o)-1)**2*(MFEXP(y+p)-       &
     1)*(a*d-b*c))

      GY = MFEXP(-x-y-o-p)*(d*MFEXP(2*x+2*o)*(MFEXP(2*y+2*p)*(x*(y**2*(a*j-c*k)+  &
     a*j*p*y+a*j*p)-c*k*o*y**2)+MFEXP(y+p)*(x*(2*c*k*y**2-a*j*p)+2*c*k*o       &
     *y**2)-c*k*y**2*(x+o))+a*j*x*MFEXP(x+y+o+p)*(b*p*MFEXP(2*y+2*p)+MFEXP       &
     (y+p)*(y**2*(b-d)+p*y*(b-d)-p*(b+d))+d*p)-a*b*j*x*MFEXP(2*y+2*p)*(p       &
     *MFEXP(y+p)+y**2+p*y-p))/(x*y**2*(MFEXP(x+o)-1)*(MFEXP(y+p)-1)**2*(b*       &
     c-a*d))

      K=CATCH_DATA(1,PLUSGROUP,YEAR) ; o=M(1,PLUSGROUP,YEAR) ; r=F_RATIO(1,YEAR) ; C=T(1,PLUSGROUP,YEAR) ; A=1-C
      J=CATCH_DATA(2,PLUSGROUP,YEAR) ; p=M(2,PLUSGROUP,YEAR) ; s=F_RATIO(2,YEAR) ; B=T(2,PLUSGROUP,YEAR) ; D=1-B
      
      F = F + (c*MFEXP(-s*y-p)+a*MFEXP(-r*x-o))*((d*k*(r*x+o)/(r*x)/(1-MFEXP(-r*x-o))-b*j*(s*y+p)/(s*y)/(1-MFEXP(-s*y-p)))/(a*d-b*c)) 

      G = G + (b*MFEXP(-r*x-o)+d*MFEXP(-s*y-p))*((a*j*(s*y+p)/(s*y)/(1-MFEXP(-s*y-p))-c*k*(r*x+o)/(r*x)/(1-MFEXP(-r*x-o)))/(a*d-b*c)) 

     FX = FX + MFEXP(-r*x-s*y-o-p)*(c*d*k*o*s*y*MFEXP(3*r*x+3*o)*(MFEXP(s*y+p)-1)-   &
     MFEXP(2*r*x+2*o)*(a*MFEXP(2*s*y+2*p)*(r**2*x**2*(s*y*(b*j-d*k)+b*j*p)-             &
     d*k*o*r*s*x*y-d*k*o*s*y)+d*k*s*y*MFEXP(s*y+p)*(r**2*x**2*(a-c)+o*r*             &
     x*(a-c)+o*(a+c))+c*d*k*s*y*(r**2*x**2+o*r*x-o))+a*MFEXP(r*x+s*y+o+p             &
     )*(MFEXP(s*y+p)*(2*b*j*r**2*x**2*(s*y+p)-d*k*o*s*y)+d*k*o*s*y)-a*b*             &
     j*r**2*x**2*MFEXP(2*s*y+2*p)*(s*y+p))/(r*s*x**2*y*(MFEXP(r*x+o)-1)**             &
     2*(MFEXP(s*y+p)-1)*(b*c-a*d))

      FY = FY + MFEXP(-r*x-s*y-o-p)*(c*MFEXP(2*r*x+2*o)*(MFEXP(2*s*y+2*p)*(r*x*(s**2*   &
     y**2*(b*j-d*k)+b*j*p*s*y+b*j*p)-d*k*o*s**2*y**2)+MFEXP(s*y+p)*(r*x*             &
     (2*d*k*s**2*y**2-b*j*p)+2*d*k*o*s**2*y**2)-d*k*s**2*y**2*(r*x+o))+             &
     b*j*r*x*MFEXP(r*x+s*y+o+p)*(a*p*MFEXP(2*s*y+2*p)+MFEXP(s*y+p)*(s**2*y             &
     **2*(a-c)+p*s*y*(a-c)-p*(a+c))+c*p)-a*b*j*r*x*MFEXP(2*s*y+2*p)*(p*             &
     MFEXP(s*y+p)+s**2*y**2+p*s*y-p))/(r*s*x*y**2*(MFEXP(r*x+o)-1)*(MFEXP(s             &
     *y+p)-1)**2*(a*d-b*c))

      GX = GX + MFEXP(-r*x-s*y-o-p)*(c*d*k*o*s*y*MFEXP(3*r*x+3*o)*(MFEXP(s*y+p)-1)-  &
     MFEXP(2*r*x+2*o)*(b*MFEXP(2*s*y+2*p)*(r**2*x**2*(s*y*(a*j-c*k)+a*j*p)-             &
     c*k*o*r*s*x*y-c*k*o*s*y)+c*k*s*y*MFEXP(s*y+p)*(r**2*x**2*(b-d)+o*r*             &
     x*(b-d)+o*(b+d))+c*d*k*s*y*(r**2*x**2+o*r*x-o))+b*MFEXP(r*x+s*y+o+p             &
     )*(MFEXP(s*y+p)*(2*a*j*r**2*x**2*(s*y+p)-c*k*o*s*y)+c*k*o*s*y)-a*b*             &
     j*r**2*x**2*MFEXP(2*s*y+2*p)*(s*y+p))/(r*s*x**2*y*(MFEXP(r*x+o)-1)**             &
     2*(MFEXP(s*y+p)-1)*(a*d-b*c))

      GY = GY + MFEXP(-r*x-s*y-o-p)*(d*MFEXP(2*r*x+2*o)*(MFEXP(2*s*y+2*p)*(r*x*(s**2*  &
     y**2*(a*j-c*k)+a*j*p*s*y+a*j*p)-c*k*o*s**2*y**2)+MFEXP(s*y+p)*(r*x*            &
     (2*c*k*s**2*y**2-a*j*p)+2*c*k*o*s**2*y**2)-c*k*s**2*y**2*(r*x+o))+            &
     a*j*r*x*MFEXP(r*x+s*y+o+p)*(b*p*MFEXP(2*s*y+2*p)+MFEXP(s*y+p)*(s**2*y            &
     **2*(b-d)+p*s*y*(b-d)-p*(b+d))+d*p)-a*b*j*r*x*MFEXP(2*s*y+2*p)*(p*            &
     MFEXP(s*y+p)+s**2*y**2+p*s*y-p))/(r*s*x*y**2*(MFEXP(r*x+o)-1)*(MFEXP(s            &
     *y+p)-1)**2*(b*c-a*d))
      
     F = N1 - F ; G = N2 - G ; FX = -1*FX ; FY = -1*FY ; GX=-1*GX ; GY = -1*GY
     FVEC(1)=F ; FVEC(2)=G ; DF(1,1)=FX ; DF(1,2)=FY; DF(2,1)=GX ; DF(2,2)=GY
     
     RETURN ; END


!-----------------------------------------------------------------
      SUBROUTINE DIFFUSION_CATCH(np,Xvec,Fvec,DF,LOOP)
!-----------------------------------------------------------------
      USE DATUM, ONLY : CATCH_DATA ; USE STATISTICS, ONLY : N_AREA,M,T,YEAR,AGE,BOX
      IMPLICIT NONE
      INTEGER LOOP,NP
      REAL (KIND = 8) ::  X,m1,z,S,a,K,xvec(np),fvec(np),df(np,np)
!      if(loop>0) write(*,*) loop,'catch',x
      IF(LOOP<=1) THEN
       K=CATCH_DATA(BOX,AGE,YEAR) ; m1=M(BOX,AGE,YEAR); a=N_AREA(BOX,AGE+1,YEAR+1)
       IF(LOOP<1) THEN ; X=(m1+0.05)*K/a/(EXP(m1+0.05)-1) ; XVEC(1)= MIN(3.0,MAX(x,1.0D-7)) ; RETURN ; ENDIF
      ENDIF
!     x=min(20.0,max(x,0.00000001))
      X=XVEC(1)
      Z=m1+x ; S = EXP(z) ; Fvec(1) = k-x*a*(S-1)/z 
      dF(1,1) =a*(m1-S*(z*x+m1))/z**2
      
      RETURN ; END

!-----------------------------------------------------------------
      SUBROUTINE DIFFUSION_PLUS(np,Xvec,Fvec,DF,LOOP)
!-----------------------------------------------------------------
      USE DATUM, ONLY : CATCH_DATA ; USE STATISTICS, ONLY : N_AREA,M,T,YEAR,AGE,BOX,F_RATIO,PLUSGROUP
      IMPLICIT NONE
      INTEGER LOOP,NP
      REAL (KIND = 8) ::  X,Y,F,FX,xvec(np),fvec(np),df(np,np),O,p,z1,z2,S1,s2,a,b,r,K,bp,aor
      
!      if(loop.gt.0) write(*,*) loop,'plus',x
      IF(LOOP<=1) THEN
       A=CATCH_DATA(BOX,PLUSGROUP-1,YEAR) ; o=M(BOX,PLUSGROUP-1,YEAR) ; K=N_AREA(BOX,PLUSGROUP,YEAR+1)
       B=CATCH_DATA(BOX,PLUSGROUP,YEAR) ; p=M(BOX,PLUSGROUP,YEAR) ; r=F_RATIO(BOX,YEAR) 
       IF(LOOP<1) THEN ; X=O/(K*(EXP(O+0.02)-1)/(A+B+1)) ; XVEC(1)=MIN(3.0,MAX(x,1.D-7)) ; RETURN ; ENDIF
      ENDIF
!      x=min(20.0,max(x,0.000001))
      X=XVEC(1)
      z1=x+o ; S1 = EXP(z1) ; y=x*r ; z2=y+p ; S2 = EXP(z2); bp=b*p ; aor=a*o*r
      F = k-a*z1/(x*(S1-1))-b*z2/(r*x*(S2-1))

      FX = (b*S1*S1*(S2*(r**2*x**2+p*r*x+p)-p)+S1*(a*r*S2*S2*(x**2+o*x+o)-2*S2*(r*x**2*(a+b*r)+r*x*(   &
           a*o+bp)+aor+bp)+a*r*x**2+aor*x+aor+2*bp)-aor*S2*S2+S2*(b*r**2*x**2+bp*r*x+2*aor+bp) &
           -aor-bp)/(r*x**2*(S1-1)**2*(S2-1)**2)

      FVEC(1)=F ; DF(1,1)=FX
      
      RETURN ; END

!--------------------------------------------------------------------------      
      SUBROUTINE RESULTS(THE_END)
!--------------------------------------------------------------------------      
      USE STATISTICS; USE DATUM; USE PARAMETERS ; USE LOGLIKELIHOODS
      IMPLICIT NONE
      INTEGER A,G,Y,OTH,I,THE_END,Y_REL,A_REL,B_REC,J
      REAL (KIND=8) :: SIG2,OBS,PRED,X(1000),RAN1,DEV_EFFORT(Gs,Bs),DEV_TAG(Bs),DEV_EFF,DEVIANCE,mOBS(1000),mPRED(1000), &
                       AIC,AICC,BIC,FULL,LESSCONSTANTS
      CHARACTER ch2*21
      CHARACTER (LEN=10) DUMMY(3)
      EXTERNAL OBJECTIVE, RAN1
      CALL DATE_AND_TIME(DUMMY(1),DUMMY(2),DUMMY(3),DATE_VALUES)
      OTH=2 ; EFFORT_DISCREPANCY=0 ; SUM_EFFORT_DISCREPANCY=0 ; DEV_EFF=0 ; DISCREPANCY=0

!     Compute probability of the data given the data to use for computing deviance
      I=0 ; CALL OBJECTIVE(X,SIG2,I) ; DEVIANCE=LIKELIHOODS
      DO BOX=1,NBOX
       DEV_TAG(box)=LIKE_TAG(box)
       DO G=1,NGEARS(BOX)
         DEV_EFFORT(g,box)=LIKE_EFFORT(g,box) ; DEV_EFF=DEV_EFF+LIKE_EFFORT(g,box)
      END DO ; END DO

!     Compute information criteria using full likelihood expression (with constants)
      CONSTANTS='Y'; CALL OBJECTIVE(PARM_EST,PARM_EST(PARM_KEY(0)+1),PARM_KEY(0))
      FULL=-LIKELIHOODS-POSTERIORS-CONSTRAINTS+PENALTY
      AIC  = 2.0*(FULL+DBLE(PARM_KEY(0)+N_Qs))
      AICc = 2.0*FULL+2.0*DBLE(PARM_KEY(0)+N_Qs)*(1.0+DBLE(PARM_KEY(0)+N_Qs+1)/DBLE(N_DATA-PARM_KEY(0)-N_Qs-1.0))
      BIC  = 2.0*FULL+(PARM_KEY(0)+N_Qs)*LOG(DBLE(N_DATA))

!     Drop constants in likelihood expressions for other statistics
      CONSTANTS='N'; CALL OBJECTIVE(PARM_EST,PARM_EST(PARM_KEY(0)+1),PARM_KEY(0))
      LESSCONSTANTS=-LIKELIHOODS-POSTERIORS-CONSTRAINTS+PENALTY
      DEVIANCE=-2.0*(LIKELIHOODS-DEVIANCE)
      CALL WRITE_PARAMETERS
      WRITE(*,*) 'WRITING TO FILE ',INFILE(3)
      DO BOX=1,NBOX
       DO G=1,NGEARS(BOX)
        IF(PDF_EFFORT(G,box)>0) THEN
          DO Y=1,LASTYEAR
            IF((PDF_EFFORT(g,box)>1 .AND. EFFORT_DATA(g,box,y)>-9) .OR. (PDF_EFFORT(g,box)==1 .AND. EFFORT_DATA(g,box,y)>0)) THEN
               CALL GETVARIANCE(EFFORT_DATA(g,box,y),EFFORT(g,box,y),SIGMA_EFFORT(g,box,y),SCALE_VARIANCE(g,box), &
                                PDF_EFFORT(g,box),SIG2)
               CALL GETDISCREPANCY(EFFORT_DATA(g,box,y),EFFORT(g,box,y),SIG2,PDF_EFFORT(g,box),EFFORT_DISCREPANCY(g,box,y))
               SUM_EFFORT_DISCREPANCY(g,box)=SUM_EFFORT_DISCREPANCY(g,box)+EFFORT_DISCREPANCY(g,box,y)
               DISCREPANCY=DISCREPANCY+EFFORT_DISCREPANCY(g,box,y)
            ENDIF
          END DO
        ENDIF
      END DO ; END DO
      OPEN(7,file=infile(3),status='unknown')
      OPEN(17,file=infile(5),status='unknown')
      WRITE(7,'(16X,A40)') '****************************************'
      WRITE(7,'(16X,A40)') '                VPA-2BOX                '
      WRITE(7,'(16X,A40)') 'SUMMARY STATISTICS AND DIAGNOSTIC OUTPUT'
      WRITE(7,'(16X,A40)') '****************************************'
      WRITE(7,*)
      WRITE(7,'(A50)') INFILE(0)     ! Title of run
      A=DATE_VALUES(6)/10 ; G=DATE_VALUES(6)-A*10 ; Y=DATE_VALUES(1) ;
      WRITE(7,'(I2,A1,I1,I1,A2,I2,A1,A9,A1,I4)') DATE_VALUES(5),':',A,G,', ',DATE_VALUES(3),' ',MONTH_NAME(DATE_VALUES(2)),' ',Y
      WRITE(7,*)
      WRITE(7,*) '======================================================================='
      OBS=0 ; PRED=0
      DO 20  BOX=1,NBOX
       DEV_TAG(box)=-2.0*(LIKE_TAG(box)-DEV_TAG(box))
       DO 20  G=1,NGEARS(BOX)
        OBS=OBS+LIKE_Q(g,box) ; PRED=PRED+LIKE_EFFORT(g,box) ; DEV_EFFORT(g,box)=-2.0*(LIKE_EFFORT(g,box)-DEV_EFFORT(g,box))
20    CONTINUE
      IF((PARM_KEY(0)+N_Qs)/N_DATA>0.25) THEN
        WRITE(*,'(/)')
        WRITE(*,*) 'WARNING: The number of parameters (',PARM_KEY(0)+N_Qs,')'
        WRITE(*,*) '         exceeds 25% of the number of data points (',N_DATA,')'
        wRITE(*,*) '         Model performance is likely questionable and bootstrap'
        wRITE(*,*) '         estimates of variance will probably be biased low'
        WRITE(*,'(/)')
        WRITE(7,'(/)')
        WRITE(7,*) 'WARNING: The number of parameters (',PARM_KEY(0)+N_Qs,')'
        WRITE(7,*) '         exceeds 25% of the number of data points (',N_DATA,')'
        wRITE(7,*) '         Model performance is likely questionable and bootstrap'
        wRITE(7,*) '         estimates of variance will probably be biased low'
        WRITE(7,'(/)')
      ENDIF
      
      OPEN(97,FILE='mix-dump.txt',STATUS='UNKNOWN')
      WRITE(97,*) T(1,2,1),T(2,2,1),LESSCONSTANTS
      CLOSE(97)

      CALL WRITE_CHAR(7,26,'Total objective function =',SPACE(1)) ; CALL WRITE_REAL(7,12,LESSCONSTANTS,ENDL)
      CALL WRITE_CHAR(7,26,'      (with constants)   =',SPACE(1)) ; CALL WRITE_REAL(7,12,FULL,ENDL)
      CALL WRITE_CHAR(7,26,'Number of parameters (P) =',SPACE(1)) ; CALL WRITE_INT(7,12,PARM_KEY(0)+N_Qs,ENDL)
      CALL WRITE_CHAR(7,26,'Number of data points (D)=',SPACE(1)) ; CALL WRITE_INT(7,12,N_DATA,ENDL)
      CALL WRITE_CHAR(7,26,'AIC : 2*objective+2P     =',SPACE(1)) ; CALL WRITE_REAL(7,12,AIC,ENDL)
      CALL WRITE_CHAR(7,26,'AICc: 2*objective+2P(...)=',SPACE(1)) ; CALL WRITE_REAL(7,12,AICc,ENDL)
      CALL WRITE_CHAR(7,26,'BIC : 2*objective+Plog(D)=',SPACE(1)) ; CALL WRITE_REAL(7,12,BIC,ENDL)
      CALL WRITE_CHAR(7,26,'Chi-square discrepancy   =',SPACE(1)) ; CALL WRITE_REAL(7,12,DISCREPANCY,ENDL)
      WRITE(7,*)
      CALL WRITE_CHAR(7,26,'Loglikelihoods (deviance)=',SPACE(1)) ; CALL WRITE_REAL(7,12,LIKELIHOODS,SPACE(1))
           CALL WRITE_CHAR(7,1,'(',SPACE(0)) ; CALL WRITE_REAL(7,12,DEVIANCE,SPACE(0)) ; CALL WRITE_CHAR(7,1,')',ENDL) ;
      CALL WRITE_CHAR(7,26,'   effort data           =',SPACE(1)) ; CALL WRITE_REAL(7,12,pred,SPACE(1)) ;
           CALL WRITE_CHAR(7,1,'(',SPACE(0)) ; CALL WRITE_REAL(7,12,-2.0*(PRED-DEV_EFF),SPACE(0)) ; CALL WRITE_CHAR(7,1,')',ENDL) ;
      IF(PDF_TAG>0) THEN
        CALL WRITE_CHAR(7,26,'   tagging data          =',SPACE(1)) ; CALL WRITE_REAL(7,12,LIKE_TAG(1)+LIKE_TAG(2),SPACE(1))
           CALL WRITE_CHAR(7,1,'(',SPACE(0)); CALL WRITE_REAL(7,12,DEV_TAG(1)+DEV_TAG(2),SPACE(0)); CALL WRITE_CHAR(7,1,')',ENDL);
      ENDIF
      WRITE(7,*)
      CALL WRITE_CHAR(7,26,'Log-posteriors           =',SPACE(1)) ; CALL WRITE_REAL(7,12,POSTERIORS         ,ENDL)
      CALL WRITE_CHAR(7,26,'   catchability          =',SPACE(1)) ; CALL WRITE_REAL(7,12,obs                ,ENDL)
      CALL WRITE_CHAR(7,26,'   f-ratio               =',SPACE(1)) ; CALL WRITE_REAL(7,12,LIKE_F(1)+LIKE_F(2),ENDL)
      CALL WRITE_CHAR(7,26,'   natural mortality     =',SPACE(1)) ; CALL WRITE_REAL(7,12,LIKE_M(1)+LIKE_M(2),ENDL)
      CALL WRITE_CHAR(7,26,'   mixing coeff.         =',SPACE(1)) ; CALL WRITE_REAL(7,12,LIKE_T(1)+LIKE_T(2),ENDL)
      IF(PDF_TAG>0) THEN
       CALL WRITE_CHAR(7,26,'   initial tag survival  =',SPACE(1)); CALL WRITE_REAL(7,12,LIKE_TAG_SURV,ENDL)
       CALL WRITE_CHAR(7,26,'   tag shedding rate     =',SPACE(1)); CALL WRITE_REAL(7,12,LIKE_TAG_LOSS,ENDL)
       CALL WRITE_CHAR(7,26,'   tag reporting rate    =',SPACE(1)); CALL WRITE_REAL(7,12,LIKE_TAG_REPORT,ENDL)
       CALL WRITE_CHAR(7,26,'   tag nomixing factor   =',SPACE(1)); CALL WRITE_REAL(7,12,LIKE_TAG_nomix,ENDL)
      ENDIF
      WRITE(7,*)
      CALL WRITE_CHAR(7,26,'Constraints              =',SPACE(1));CALL WRITE_REAL(7,12,CONSTRAINTS                            ,ENDL)
      CALL WRITE_CHAR(7,26,'   terminal F            =',SPACE(1));CALL WRITE_REAL(7,12,LIKE_terminal(1)+LIKE_terminal(2)      ,ENDL)
      CALL WRITE_CHAR(7,26,'   stock-rec./sex ratio  =',SPACE(1));CALL WRITE_REAL(7,12,LIKE_recruitment(1)+LIKE_recruitment(2),ENDL)
      WRITE(7,*)
      CALL WRITE_CHAR(7,26,'Out of bounds penalty    =',SPACE(1)) ; CALL WRITE_REAL(7,12,-PENALTY,ENDL)
      WRITE(7,*) '======================================================================='
      Write(7,'(/)')

      DO 30 BOX=1,NBOX
        IF(NBOX>1) THEN
          IF(BOX==2) THEN ; OTH=1 ; ELSE ; OTH=2 ; ENDIF
          WRITE(7,'(50I1)') (BOX,I=1,50) ; WRITE(7,'(A20)') TITLE_FISHERY(box) ; WRITE(7,'(50I1)') (BOX,I=1,50) ; WRITE(7,'(/)')
        ENDIF
!       get stock abundance at the beginning of the year after last 
        DO 35 A=FIRSTAGE,EXPANDED_PLUSGROUP
          AGE=A; IF(AGE>PLUSGROUP) AGE=PLUSGROUP
          IF(MODEL_TYPE==1) THEN ! diffusion and no-mixing models
              N_AREA(box,a+1,lastyear+1)=(N_AREA(box,a,lastyear)*(1-T(box,age,lastyear))+N_AREA(oth,a,lastyear)*T(oth,age,lastyear))* &
                                  EXP(-F(box,age,lastyear)-M(box,age,lastyear))
              IF(A==EXPANDED_PLUSGROUP) N_AREA(box,a,lastyear+1)=N_AREA(box,a,lastyear+1)+N_AREA(box,a+1,lastyear+1)
            ELSE ! overlap model
              N_STOCK(box,a+1,lastyear+1)=N_STOCK(box,a,lastyear)*((1-T(box,age,lastyear))*EXP(-F(box,age,lastyear)-M(box,age,lastyear))+&
                                                     T(box,age,lastyear)*EXP(-F(oth,age,lastyear)-M(oth,age,lastyear)))
              IF(A==EXPANDED_PLUSGROUP) N_STOCK(box,a,lastyear+1)=N_STOCK(box,a,lastyear+1)+N_STOCK(box,a+1,lastyear+1)
          ENDIF
35      CONTINUE
!       get stock-specific (model 1) or area-specific (model 2) abundances
        IF(NBOX>1) THEN    
          IF(MODEL_TYPE==1) THEN ! get stock specific abundance
              CALL GET_STOCK
            ELSE ! get area specific abundance
            DO Y=FIRSTYEAR,LASTYEAR
              DO A=FIRSTAGE,LASTAGE
                AGE=A; IF(AGE>PLUSGROUP) AGE=PLUSGROUP
                N_AREA(box,a,y)=N_STOCK(box,a,y)*(1-T(box,age,y))+N_STOCK(oth,a,y)*T(oth,age,y)
              END DO
            END DO
          ENDIF
        ENDIF
        IF(NBOX>1) THEN
            write(7,'(A8,I1,A28,A50)') 'TABLE 1.',box,' FISHING MORTALITY RATE FOR ',TITLE_FISHERY(box)
          ELSE
            write(7,'(A36,A50)') 'TABLE 1. FISHING MORTALITY RATE FOR ',TITLE_FISHERY(box)
        ENDIF
        write(7,311) ('=======',A=FIRSTAGE,LASTAGE+1)
        write(17,*) 'F'
        write(7,301) (A,A=FIRSTAGE,LASTAGE)
        write(7,311) ('-------',A=FIRSTAGE,LASTAGE+1)
        DO Y=FIRSTYEAR,LASTYEAR ; write(7,300) Y+DISPLAYYEAR-1,(F(box,a,y),A=FIRSTAGE,LASTAGE) ; END DO
        DO Y=1,99
          IF(Y<=LASTYEAR) THEN ; write(17,300) Y+DISPLAYYEAR-1,(F(box,a,y),A=FIRSTAGE,LASTAGE) ; ELSE ; WRITE(17,*) ; ENDIF
        END DO
        write(7,311) ('=======',A=FIRSTAGE,LASTAGE+1)
        WRite(7,'(/)')
        IF(NBOX>1) THEN
            write(7,'(A8,I1,A54,A50)') 'TABLE 2.',box,' ABUNDANCE AT THE BEGINNING OF THE YEAR [BY AREA] FOR ',TITLE_FISHERY(box)
          ELSE
            write(7,'(A62,A50)') 'TABLE 2. ABUNDANCE AT THE BEGINNING OF THE YEAR [BY AREA] FOR ',TITLE_FISHERY(box)
        ENDIF
        write(7,211) ('=============',A=FIRSTAGE,EXPANDED_PLUSGROUP+1)
        write(17,*) 'N by area'
        write(7,*)        
        write(7,201) (A,A=FIRSTAGE,EXPANDED_PLUSGROUP)
        write(7,211) ('-------------',A=FIRSTAGE,EXPANDED_PLUSGROUP+1)
        DO Y=FIRSTYEAR,LASTYEAR ; write(7,200) Y+DISPLAYYEAR-1,(N_AREA(box,a,y),A=FIRSTAGE,EXPANDED_PLUSGROUP) ; END DO
          write(7,202) LASTYEAR+DISPLAYYEAR,(N_AREA(box,a,y),A=FIRSTAGE+1,EXPANDED_PLUSGROUP)
        DO Y=1,99
          IF(Y<=LASTYEAR) THEN ; write(17,200) Y+DISPLAYYEAR-1,(N_AREA(box,a,y),A=FIRSTAGE,EXPANDED_PLUSGROUP)
            ELSE ; WRITE(17,*) ; ENDIF
        END DO
        write(7,211) ('=============',A=FIRSTAGE,EXPANDED_PLUSGROUP+1)
        WRite(7,'(/)')
        IF(NBOX==1 .OR. (LINK_STOCK < 1 .AND. MODEL_TYPE < 2)) GOTO 55
        write(17,*) 'N by stock'
        DO Y=1,99
          IF(Y<=LASTYEAR) THEN ; write(17,200) Y+DISPLAYYEAR-1,(N_STOCK(box,a,y),A=FIRSTAGE,EXPANDED_PLUSGROUP) ; ELSE ; WRITE(17,*) ; ENDIF
        END DO
        IF(NBOX>1) THEN
            write(7,'(A9,I1,A55,A50)') 'TABLE 2a.',box,' ABUNDANCE AT THE BEGINNING OF THE YEAR [BY STOCK] FOR ',TITLE_FISHERY(box)
          ELSE
            write(7,'(A64,A50)') 'TABLE 2a. ABUNDANCE AT THE BEGINNING OF THE YEAR [BY STOCK] FOR ',TITLE_FISHERY(box)
        ENDIF
        write(7,211) ('=============',A=FIRSTAGE,EXPANDED_PLUSGROUP+1)
        write(7,201) (A,A=FIRSTAGE,EXPANDED_PLUSGROUP)
        write(7,211) ('-------------',A=FIRSTAGE,EXPANDED_PLUSGROUP+1)
        DO Y=FIRSTYEAR,LASTYEAR ; write(7,200) Y+DISPLAYYEAR-1,(N_STOCK(box,a,y),A=FIRSTAGE,EXPANDED_PLUSGROUP) ; END DO
           write(7,202) LASTYEAR+DISPLAYYEAR,(N_STOCK(box,a,y),A=FIRSTAGE+1,EXPANDED_PLUSGROUP)
        IF(MODEL_TYPE==1 .AND. PLUSGROUP>1) WRITE(7,*) '     Note: The last age in this table does not represent a plus-group.'
        write(7,211) ('=============',A=FIRSTAGE,EXPANDED_PLUSGROUP+1)
        WRite(7,'(/)')
55      IF(NBOX>1) THEN
            write(7,'(A8,I1,A10,A50)') 'TABLE 3.',box,' CATCH OF ',TITLE_FISHERY(box)
          ELSE
            write(7,'(A18,A50)') 'TABLE 3. CATCH OF ',TITLE_FISHERY(box)
        ENDIF
        write(7,211) ('=============',A=FIRSTAGE,LASTAGE+1)
        write(17,*) 'CATCH'
        write(7,201) (A,A=FIRSTAGE,LASTAGE)
        write(7,211) ('-------------',A=FIRSTAGE,LASTAGE+1)
        DO Y=FIRSTYEAR,LASTYEAR ; write(7,200) Y+DISPLAYYEAR-1,(CATCH_DATA(box,a,y),A=FIRSTAGE,LASTAGE) ; END DO
        write(7,211) ('=============',A=FIRSTAGE,LASTAGE+1)
        DO Y=1,99
          IF(Y<=LASTYEAR) THEN ; write(17,200) Y+DISPLAYYEAR-1,(CATCH_DATA(box,a,y),A=FIRSTAGE,LASTAGE) ; ELSE ; WRITE(17,*) ; ENDIF
        END DO
        write(17,*) 'WEIGHT FOR SSB'
        DO Y=1,99
          IF(Y<=LASTYEAR) THEN ; write(17,200) Y+DISPLAYYEAR-1,(WEIGHT_SSB(box,a,y),A=FIRSTAGE,LASTAGE) ; ELSE ; WRITE(17,*) ; ENDIF
        END DO
        WRite(7,'(/)')
        IF(NBOX>1) THEN
            write(7,'(A8,I1,A45,A50)') 'TABLE 4.',box,' SPAWNING STOCK FECUNDITY AND RECRUITMENT OF ',TITLE_FISHERY(box)
          ELSE
            write(7,'(A53,A50)') 'TABLE 4. SPAWNING STOCK FECUNDITY AND RECRUITMENT OF ',TITLE_FISHERY(box)
        ENDIF
        WRITE(7,*)     '============================================================='
        IF(PDF_STOCKRECRUIT/=0) THEN
            write(7,'(1x,a4,2x,3(a12,1X))')  '    ','  measure of','    recruits','   recruits'
            write(7,'(1x,a4,2x,3(a12,1X))')  'year','   fecundity','    from VPA','  S/R curve'
            WRITE(7,*) '---------------------------------------------'
            DO 70 Y=FIRSTYEAR,LASTYEAR
              IF(Y.GT.FIRSTAGE) THEN ; PRED=STOCK_RECRUIT(box,1)*SSB(box,Y-firstage)/(STOCK_RECRUIT(box,2)+SSB(box,Y-firstage))
                ELSE ; PRED=99999999999999999.0 ; ENDIF
              IF(MODEL_TYPE==1) THEN ; WRITE(7,200) Y+DISPLAYYEAR-1,SSB(BOX,Y),N_AREA(BOX,FIRSTAGE,Y),PRED
                ELSE ; WRITE(7,200) Y+DISPLAYYEAR-1,SSB(BOX,Y),N_STOCK(BOX,FIRSTAGE,Y),PRED ; ENDIF
70          CONTINUE
            write(7,*) '    std. dev. of random component',STOCK_RECRUIT(box,5)
            write(7,*) '    std. dev. of recruitment     ',STOCK_RECRUIT(box,5)/sqrt(1-STOCK_RECRUIT(box,4)**2)
          ELSE  
            write(7,'(1x,a4,2x,3(a12,1X))')  '    ','    spawning','    recruits'
            write(7,'(1x,a4,2x,3(a12,1X))')  'year','     biomass','    from VPA'
            WRITE(7,*) '--------------------------------'
            DO 71 Y=FIRSTYEAR,LASTYEAR
              IF(MODEL_TYPE==1) THEN ; WRITE(7,200) Y+DISPLAYYEAR-1,SSB(BOX,Y),N_AREA(BOX,FIRSTAGE,Y)
                ELSE ; WRITE(7,200) Y+DISPLAYYEAR-1,SSB(BOX,Y),N_STOCK(BOX,FIRSTAGE,Y) ; ENDIF
71          CONTINUE           
        ENDIF
        WRITE(7,*)     '======================================================='
        WRite(7,'(/)')
        
        IF(NBOX>1) THEN
            write(7,'(A8,I1,A24,A50)') 'TABLE 5.',box,' FITS TO INDEX DATA ',TITLE_FISHERY(box)
          ELSE
            write(7,'(A32,A50)') 'TABLE 5. FITS TO INDEX DATA FOR ',TITLE_FISHERY(box)
        ENDIF
        WRITE(7,*) '======================================================================='
        DO 80 g=1,NGEARS(box)
          WRITE(7,*)
          WRITE(7,'(4X,A20)') '--------------------'
          IF(NBOX>1) THEN
              IF(G<10) THEN ; WRITE(7,'(4X,A2,I1,A1,I1,A1,A50)') '5.',BOX,'.',G,' ',TITLE_EFFORT(g,box)
                ELSE ; WRITE(7,'(4X,A2,I1,A1,I2,A1,A50)') '5.',BOX,'.',G,' ',TITLE_EFFORT(g,box) ; ENDIF
            ELSE
              IF(G<10) THEN ; WRITE(7,'(4X,A2,I1,A1,A50)') '5.',G,' ',TITLE_EFFORT(g,box)
                ELSE ; WRITE(7,'(4X,A2,I2,A1,A50)') '5.',G,' ',TITLE_EFFORT(g,box) ; ENDIF
          ENDIF
          WRITE(7,'(4X,A20)') '--------------------'
          WRITE(7,'(4X,A20)') PDFNAME(PDF_EFFORT(g,box))
          IF(PDF_EFFORT(G,box).LE.0) GOTO 80
          SELECT CASE (BIO_EFFORT(g,box))
              CASE (2) ;     CH2='biomass               '
              CASE (3) ;     CH2='MAX(F) - F            '
              CASE (4) ;     CH2='MAX(Z) - Z            '
              CASE (5) ;     CH2='fishing mortality rate'
              CASE (6) ;     CH2='total mortality rate  '
              CASE (7) ;     CH2='proportion- numbers   '
              CASE (8) ;     CH2='proportion- weight    '
              CASE DEFAULT ; CH2='numbers               '
          END SELECT
          IF(SEASON_EFFORT(g,box)<0) THEN
              WRITE(7,'(4X,A8,A21)') 'average ',CH2
            ELSE
              WRITE(7,'(4X,A5,1X,I2,1X,A7)') 'month',SEASON_EFFORT(g,box),CH2
          ENDIF
          write(7,'(4x,a5,i2,1x,a1,1x,i2)') 'Ages ',age_EFFORT(1,g,box),'-',age_EFFORT(2,g,box)
          write(7,'(4x,a21,F12.2)') 'log-likelihood     = ',LIKE_EFFORT(g,box)
          write(7,'(4x,a21,F12.2)') 'deviance           = ',DEV_EFFORT(g,box)
          write(7,'(4x,a21,F12.2)') 'Chi-sq. discrepancy= ',SUM_EFFORT_DISCREPANCY(g,box)
          write(7,*)
          IF (PDF_EFFORT(g,box)==1) then
             write(7,401) '    ','           ','           ',' Residuals ',' Standard  ','     Q     ','Untransfrmd','Untransfrmd',&
                                 'Chi-square '
             write(7,401) 'Year','  Observed ',' Predicted ',' (Obs-pred)',' Deviation ',' Catchabil.','  Observed ',' Predicted ',&
                                 'Discrepancy'
             write(7,401) '----','-----------','-----------','-----------','-----------','-----------','-----------','-----------',&
                                 '-----------'
           ELSE
             write(7,401) '    ','           ','           ',' Residuals ',' Standard  ','     Q     ','Chi-square '
             write(7,401) 'Year','  Observed ',' Predicted ',' (Obs-pred)',' Deviation ',' Catchabil.','Discrepancy'
             write(7,401) '----','-----------','-----------','-----------','-----------','-----------','-----------'
          ENDIF
          DO Y=1,LASTYEAR
            IF((PDF_EFFORT(g,box)>1 .AND. EFFORT_DATA(g,box,y)>-9) .OR. (PDF_EFFORT(g,box)==1 .AND. EFFORT_DATA(g,box,y)>0)) THEN
               CALL GETVARIANCE(EFFORT_DATA(g,box,y),EFFORT(g,box,y),SIGMA_EFFORT(g,box,y),SCALE_VARIANCE(g,box), &
                                PDF_EFFORT(g,box),SIG2)
               IF(PDF_EFFORT(g,box)==1) THEN
                 OBS=DLOG(EFFORT_DATA(g,box,y))
                 IF(EFFORT(g,box,y)<=0) THEN ; PRED=-999999999999999999.0 ; ELSE ; PRED=DLOG(EFFORT(g,box,y)) ; ENDIF
                 WRITE(7,400) y+DISPLAYYEAR-1,OBS,PRED,OBS-PRED,dSQRT(SIG2),Q_EFFORT(g,box,y)*INDEX_MEAN(g,box), &
                              EFFORT_DATA(g,box,y)*INDEX_MEAN(g,box),EFFORT(g,box,y)*INDEX_MEAN(g,box),EFFORT_DISCREPANCY(g,box,y)
                ELSE
                 OBS=EFFORT_DATA(g,box,y) ; PRED=EFFORT(g,box,y)
                 WRITE(7,400) y+DISPLAYYEAR-1,OBS,PRED,OBS-PRED,dSQRT(SIG2),Q_EFFORT(g,box,y)*INDEX_MEAN(g,box), &
                              EFFORT_DISCREPANCY(g,box,y)
               ENDIF
            ENDIF
          END DO
          WRITE(7,*)
          write(7,'(4x,a34)') 'Selectivities by age              '
          write(7,501) 'Year',(a,a=AGE_EFFORT(1,g,box),AGE_EFFORT(2,g,box))
          write(7,499) '----',('-----',a=AGE_EFFORT(1,g,box),AGE_EFFORT(2,g,box))
          DO 82 YEAR=FIRSTYEAR,LASTYEAR
           IF((PDF_EFFORT(g,box)>1 .AND. EFFORT_DATA(g,box,year)>-9) .OR. (PDF_EFFORT(g,box)==1 .AND. EFFORT_DATA(g,box,year)>0)) &
             write(7,500) YEAR+DISPLAYYEAR-1,(SEL_EFFORT(g,box,a,year),a=AGE_EFFORT(1,g,box),AGE_EFFORT(2,g,box))
82        CONTINUE
          write(7,*)
80      CONTINUE        
        WRITE(7,*) '======================================================================='

        write(17,*) 'Indices (paired columns with observed and predicted values)'
        WRITE(17,251) '# ',(G,G,G=1,NGEARS(box))
        DO Y=1,98
          IF(Y<=LASTYEAR) THEN
            DO G=1,NGEARS(box)
              IF(EFFORT_DATA(g,box,y)<0.OR.(PDF_EFFORT(g,box)==1.and.EFFORT_DATA(g,box,y)<=0.0)) THEN
                  mOBS(g)=-9.0; mPRED(g)=-9.0
                ELSE
                  IF(PDF_EFFORT(g,box)==1) THEN
                      mOBS(g)=DLOG(EFFORT_DATA(g,box,y));
                      IF(EFFORT(g,box,y)>0) THEN; mPRED(g)=DLOG(EFFORT(g,box,y)); ELSE; mPRED(g)=-999; ENDIF
                    ELSE
                      mOBS(g)=EFFORT_DATA(g,box,y) ; mPRED(g)=EFFORT(g,box,y)
                  ENDIF
              ENDIF
            END DO
            write(17,250) Y+DISPLAYYEAR-1,(mOBS(g),mPRED(g),g=1,ngears(box))
           ELSE ; WRITE(17,*) ; ENDIF
        END DO

      
        IF(PDF_TAG>=1) THEN
          WRITE(7,'(/)')
          IF(NBOX>1) THEN
              write(7,'(A8,I1,A40,A50)') 'TABLE 6.',box,' FITS TO TAGS RELEASED IN HOME RANGE OF ',TITLE_FISHERY(box)
            ELSE
              write(7,'(A48,A50)') 'TABLE 6. FITS TO TAGS RELEASED IN HOME RANGE OF ',TITLE_FISHERY(box)
          ENDIF
          WRITE(7,'(4X,A20)') PDFNAME(PDF_TAG)
          WRITE(7,'(4x,a17,F12.2)') 'log-likelihood = ',LIKE_TAG(box)
          WRITE(7,'(4x,a17,F12.2)') 'deviance       = ',DEV_TAG(box)
          WRITE(7,*)
          WRITE(7,'(2x,3(A9,1x),a9,6x,a10,2x,a52)') ' Cohort  ','Release  ',' Number  ','Recovery ','Recoveries','Recoveries from fishery (by year following release)'
          WRITE(7,'(3x,A6,2x,A4,1x,A4,1x,A9,4x,a4,9x,a11,1000(i4,4x))') 'ID No.','Year',' Age','released','Area','at end exp.',(y,y=firstyear,lastyear)
          WRITE(7,*) '======================================================================='
          DO J=1,N_tag_cohorts
            IF(TC(j)%Br==box .AND. TC(j)%N>0) THEN
              DO B_REC=1,NBOX
                WRITE(7,'(4x,i4,3x,i4,1x,i4,2x,f6.0,4x,i4,6x,a4,2x,F7.2,2x,900(F7.2,1x))') j,TC(j)%Yr+DISPLAYYEAR-1,TC(j)%Ar,&
                          TC(j)%N,b_rec,'obsd',RECOVERY_DATA(J,b_rec),(RECAPTURE_DATA(J,b_rec,y),y=TC(j)%Yr,TC(j)%Ye)
                WRITE(7,'(4x,i4,3x,i4,1x,i4,2x,f6.0,4x,i4,6x,a4,2x,F7.2,2x,900(F7.2,1x))') j,TC(j)%Yr+DISPLAYYEAR-1,TC(j)%Ar,&
                          TC(j)%N,b_rec,'pred',RECOVERIES(J,b_rec),(RECAPTURES(J,b_rec,y),    y=TC(j)%Yr,TC(j)%Ye)
              END DO ! B_REC
            ENDIF
          END DO ! tag cohorts
          WRITE(7,*) '======================================================================='

        ENDIF
        WRITE(7,'(/)')

30    CONTINUE
      write(7,*) 'TOTAL NUMBER OF FUNCTION EVALUATIONS = ',ICALL

300   format(1x,i4,2x,50(f6.3,1X)) ; 301   format(5x,50(I6,1X)) ; 311   format(50(A7))
200   format(1x,i4,2x,50(f12.0,1X)) ; 201   format(3x,50(i12,1x)) ; 202   format(1x,i4,15x,50(f12.0,1X))
211   format(50(A13))
250   format(1x,i4,2x,1000(E10.3,1X))
251   format(7x,a2,1000(I6,5X))
400   format(4x,i4,4(2x,f11.3),2X,E11.3,3(2X,F11.3))
401   format(4x,a4,10(2x,a11))
499   format(4x,a4,3x,50(a6,1x)) ; 500   format(4x,i4,3x,50(f6.3,1x)) ; 501   format(4x,a4,3x,50(i4,3x))

1000  FORMAT(A27,F12.2) ; 1001 FORMAT(A27,D12.5)
      
      CLOSE(7) ; CLOSE(17)
      RETURN ; END

!------------------------------------------------------------------------------
      SUBROUTINE GET_STOCK
! Computes the abundance of each stock under model 1
!------------------------------------------------------------------------------
      USE STATISTICS ; USE PARAMETERS ;
      IMPLICIT NONE
      INTEGER o,b
      REAL (KIND=8) :: z1,z2,S1,S2

      N_temp=0
!     {the recruitment to each stock is set to that of the corresponding area}          
      DO b=1,2
        DO YEAR=FIRSTYEAR,LASTYEAR
          N_temp(b,b,firstage,year)=N_AREA(b,firstage,year)
      END DO ; END DO

      DO YEAR = FIRSTYEAR,LASTYEAR
        DO 210 b=1,2
         IF(b.EQ.1) THEN ; o=2 ; ELSE ; o=1 ; ENDIF
         DO 220 AGE=FIRSTAGE,LASTAGE
!          {abundance after mixing} 
           N_temp(b,b,age+1,year+1)=N_temp(b,b,age,year)*(1-T(b,age,year))+N_temp(b,o,age,year)*T(o,age,year)
           N_temp(b,o,age+1,year+1)=N_temp(b,b,age,year)*T(b,age,year)+N_temp(b,o,age,year)* &
                                             (1-T(o,age,year))
           Z1=M(b,age,year)+F(b,age,year) ; S1=DEXP(-Z1)
           Z2=M(o,age,year)+F(o,age,year) ; S2=DEXP(-Z2)

!          {compute the tag abundance at beginning of next year (but prior to mixing)}
           N_temp(b,b,age+1,year+1)=N_temp(b,b,age+1,year+1)*S1
           N_temp(b,o,age+1,year+1)=N_temp(b,o,age+1,year+1)*S2
220      CONTINUE
         IF(PLUSGROUP>1) THEN  
!         {accumulate older fish in the last age group}          
          N_temp(b,b,lastage,year+1)=N_temp(b,b,lastage,year+1)+N_temp(b,b,lastage+1,year+1)
          N_temp(b,o,lastage,year+1)=N_temp(b,o,lastage,year+1)+N_temp(b,o,lastage+1,year+1)
         ENDIF
210     CONTINUE
      END DO ! year
      DO YEAR=FIRSTYEAR,LASTYEAR+1
        DO b=1,2
         IF(b.EQ.1) THEN ; o=2 ; ELSE ; o=1 ; ENDIF
         DO AGE=FIRSTAGE,LASTAGE
           N_STOCK(b,age,year)=N_temp(b,b,age,year)+N_temp(b,o,age,year)
      END DO ; END DO ; END DO
      RETURN ; END

!----------------------------------------------------------------------------
      SUBROUTINE BOOTSTRAP
!----------------------------------------------------------------------------
      USE DATUM ; USE STATISTICS ; USE PARAMETERS 
      IMPLICIT NONE
      INTEGER NE(GS,BS),I,II,j,A,G,Y,B,OTH,IRECL,IRECL2,IRECL4,LOOP,NBOOTS,loopkill,MCHECK,COUNT,NAGES, &
              jmax,imax,JMIN,HMIN,HMAX,FOOTNOTE
      REAL (KIND=8) :: EFFORT_STORE(Gs,Bs,0:Ys),E(Gs,Bs,0:Ys),F_RATIO_STORE(Bs,Ys),CV,RAN1, &
                       PARM_BOOT_COVAR(NPAR,NPAR),PARM_BOOT_avg(NPAR),PARM_BOOT_VAR(NPAR),X(NPAR),  &
                       SEL_TERMINAL_STORE(Bs,0:As,2), M_STORE(Bs,0:As,0:Ys),T_STORE(Bs,0:As,0:Ys),CATCH_DATA_STORE(Bs,0:As,0:Ys), &
                       BOOTERROR,E_s(Gs,Bs,0:Ys), PARM_EST_STORE(NPAR), OBJ,N_STORE(Bs,0:As,0:Ys),F_STORE(Bs,0:As,0:Ys), &
                       STOCK_RECRUIT_STORE(Bs,1:5),scale_variance_store(Gs,Bs),aSV(Gs,Bs),varSV(Gs,Bs), &
                       aN(Bs,0:As,0:Ys),aF(Bs,0:As,0:Ys),aC(Bs,0:As,0:Ys),aveM(Bs,0:As,0:Ys),aT(Bs,0:As,0:Ys),aR(Bs,5), &
                       aI(Gs,Bs,0:Ys),aterm(Bs,0:As),N_correct(Bs,0:As,0:Ys),aF_ratio(Bs,0:Ys), &
                       varN(Bs,0:As,0:Ys),varF(Bs,0:As,0:Ys),varM(Bs,0:As,0:Ys),varT(Bs,0:As,0:Ys),varR(Bs,5), &
                       varterm(Bs,0:As),varF_ratio(Bs,0:Ys),sig2,aQ(Gs,Bs,0:Ys),varQ(Gs,Bs,0:Ys),Q_STORE(Gs,Bs,0:Ys)
      REAL (KIND=4) :: sNAA(Bs,0:As,0:Ys),sFAA(Bs,0:As,0:Ys),sCAA(Bs,0:As,0:Ys),sM(Bs,0:As,0:Ys),sT(Bs,0:As,0:Ys),sR(Bs,Ys+6), &
                       sI(Gs,Bs,0:Ys),sTERM(Bs,0:As)
      EXTERNAL BOOTERROR,OBJECTIVE,RAN1
      CHARACTER (LEN=1) :: CONT
      IF(OPTION(3)==0) RETURN
      COUNT=0 ; aN=0 ; aF=0 ; aC=0 ; aveM=0 ; aT=0 ; aR=0 ; ai=0 ; aterm=0 ; aSV=0 ; N_correct=0 ; aF_ratio=0
      varN=0 ; varF=0 ; varM=0 ; varT=0 ; varterm=0 ; varF_ratio=0 ; varR=0 ; varSV=0  ; parm_boot_covar=0
      NBOOTS=ABS(OPTION(3)) ; NAGES=EXPANDED_PLUSGROUP-FIRSTAGE+1
      IRECL=NAGES*NYEARS*NBOX*4 ; IRECL2=NYEARS*NGEARS(1)*4 ; IRECL4=NBOX*NAGES*4 ; IF(NBOX==2) IRECL2=IRECL2+NYEARS*NGEARS(2)*4
      IF(BOOT_FILE_TYPE>0) THEN
          OPEN(62,FILE='NAA.OUT',ACCESS='DIRECT',STATUS='UNKNOWN',RECL=IRECL)
          OPEN(63,FILE='FAA.OUT',ACCESS='DIRECT',STATUS='UNKNOWN',RECL=IRECL)
          OPEN(64,FILE='CAA.OUT',ACCESS='DIRECT',STATUS='UNKNOWN',RECL=IRECL)
          OPEN(65,FILE='MAA.OUT',ACCESS='DIRECT',STATUS='UNKNOWN',RECL=IRECL)
          IF(NBOX>1) OPEN(66,FILE='TAA.OUT',ACCESS='DIRECT',STATUS='UNKNOWN',RECL=IRECL)
          OPEN(67,FILE='IND.OUT',ACCESS='DIRECT',STATUS='UNKNOWN',RECL=IRECL2)
          OPEN(68,FILE='TERM.OUT',ACCESS='DIRECT',STATUS='UNKNOWN',RECL=IRECL4)
          IF(PDF_STOCKRECRUIT.NE.0) THEN
            IRECL=(5+NYEARS)*NBOX*4
            OPEN(71,FILE='SR.OUT',ACCESS='DIRECT',STATUS='UNKNOWN',RECL=IRECL)
          ENDIF
        ELSE
          OPEN(62,FILE='NAA.OUT',STATUS='UNKNOWN')
          OPEN(63,FILE='FAA.OUT',STATUS='UNKNOWN')
          OPEN(64,FILE='CAA.OUT',STATUS='UNKNOWN')
          OPEN(65,FILE='MAA.OUT',STATUS='UNKNOWN')
          IF(NBOX>1) OPEN(66,FILE='TAA.OUT',STATUS='UNKNOWN')
          OPEN(67,FILE='IND.OUT',STATUS='UNKNOWN')
          OPEN(68,FILE='TERM.OUT',STATUS='UNKNOWN')
          IF(PDF_STOCKRECRUIT.NE.0) THEN
            OPEN(71,FILE='SR.OUT',STATUS='UNKNOWN')
          ENDIF
      ENDIF
      OPEN(69,FILE='BAD.OUT',STATUS='UNKNOWN')
      OPEN(70,FILE='OBJ.OUT',STATUS='UNKNOWN')

! \\\\\\\\\\\\\\\\\\\\\\\\\      
! User interactive sequence
! /////////////////////////
      WRITE(*,'(/,/,/,/,/)')
      WRITE(*,*) '*******************************************'
      IF(OPTION(3)>0) THEN
        WRITE(*,*) 'NONPARAMETRIC BOOTSTRAPS OF ABUNDANCE INDICES'
       ELSE
        WRITE(*,*) 'PARAMETRIC BOOTSTRAPS OF ABUNDANCE INDICES'
      ENDIF
      IF(PDF_CATCH(1)>0) WRITE(*,*) 'PARAMETRIC BOOTSTRAPS OF CATCH AT AGE MATRIX'
      WRITE(*,*) '*******************************************'
      WRITE(*,'(/)')
      IF((PARM_KEY(0)+N_Qs)/N_DATA>0.25) THEN
        WRITE(*,*) 'WARNING: The number of parameters (',PARM_KEY(0)+N_Qs,')'
        WRITE(*,*) '         exceeded 25% of the number of data points (',N_DATA,')'
        wRITE(*,*) '         Model performance is likely questionable and bootstrap'
        wRITE(*,*) '         estimates of variance are likely to be biased low'
        WRITE(12,*) 'WARNING: The number of parameters (',PARM_KEY(0)+N_Qs,')'
        WRITE(12,*) '         exceeded 25% of the number of data points (',N_DATA,')'
        wRITE(12,*) '         Model performance is likely questionable and bootstrap'
        wRITE(12,*) '         estimates of variance are likely to be biased low'
        WRITE(*,*)
        WRITE(*,*) 'CONTINUE (Y/N)? '; READ(*,*) CONT; IF (CONT=='N' .OR. CONT=='n') STOP
        WRITE(*,'(/)')
      ENDIF

! \\\\\\\\\\\\\\\\\\\\\\\\\      
! Preserve parameters of base run
! /////////////////////////

!     Indices of abundance 
      DO B=1,NBOX
       DO G = 1,NGEARS(b)
        IF(PDF_EFFORT(g,b)/=0) THEN
          NE(g,b)=0
          DO y = FirstYear,LastYear
            IF(EFFORT_DATA(g,b,y)>0 .or. (EFFORT_DATA(g,b,y)>-9 .and. PDF_EFFORT(g,b)>1)) THEN ! {remember-can't take log of zero}
              NE(g,b)=NE(g,b)+1 ; EFFORT_STORE(g,b,y)=EFFORT(g,b,y)
              CALL GETVARIANCE(EFFORT_DATA(g,b,y),EFFORT(g,b,y),SIGMA_EFFORT(g,b,y),SCALE_VARIANCE(g,b), &
                                  PDF_EFFORT(g,b),SIG2)
              E_s(g,b,y)=SQRT(SIG2)
              IF(PDF_EFFORT(g,b)==1) THEN ! get lognormally-distributed residuals normalized by their std. error
                  E(g,b,NE(g,b))=(LOG(EFFORT_DATA(g,b,y))-LOG(EFFORT(g,b,y)))/E_s(g,b,y)
                ELSE ! get residuals normalized by their std. dev.
                  E(g,b,NE(g,b))=(EFFORT_DATA(g,b,y)-EFFORT(g,b,y))/E_s(g,b,y)
              ENDIF
            ENDIF
          END DO ! Y
          DO y = 1,NE(g,b)
            IF(NPI(g,b)<NE(g,b)) THEN
                IF(STINE_CORR>0) E(g,b,Y)=E(g,b,Y)/SQRT(1.0-dble(NPI(g,b))/dble(NE(g,b)))
              ELSE
                WRITE(*,*) 'ERROR: There are more parameters for index,',g,' (',npi(g,b),')'
                WRITE(*,*) '       than data points (',ne(g,b),')'
                WRITE(12,*) 'ERROR: There are more parameters for index,',g,' (',npi(g,b),')'
                WRITE(12,*) '       than data points (',ne(g,b),')' ; PAUSE; STOP
            ENDIF
          END DO ! Y
        ENDIF
       END DO ! G
      END DO ! B
!     Remaining data and parameters
      DO B=1,NBOX
        IF(PDF_STOCKRECRUIT.NE.0) THEN ; DO Y=1,5 ; STOCK_RECRUIT_STORE(b,y)=STOCK_RECRUIT(b,y) ; END DO ; END IF
        DO G=1,NGEARS(b) ; SCALE_VARIANCE_STORE(g,b)=SCALE_VARIANCE(g,b) ; END DO
        DO y = FirstYear,LastYear
          F_RATIO_STORE(b,y)=F_RATIO(b,y)
          DO G=1,NGEARS(b) ; Q_STORE(g,b,y)=Q_EFFORT(g,b,y) ; END DO ! G
          DO A = Firstage,Lastage
            CATCH_DATA_STORE(B,A,Y)=CATCH_DATA(B,A,Y)
            T_STORE(b,a,y)=T(b,a,y) ; M_STORE(b,a,y)=M(b,a,y) ; F_STORE(b,a,y)=F(b,a,y)
            IF(MODEL_TYPE==1) THEN ; N_STORE(b,a,y)=N_AREA(b,a,y) ; ELSE ; N_STORE(b,a,y)=N_STOCK(b,a,y) ; ENDIF
            IF(Y==lastyear) THEN ; SEL_TERMINAL_STORE(b,a,1)=SEL_TERMINAL(b,a,1) ; SEL_TERMINAL_STORE(b,a,2)=SEL_TERMINAL(b,a,2)
              IF(MODEL_TYPE==1) THEN ; N_STORE(b,a,lastyear+1)=N_AREA(b,a,lastyear+1) 
                ELSE ; N_STORE(b,a,lastyear+1)=N_STOCK(b,a,lastyear+1) ; ENDIF 
            ENDIF
          END DO ! A
        END DO ! Y
      END DO ! B

      DO I=1,PARM_KEY(0)+1 ; PARM_EST_STORE(I)=PARM_EST(I) ; END DO

! \\\\\\\\\\\\\\\\\\\\\\\\\      
! begin bootstrap looping
! /////////////////////////

      WRITE(*,*) 'ORIG. OBJ.',PARM_EST(PARM_KEY(0)+1)
      DO 100 LOOP=1,NBOOTS+1
        LOOPKILL=0 ; MCHECK=0
        IF(LOOP==1) GOTO 53 ! STORE MLE RESULTS IN FIRST RECORD
! {get new indices}
        DO 110 b=1,NBOX
          DO 111 y = FirstYear,LastYear
            DO 111 G=1,NGEARS(b)
              IF(PDF_EFFORT(g,b)==1 .AND. EFFORT_DATA_STORE(g,b,y)>0) THEN
                  IF(OPTION(3)>0) THEN 
                      EFFORT_DATA(g,b,y)=EFFORT_STORE(g,b,y)*EXP(BOOTERROR(g,b,y,E,E_s,NE,SEED,GS,BS,YS))
                    ELSE
!                      write(*,*) g,b,y,effort_store(g,b,y),effort_data_store(g,b,y)
                      CALL GENDATA(EFFORT_STORE(g,b,y),E_S(g,b,y),PDF_EFFORT(g,b),EFFORT_DATA(g,b,y),SEED)
                  ENDIF     
                ELSEIF(PDF_EFFORT(g,b)>1 .AND. EFFORT_DATA_STORE(g,b,y)>-9) THEN
                  IF(OPTION(3)>0) THEN 
                      EFFORT_DATA(g,b,y)=EFFORT_STORE(g,b,y)+BOOTERROR(g,b,y,E,E_s,NE,SEED,GS,BS,YS)
                    ELSE
                      CALL GENDATA(EFFORT_STORE(g,b,y),E_S(g,b,y),PDF_EFFORT(g,b),EFFORT_DATA(g,b,y),SEED)
                  ENDIF
               ENDIF   
111       CONTINUE
110     CONTINUE
! {get new catches}
        DO A=FIRSTAGE,LASTAGE
          DO Y=FIRSTYEAR,LASTYEAR
            DO B=1,NBOX
              IF(PDF_CATCH(b)==1 .AND. CATCH_DATA_STORE(b,a,y)<=0) CATCH_DATA_STORE(b,a,y)=1.0
              CALL GENDATA(CATCH_DATA_STORE(b,a,y),SIGMA_CATCH(b),PDF_CATCH(b),CATCH_DATA(b,a,y),SEED)
        END DO ; END DO ; END DO
        
! {search for best parameter estimates}
        ICALL=0
        DO 50 I=1,PARM_KEY(0)
50        PARM_EST(I)=PARM_EST_STORE(I)
        CALL ESTIMATE(0)
        CALL OBJECTIVE(PARM_EST,OBJ,PARM_KEY(0))
        WRITE(*,'(1X,A10,I4,A23,F14.3)') 'BOOTSTRAP ',LOOP-1,', OBJECTIVE FUNCTION = ',OBJ

! {write REC, MIX, M, C, N and F matrices to 4-bit binary files}
53      DO 951 B=1,NBOX
         IF(B==1) OTH=2 ; IF(B==2) OTH=1
         DO 951 A=FIRSTAGE,EXPANDED_PLUSGROUP
          AGE=A; IF(AGE>PLUSGROUP) AGE=PLUSGROUP
          IF(MODEL_TYPE==1) THEN
              N_AREA(b,a+1,lastyear+1)=(N_AREA(b,a,lastyear)*(1-T(b,age,lastyear))+N_AREA(oth,a,lastyear)*T(oth,age,lastyear))* &
                                  EXP(-F(b,age,lastyear)-M(b,age,lastyear))
              IF(A==EXPANDED_PLUSGROUP) N_AREA(b,a,lastyear+1)=N_AREA(b,a,lastyear+1)+N_AREA(b,a+1,lastyear+1)
            ELSE
              N_STOCK(b,a+1,lastyear+1)=N_STOCK(b,a,lastyear)*((1-T(b,age,lastyear))*EXP(-F(b,age,lastyear)-M(b,age,lastyear))+ &
                                                     T(b,age,lastyear)*EXP(-F(oth,age,lastyear)-M(oth,age,lastyear)))
              IF(A==EXPANDED_PLUSGROUP) N_STOCK(b,a,lastyear+1)=N_STOCK(b,a,lastyear+1)+N_STOCK(b,a+1,lastyear+1)
          ENDIF
951     CONTINUE
        DO 210 b=1,NBOX
          IF(PDF_STOCKRECRUIT.NE.0) THEN
            sR(b,1)=1 ; DO G=1,2 ; sR(b,g+1)=SNGL(STOCK_RECRUIT(b,g)) ; END DO
            sR(b,4)=SNGL(STOCK_RECRUIT(b,5)) ; sR(b,5)=SNGL(STOCK_RECRUIT(b,4))
          ENDIF
          DO 200 y = FirstYear,LastYear
            sR(b,5+Y)=SNGL(ETA_SR(b,y))
            DO 201 g=1,NGEARS(b)
201           sI(g,b,y)=SNGL(EFFORT_DATA(g,b,y)) 
            DO 200 a = Firstage,EXPANDED_PLUSGROUP
              AGE=a; IF(AGE>PLUSGROUP) AGE=PLUSGROUP
              sT(b,a,y)=SNGL(T(b,age,y)) ; sM(b,a,y)=SNGL(M(b,age,y)) ; sCAA(b,a,y)=SNGL(CATCH_DATA(b,a,y)) ; sFAA(b,a,y)=SNGL(F(b,age,y))
              IF(MODEL_TYPE==1) THEN ; sNAA(b,a,y)=SNGL(N_AREA(b,a,y)) ; ELSE ; sNAA(b,a,y)=SNGL(N_STOCK(b,a,y)) ; ENDIF
              IF(sNAA(b,a,y)<0) LOOPKILL=LOOPKILL+1
              IF(Y==lastyear) then ; sTERM(b,a)=SNGL(SEL_TERMINAL(b,age,1)) ;
                IF(MODEL_TYPE==1) THEN ; sNAA(b,a,lastyear+1)=SNGL(N_AREA(b,a,lastyear+1)) 
                  ELSE ; sNAA(b,a,lastyear+1)=SNGL(N_STOCK(b,a,lastyear+1)) ; ENDIF
              ENDIF
200         CONTINUE  
210     CONTINUE        
        IF(BOOT_FILE_TYPE>0) THEN
            WRITE(62,REC=LOOP) (((sNAA(b,a,y),Y=firstyear,lastyear),A=firstage,expanded_plusgroup),b=1,NBOX)
            WRITE(63,REC=LOOP) (((sFAA(b,a,y),Y=firstyear,lastyear),A=firstage,expanded_plusgroup),b=1,NBOX)
            WRITE(64,REC=LOOP) (((sCAA(b,a,y),Y=firstyear,lastyear),A=firstage,expanded_plusgroup),b=1,NBOX)
            WRITE(65,REC=LOOP) (((sM(b,a,y),Y=firstyear,lastyear),A=firstage,expanded_plusgroup),b=1,NBOX)
            IF(nbox>1) WRITE(66,REC=LOOP) (((sT(b,a,y),Y=firstyear,lastyear),A=firstage,expanded_plusgroup),b=1,NBOX)
            WRITE(67,REC=LOOP) ((sI(g,1,y),Y=firstyear,lastyear),g=1,NGEARS(1)),((sI(g,2,y),Y=firstyear,lastyear),g=1,NGEARS(2))
            WRITE(68,REC=LOOP) ((sTerm(b,a),A=firstage,expanded_plusgroup),b=1,NBOX)
            IF(PDF_STOCKRECRUIT.NE.0) THEN
              IRECL=(6+NYEARS)*NBOX*4
              WRITE(71,REC=LOOP) ((sR(b,Y),Y=1,5+lastyear),b=1,NBOX)
            ENDIF
          ELSE
            DO b=1,NBOX;
              DO A=firstage,expanded_plusgroup
                WRITE(62,*) loop,b,a,(sNAA(b,a,y),Y=firstyear,lastyear)
                WRITE(63,*) loop,b,a,(sFAA(b,a,y),Y=firstyear,lastyear)
                WRITE(64,*) loop,b,a,(sCAA(b,a,y),Y=firstyear,lastyear)
                WRITE(65,*) loop,b,a,(sM(b,a,y),Y=firstyear,lastyear)
                IF(nbox>1) WRITE(66,*) loop,b,a,(sT(b,a,y),Y=firstyear,lastyear)
                IF(A==firstage .AND. PDF_STOCKRECRUIT.NE.0) THEN
                  WRITE(71,*) loop,b,(sR(b,Y),Y=1,5+lastyear)
                ENDIF
              END DO; ! A loop
              DO g=1,NGEARS(b); WRITE(67,*) loop,b,g,(sI(g,b,y),Y=firstyear,lastyear) ; END DO;
              WRITE(68,*) loop,b,(sTerm(b,a),A=firstage,lastage)
            END DO;
        ENDIF

! Don't use MLE as though it were a bootstrap for bias corrections
        IF(LOOP==1) GOTO 100

! {report bad records to ASC-II file 'BAD.OUT and }
        DO 51 I=1,PARM_KEY(0)
          IF(PARM_EST(I)<=PARM_SPECS(PARM_KEY(I),1) .OR. PARM_EST(I)>=PARM_SPECS(PARM_KEY(I),3)) MCHECK=MCHECK+1
51      CONTINUE
        IF(MCHECK>0) then
            write(69,*) LOOP,MCHECK,' parameters hit bounds'
          ELSE if(loopkill>0) then
            write(69,*) LOOP,LOOPKILL,' negative stock sizes were predicted'
          ELSE IF(OBJ.GT.PARM_EST_STORE(PARM_KEY(0)+1) .AND. (ABS(OBJ)/ABS(PARM_EST_STORE(PARM_KEY(0)+1)))>10) THEN ;
            write(69,*) LOOP,LOOPKILL,' very poor fit-- objective function was ',obj
            GOTO 100 ! Ignore really bad fits when computing confidence intervals and bias correction (probably a numerical problem)
        ENDIF
        COUNT=COUNT+1
        write(70,*) LOOP,OBJ ! objective function for all loops except numerical blow outs and MLE

! {determine sums and variances of bootstrap parameters}
        DO 320 b=1,NBOX
          IF(PDF_STOCKRECRUIT.NE.0) THEN
            DO Y=1,5 ; CALL VARIANCE_RECURSION(COUNT,aR(b,y),varR(b,y),STOCK_RECRUIT(b,y)) ; END DO
          ENDIF
          DO G=1,NGEARS(b) ; CALL VARIANCE_RECURSION(COUNT,aSV(g,b),varSV(g,b),SCALE_VARIANCE(g,b)) ; END DO
          DO 320 y = FirstYear,LastYear
            DO 321 g=1,NGEARS(b)
              CALL VARIANCE_RECURSION(COUNT,aQ(g,b,y),varQ(g,b,y),Q_EFFORT(g,b,y))
              aI(g,b,y)=aI(g,b,y)+dble(sI(g,b,y))
321         CONTINUE
            DO 322 a = Firstage,Lastage
              if(y==lastyear) then
                CALL VARIANCE_RECURSION(COUNT,aTERM(b,a),varTERM(b,a),dble(sTERM(b,a)))
                CALL VARIANCE_RECURSION(COUNT,aN(b,a,lastyear+1),varN(b,a,lastyear+1),dble(sNAA(b,a,lastyear+1)))
              endif
              CALL VARIANCE_RECURSION(COUNT,aT(b,a,y),varT(b,a,y),dble(sT(b,a,y)))
              CALL VARIANCE_RECURSION(COUNT,aveM(b,a,y),varM(b,a,y),dble(sM(b,a,y)))
              CALL VARIANCE_RECURSION(COUNT,aF(b,a,y),varF(b,a,y),dble(sFaa(b,a,y)))
              CALL VARIANCE_RECURSION(COUNT,aN(b,a,y),varN(b,a,y),dble(sNaa(b,a,y)))
              aC(b,a,y)=aC(b,a,y)+dble(sCAA(b,a,y))
322         CONTINUE
            CALL VARIANCE_RECURSION(COUNT,aF_RATIO(b,y),varF_RATIO(b,y),dble(sFaa(b,lastage,y)/sFaa(b,lastage-1,y)))
320     CONTINUE
        DO I=1,PARM_KEY(0)
          CALL VARIANCE_RECURSION(COUNT,PARM_BOOT_AVG(i),PARM_BOOT_VAR(i),PARM_EST(i))
          DO J=1,I
            PARM_BOOT_COVAR(I,J)=PARM_BOOT_COVAR(I,J)+PARM_EST(I)*PARM_EST(j)
          END DO
        END DO

100   CONTINUE  ! End bootstrap loop
      WRITE(69,*) 'Number of good bootstrap loops (not counting the MLE) = ',COUNT

! {determine means from sums of bootstrap parameters (computed above)}
      DO 350 b=1,NBOX
        IF(PDF_STOCKRECRUIT.NE.0) THEN ; DO Y=1,5 ; aR(b,y)=aR(b,y)/COUNT ; END DO ; ENDIF
        DO g=1,NGEARS(b) ; aSV(g,b)=aSV(g,b)/COUNT ; END DO
        DO 350 y = FirstYear,LastYear
          DO g=1,NGEARS(b) ; aI(g,b,y)=aI(g,b,y)/COUNT ; aQ(g,b,y)=aQ(g,b,y)/COUNT ; END DO
          DO 352 a = Firstage,Lastage
            if(y==lastyear) then
              aTERM(b,A)=aTERM(b,a)/COUNT ; aN(b,a,LASTYEAR+1)=aN(b,a,lastyear+1)/COUNT 
            endif
            aT(b,A,y)=aT(b,a,y)/COUNT ; aF(b,a,Y)=aF(b,a,y)/COUNT ; aveM(b,a,Y)=aveM(b,a,y)/COUNT ; aC(b,a,Y)=aC(b,a,y)/COUNT 
            aN(b,a,Y)=aN(b,a,y)/COUNT 
352       CONTINUE                
          aF_RATIO(b,y)=aF_RATIO(b,y)/COUNT
350   CONTINUE        
      PARM_BOOT_AVG=PARM_BOOT_AVG/DBLE(COUNT)
      DO i=1,PARM_KEY(0) ; DO J=1,i
        PARM_BOOT_COVAR(i,J)=(PARM_BOOT_COVAR(i,J)/DBLE(COUNT)-PARM_BOOT_AVG(i)*PARM_BOOT_AVG(j))
      END DO ; END DO

      DO i=1,7 ; CLOSE(61+i) ; END DO ; CLOSE(69)
      
! {search for best parameter estimates using 'Average' data}
      DO 430 b=1,NBOX
        DO G=1,NGEARS(b) ; SCALE_VARIANCE(g,b)=aSV(g,b) ; END DO
        DO 431 y = FirstYear,LastYear
          F_RATIO(b,y)=aF_RATIO(b,y)
          DO 432 g=1,NGEARS(b)
432         EFFORT_DATA(g,b,y)=aI(g,b,y)
          DO 431 a = Firstage,Lastage
            T(b,a,y)=aT(b,a,y) ; M(b,a,y)=aveM(b,a,y) ; Catch_Data(b,a,y)=aC(b,a,y) ; SEL_TERMINAL(b,a,1)=aTERM(b,a)
431     CONTINUE
430   CONTINUE
      DO 433 I=1,PARM_KEY(0)
433     PARM_EST(I)=PARM_EST_STORE(I)
      CALL ESTIMATE(0)
      CALL OBJECTIVE(PARM_EST,OBJ,PARM_KEY(0))
      write(*,*) 'obj of ave',obj
      DO 451 B=1,NBOX
       IF(B==1) THEN ; OTH=2 ; ELSEIF(B==2) THEN ; OTH=1 ; ENDIF
       DO 451 A=FIRSTAGE,LASTAGE
         IF(MODEL_TYPE==1) THEN
            N_AREA(B,a+1,lastyear+1)=(N_AREA(B,a,lastyear)*(1-T(B,a,lastyear))+N_AREA(oth,a,lastyear)*T(oth,a,lastyear))* &
                                EXP(-F(B,a,lastyear)-M(B,a,lastyear))
            IF(A==LASTAGE) N_AREA(B,a,lastyear+1)=N_AREA(B,a,lastyear+1)+N_AREA(B,a+1,lastyear+1)
            N_correct(b,a,lastyear+1)=aN(b,a,lastyear+1)-N_AREA(b,a,lastyear+1)
          ELSE
            N_STOCK(B,a+1,lastyear+1)=N_STOCK(B,a,lastyear)*((1-T(B,a,lastyear))*EXP(-F(B,a,lastyear)-M(B,a,lastyear))+ &
                                                   T(B,a,lastyear)*EXP(-F(oth,a,lastyear)-M(oth,a,lastyear)))
            IF(A==LASTAGE) N_STOCK(B,a,lastyear+1)=N_STOCK(B,a,lastyear+1)+N_STOCK(B,a+1,lastyear+1)
            N_correct(b,a,lastyear+1)=aN(b,a,lastyear+1)-N_STOCK(b,a,lastyear+1)
         ENDIF
451    CONTINUE

! {Compute bias with Efron's (1992) estimator, bias = average of P(bootstrap) - P(average bootstrap)}
! {except for scale_variance, which will be near zero since residuals of average case will be small by definition, 
! {so the formula used for scale_variance is the standard one, bias = average of P(bootstrap)- MLE }
      DO 440 B=1,NBOX
        DO G=1,NGEARS(b) ; SCALE_VARIANCE(g,b)=aSV(g,b)-SCALE_VARIANCE_store(g,b); END DO
        DO 440 y = FirstYear,LastYear
          F_RATIO(b,y)=aF_RATIO(b,y)-F_RATIO(b,y)
         DO 441 g=1,ngears(b)
441        Q_effort(g,b,y)=aQ(g,b,y)-Q_effort(g,b,y)
         DO 440 A = Firstage,Lastage
           IF(Y==lastyear) SEL_TERMINAL(b,a,1)=aTERM(b,a)-SEL_TERMINAL(b,a,1)
           T(b,a,y)=aT(b,a,y)-T(b,a,y) ; M(b,a,y)=aveM(b,a,y)-M(b,a,y)
           F(b,a,y)=aF(b,a,y)-F(b,a,y)
           IF(MODEL_TYPE==1) THEN ; N_correct(b,a,y)=aN(b,a,y)-N_AREA(b,a,y)
            ELSE ; N_correct(b,a,y)=aN(b,a,y)-N_STOCK(b,a,y)
           ENDIF
440   CONTINUE

      OPEN(51,FILE='BOOTSTRP.OUT',STATUS='UNKNOWN')
      WRITE(51,'(16X,A40)') '****************************************'
      WRITE(51,'(16X,A40)') '                VPA-2BOX                '
      WRITE(51,'(16X,A40)') 'BOOTSTRAP ESTIMATES OF BIAS AND VARIANCE'
      WRITE(51,'(16X,A40)') '****************************************'
      write(51,'(/)')
      WRITE(51,'(A50)') INFILE(0)     ! Title of run
      A=DATE_VALUES(6)/10 ; G=DATE_VALUES(6)-A*10 ; Y=DATE_VALUES(1) ;
      WRITE(51,'(I2,A1,I1,I1,A2,I2,A1,A9,A1,I4)') DATE_VALUES(5),':',A,G,', ',DATE_VALUES(3),' ',MONTH_NAME(DATE_VALUES(2)),' ',Y
      write(51,'(/)')
      WRITE(51,*) '*************************************************'
      WRITE(51,*) 'BIAS AND STANDARD ERROR ESTIMATES'
      WRITE(51,*) '*************************************************'
      write(51,'(/)')
      DO 450 B=1,NBOX
        write(51,'(a6,i1,A38,A50)') 'TABLE ',b,'A. BIAS OF FISHING MORTALITY RATE FOR ',TITLE_FISHERY(b)
        write(51,1311) ('=======',A=FIRSTAGE,LASTAGE+1)
        write(51,1301) (A,A=FIRSTAGE,LASTAGE)
        write(51,1311) ('-------',A=FIRSTAGE,LASTAGE+1)
        DO 460 Y=FIRSTYEAR,LASTYEAR
460       write(51,1300) Y+DISPLAYYEAR-1,(F(b,a,y),A=FIRSTAGE,LASTAGE)
        write(51,1311) ('=======',A=FIRSTAGE,LASTAGE+1)
        write(51,'(/)')
        write(51,'(a6,i1,A48,A50)') 'TABLE ',b,'B. STANDARD ERROR OF FISHING MORTALITY RATE FOR ',TITLE_FISHERY(b)
        write(51,1211) ('===========',A=FIRSTAGE,LASTAGE),'======     '
        write(51,1203) (A,A=FIRSTAGE,LASTAGE)
        write(51,1211) ('-----------',A=FIRSTAGE,LASTAGE),'------     '
        DO 461 Y=FIRSTYEAR,LASTYEAR
461       write(51,1200) Y+DISPLAYYEAR-1,(sqrt(varF(b,a,y)),A=FIRSTAGE,LASTAGE)
        write(51,1211) ('===========',A=FIRSTAGE,LASTAGE),'======     '
        write(51,'(/)')
        IF(MODEL_TYPE==1) THEN
          write(51,'(a6,i1,A34,A50)') 'TABLE ',b,'C. BIAS OF ABUNDANCE (BY AREA) OF ',TITLE_FISHERY(b)
         ELSE
          write(51,'(a6,i1,A35,A50)') 'TABLE ',b,'C. BIAS OF ABUNDANCE (BY STOCK) OF ',TITLE_FISHERY(b)
        ENDIF
        write(51,1211) ('===========',A=FIRSTAGE,LASTAGE),'======     '
        write(51,1201) (A,A=FIRSTAGE,LASTAGE)
        write(51,1211) ('-----------',A=FIRSTAGE,LASTAGE),'------     '
        DO 470 Y=FIRSTYEAR,LASTYEAR
470         write(51,1200) Y+DISPLAYYEAR-1,(N_correct(b,a,y),A=FIRSTAGE,LASTAGE)
            write(51,1210) Y+DISPLAYYEAR-1,(N_correct(b,a,y),A=FIRSTAGE+1,LASTAGE)
        write(51,1211) ('===========',A=FIRSTAGE,LASTAGE),'======     '
        write(51,'(/)')
        IF(MODEL_TYPE==1) THEN
          write(51,'(a6,i1,A44,A50)') 'TABLE ',b,'D. STANDARD ERROR OF ABUNDANCE (BY AREA) OF ',TITLE_FISHERY(b)
         ELSE
          write(51,'(a6,i1,A45,A50)') 'TABLE ',b,'D. STANDARD ERROR OF ABUNDANCE (BY STOCK) OF ',TITLE_FISHERY(b)
        ENDIF
        write(51,1211) ('===========',A=FIRSTAGE,LASTAGE),'======     '
        write(51,1203) (A,A=FIRSTAGE,LASTAGE)
        write(51,1211) ('-----------',A=FIRSTAGE,LASTAGE),'------     '
        DO 471 Y=FIRSTYEAR,LASTYEAR
471         write(51,1200) Y+DISPLAYYEAR-1,(sqrt(varN(b,a,y)),A=FIRSTAGE,LASTAGE)
            write(51,1210) Y+DISPLAYYEAR-1,(sqrt(varN(b,a,y)),A=FIRSTAGE+1,LASTAGE)
        write(51,1211) ('===========',A=FIRSTAGE,LASTAGE),'======     '
        WRITE(51,'(/)')
        write(51,'(a6,i1,A27,A50)') 'TABLE ',b,'E. PARAMETER ESTIMATES FOR ',TITLE_FISHERY(b)
        write(51,1211) ('===========',A=1,6)
        WRITE(51,*)
        WRITE(51,*) 'TERMINAL AGE STRUCTURE OF POPULATION ABUNDANCE'
        WRITE(51,1204) '       ','          ','Average of','          ','   Std.   ','          '
        WRITE(51,1204) '  Age  ','    MLE   ','bootstraps','   Bias   ','   Error  ','% CV      '
        WRITE(51,1204) '  ---  ','----------','----------','----------','----------','-----     '
        DO 480 A=FIRSTAGE+1,LASTAGE
          IF(aN(b,a,lastyear+1).NE.0.0) THEN ; CV=100*SQRT(varN(b,a,lastyear+1))/aN(b,a,lastyear+1) ; ELSE ; CV=0 ; ENDIF
480       WRITE(51,1205) A,N_store(b,a,Lastyear+1),aN(b,a,lastyear+1),N_correct(b,a,lastyear+1),SQRT(varN(b,a,lastyear+1)),CV
        WRITE(51,*)
        WRITE(51,*) 'TERMINAL AGE STRUCTURE OF FISHING MORTALITY RATE'
        WRITE(51,1204) '       ','          ','Average of','          ','   Std.   ','          '
        WRITE(51,1204) '  Age  ','    MLE   ','bootstraps','   Bias   ','   Error  ','% CV      '
        WRITE(51,1204) '  ---  ','----------','----------','----------','----------','-----     '
        DO 485 A=FIRSTAGE,LASTAGE
          IF(aF(b,a,lastyear).NE.0.0) THEN ; CV=100*SQRT(varF(b,a,lastyear))/aF(b,a,lastyear) ; ELSE ; CV=0 ; ENDIF
485       WRITE(51,1205) A,F_store(b,a,Lastyear),aF(b,a,lastyear),F(b,a,lastyear),SQRT(varF(b,a,lastyear)),CV
        WRITE(51,*)
        WRITE(51,*) 'RATIO OF FISHING MORTALITY RATE ON LAST TWO AGES (F-RATIO)'
        WRITE(51,1204) '       ','          ','Average of','          ','   Std.   ','          '
        WRITE(51,1204) ' Year  ','    MLE   ','bootstraps','   Bias   ','   Error  ','% CV      '
        WRITE(51,1204) ' ----  ','----------','----------','----------','----------','-----     '
        DO 490 Y=FIRSTYEAR,LASTYEAR
          IF(aF_RATIO(b,y).NE.0.0) THEN ; CV=100*SQRT(varF_RATIO(b,y))/aF_RATIO(b,y) ; ELSE ; CV=0 ; ENDIF
          CV=100*SQRT(varF_RATIO(b,y))/aF_RATIO(b,y)
490       WRITE(51,1205) Y+displayyear-1,F_RATIO_STORE(b,y),aF_RATIO(b,y),F_RATIO(b,y),SQRT(varF_RATIO(b,y)),CV
        WRITE(51,*)
        WRITE(51,*) 'NATURAL MORTALITY RATE'
        WRITE(51,1204) '       ','          ','Average of','          ','   Std.   ','          '
        WRITE(51,1204) '  Age  ','    MLE   ','bootstraps','   Bias   ','   Error  ','% CV      '
        WRITE(51,1204) '  ---  ','----------','----------','----------','----------','-----     '
        DO 500 A=FIRSTAGE,LASTAGE
          IF(aveM(b,a,lastyear).NE.0.0) THEN ; CV=100*SQRT(varM(b,a,lastyear))/aveM(b,a,lastyear)  ; ELSE ; CV=0 ; ENDIF
500       WRITE(51,1205) A,M_store(b,a,Lastyear),aveM(b,a,lastyear),M(b,a,lastyear),SQRT(varM(b,a,lastyear)),CV
        WRITE(51,*)
        WRITE(51,*) 'TRANSFER COEFFICIENTS'
        WRITE(51,1204) '       ','          ','Average of','          ','   Std.   ','          '
        WRITE(51,1204) '  Age  ','    MLE   ','bootstraps','   Bias   ','   Error  ','% CV      '
        WRITE(51,1204) '  ---  ','----------','----------','----------','----------','-----     '
        DO 510 A=FIRSTAGE,LASTAGE
          IF(aT(b,a,lastyear).NE.0.0) THEN ; CV=100*SQRT(varT(b,a,lastyear))/aT(b,a,lastyear)  ; ELSE ; CV=0 ; ENDIF
510       WRITE(51,1205) A,T_store(b,a,Lastyear),aT(b,a,lastyear),T(b,a,lastyear),SQRT(varT(b,a,lastyear)),CV
        WRITE(51,*)
        IF(PDF_STOCKRECRUIT.NE.0) THEN
          WRITE(51,*)
          WRITE(51,*) 'STOCK-RECRUITMENT COEFFICIENTS'
          WRITE(51,1204) ' Param.','          ','Average of','          ','   Std.   ','          '
          WRITE(51,1204) '  ID#  ','    MLE   ','bootstraps','   Bias   ','   Error  ','% CV      '
          WRITE(51,1204) ' ----- ','----------','----------','----------','----------','-----     '
          DO a=1,5
            IF(aR(b,a).NE.0.0) THEN ; CV=100*SQRT(varR(b,a))/aR(b,a)  ; ELSE ; CV=0 ; ENDIF
            WRITE(51,1205) A,stock_recruit_store(b,a),aR(b,a),stock_recruit(b,a),SQRT(varR(b,a)),CV
          END DO
          WRITE(51,*)
        ENDIF

        WRITE(51,*)
        WRITE(51,*) 'VARIANCE SCALING PARAMETERS'
        WRITE(51,1204) '       ','          ','Average of','          ','   Std.   ','          '
        WRITE(51,1204) ' Index ','    MLE   ','bootstraps','   Bias   ','   Error  ','% CV      '
        WRITE(51,1204) ' ----- ','----------','----------','----------','----------','-----     '
        DO G=1,ngears(b)
          IF(aSV(g,b).NE.0.0) THEN ; CV=100*SQRT(varSV(g,b))/aSV(g,b)  ; ELSE ; CV=0 ; ENDIF
          WRITE(51,1205) G,SCALE_VARIANCE_store(g,b),aSV(g,b),SCALE_VARIANCE(g,b),SQRT(varSV(g,b)),CV
        END DO
        WRITE(51,*)
        WRITE(51,*) 'CATCHABILITY COEFFICIENTS'
        WRITE(51,1216) '       ','       ','          ','Average of','          ','   Std.   ','          '
        WRITE(51,1216) ' Index ',' Year  ','    MLE   ','bootstraps','   Bias   ','   Error  ','% CV      '
        WRITE(51,1216) ' ----- ',' ----  ','----------','----------','----------','----------','-----     '
        DO 520 G=1,ngears(b)
          DO 520 y=FIRSTYEAR,LASTYEAR
          IF(aQ(g,b,y).NE.0.0) THEN ; CV=100*SQRT(varQ(g,b,y))/aQ(g,b,y)  ; ELSE ; CV=0 ; ENDIF
520       WRITE(51,1206) G,Y+displayyear-1,Q_store(g,b,y),aQ(g,b,y),Q_effort(g,b,y),SQRT(varQ(g,b,y)),CV
        WRITE(51,*)

450    CONTINUE

       WRITE(51,'(/)')
       WRITE(51,*) '*************************************************'
       WRITE(51,*) 'CORRELATION AND COVARIANCE MATRICES OF PARAMETERS'
       WRITE(51,*) '*************************************************'
       WRITE(51,'(/)')
       write(51,'(a6,i1,A21)') 'TABLE ',1+NBOX,'A. CORRELATION MATRIX'
       write(51,'(22(A5))') ('=====',I=1,MIN(PARM_KEY(0),20)+1),'===  '
       IMAX=0 ; JMAX=0
       DO WHILE (IMAX<PARM_KEY(0))
          IMAX=MIN(PARM_KEY(0),20+IMAX) ; JMIN=JMAX+1
          DO II=0,(IMAX-1)/20
            HMIN=20*II+1 ; HMAX=HMIN+19  ; JMAX=JMIN+19 ; IF(JMAX>IMAX) JMAX=IMAX ; FOOTNOTE=0
            DO I=JMIN,JMAX
              DO J=HMIN,I
                IF(PARM_BOOT_COVAR(I,I)<0 .OR. PARM_BOOT_COVAR(J,J)<0) THEN
                    X(J)=9 ; FOOTNOTE=1
                  ELSEIF(PARM_BOOT_COVAR(I,J)==0.0) THEN
                    X(J)=0 ! zero covariances imply zero correlations no matter what
                  ELSEIF(PARM_BOOT_COVAR(I,I)==0 .OR. PARM_BOOT_COVAR(J,J)==0) THEN
                    X(J)=8 ; FOOTNOTE=1 ! nonzero covariance with zero variances make no sense
                  ELSE ; X(J)=PARM_BOOT_COVAR(I,J)/SQRT(PARM_BOOT_COVAR(I,I))/SQRT(PARM_BOOT_COVAR(J,J))
                         IF(ABS(X(J))>1.0 .and. ABS(X(J))<1.1) THEN ; X(J)=1.0
                           ELSEIF(ABS(X(J))>1.1) THEN ; X(J)=999999999999999999.0 ; ENDIF
                ENDIF
              END DO
              IF(I==JMIN) WRITE(51,'(7X,20(I4,1X))') (J,J=HMIN,JMAX)
              WRITE(51,'(1X,I4,A3,20(F4.1,1X))') I,' : ',(X(J),J=HMIN,MIN(HMAX,I))
            END DO
            IF(FOOTNOTE>0) THEN
              WRITE(51,*) '    Values of 9.0 indicate that there were negative variances'
              WRITE(51,*) '    Values of 8.0 indicate that the covariance was nonzero but'
              WRITE(51,*) '      one or both variances were zero'
            ENDIF
            WRITE(51,*)
          END DO
          WRITE(51,*)
       END DO

       write(51,'(a6,i1,A20)') 'TABLE ',1+NBOX,'B. COVARIANCE MATRIX'
       write(51,'(22(A13))') ('=============',I=1,MIN(PARM_KEY(0),20)),'=======      '
       IMAX=0 ; JMAX=0
       DO WHILE (IMAX<PARM_KEY(0))
          IMAX=MIN(PARM_KEY(0),20+IMAX) ; JMIN=JMAX+1
          DO II=0,(IMAX-1)/20
            HMIN=20*II+1 ; HMAX=HMIN+19  ; JMAX=JMIN+19 ; IF(JMAX>IMAX) JMAX=IMAX
            DO I=JMIN,JMAX
              DO J=HMIN,I
                  X(J)=PARM_BOOT_COVAR(I,J)
              END DO
              IF(I==JMIN) WRITE(51,'(7X,20(I8,5X))') (J,J=HMIN,JMAX)
              WRITE(51,'(1X,I4,A3,20(d12.4,1X))') I,' : ',(X(J),J=HMIN,MIN(HMAX,I))
            END DO
            WRITE(51,*)
          END DO
          WRITE(51,*)
       END DO

1300   format(1x,i5,2x,50(f6.3,1X))
1301   format(3x,50(i6,1x)) ; 1311   format(50(a7))
1200   format(1x,i5,1x,50(E10.3,1X)) ; 1210   format(1x,i5,12x,50(E10.3,1X))
1201   format(3x,50(i10,1x)) ; 1211   format(50(a11))
1203   format(3x,50(i10,1x))
1204   format(A7,40(A10,1x),A7)
1205   format(1x,i4,2x,4(E10.3,1X),F5.1)
1206   format(2(1x,i4,2x),4(E10.3,1X),F5.1) ; 1216   format(2(a7),4(a10,1X),a7)
       CLOSE(51) ; RETURN ; END

!----------------------------------------------------------------------------
      REAL (KIND=8) FUNCTION BOOTERROR(ind,iar,Y,E,E_s,NE,ISEED,GS,BS,YS)
!----------------------------------------------------------------------------
      INTEGER GS,BS,YS,IND,IAR,ISEED
      REAL (KIND=8) :: E(GS,BS,0:YS),E_s(Gs,Bs,0:Ys),RAN1,U,RU
      INTEGER NE(GS,BS),I,Y
      EXTERNAL RAN1
      
      U=RAN1(ISEED)
      RU=DBLE(NE(IND,IAR))*U
      DO 10 I=1,NE(IND,IAR)
        IF(RU.LE.FLOAT(I)) THEN
          BOOTERROR=E(IND,IAR,I)*E_s(ind,iar,y)
          RETURN
        ENDIF
10    CONTINUE

      END

!----------------------------------------------------------------------------
      REAL (KIND=8) FUNCTION MFEXP(X)
!----------------------------------------------------------------------------
! approximates intrinsic function EXP() for large arguments that otherwise would cause numeric overflows
! based on function defined in AD Model Builder (Fournier et al., Otter research Lmtd)
      REAL (KIND=8) :: X,XMAX=60

      IF (X<XMAX) THEN
        MFEXP=EXP(X)
      ELSE
        MFEXP=EXP(XMAX)*(1.0+2.0*(X-XMAX))/(1.0+X-XMAX)
      ENDIF
      RETURN
      END



