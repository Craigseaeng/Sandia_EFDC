SUBROUTINE WQSKE1  

  ! *** ORGINALLY CODED BY K.-Y. PARK 
  ! *** WQSKE1 - CE-QUAL-ICM KINETICS WITH UPDATES
  !
  !----------------------------------------------------------------------!  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2011-09           Paul M. Craig    Fixed Reaeration
  ! 2011-07           Paul M. Craig    Rewritten to F90
  ! 2008              SCOTT JAMES      ADDED CARBON DIOXIDE
  ! 2006-01-12        PAUL M. CRAIG    MAJOR REWRITE

  USE GLOBAL  

  IMPLICIT NONE
  
  real::wqv10last

  INTEGER*4 NQ,NW,NS,IZ,IMWQZ,NSTPTMP
  INTEGER*4 K,L,CITER

  REAL WQAVGIO,CNS1,RMULTMP,TIME,RLIGHT1,RLIGHT2
  REAL WQGNC,WQGND,WQGNG,WQGNM,WQGPM,WQF1NM,WQGPC,WQGPD,WQGPG
  REAL WQKESS,XMRM,YMRM,WQTT1
  REAL WQFDI0,WQTTT
  REAL SADWQ,WQGSD,WQTTB,WQISM,WQFDM,WQF2IM
  REAL UMRM,VMRM,WQVEL,WQLVF,WQF4SC,WQKDOC,WQKHP,WQTTS
  REAL WQKHN,WQTTM,TVAL1,TVAL2,TVAL3,TVAL4,TVAL5
  REAL RLNSAT1,RLNSAT2,XNUMER,XDENOM,WQLDF,WQTTC,WQTTD,WQTTG
  REAL WINDREA,WQWREA,WQVREA,WQA1C,WQVA1C,WQR1C,WQA2D
  REAL WQR2D,WQA3G,WQR3G,WQB4,WQA4,WQR4,WQC5,WQA5,WQR5
  REAL WQD6,WQA6C,WQA6D,WQA6G,WQA6,WQA6M,WQR6
  REAL WQE7,WQA7C,WQA7D,WQA7G,WQA7,WQR7
  REAL WQF8,WQA8C,WQA8D,WQA8G,WQA8,WQR8
  REAL WQF9,WQA9C,WQA9D,WQA9G,WQA9,WQR9
  REAL WQA10C,WQA10D,WQA10G,WQR10,WQKKL
  REAL WQI11,WQA11C,WQA11D,WQA11G,WQA11,WQR11
  REAL WQJ12,WQA12C,WQA12D,WQA12G,WQA12,WQR12
  REAL WQF13,WQA13C,WQA13D,WQA13G,WQA13,WQR13
  REAL WQR14,WQF14,WQA14C,WQA14D,WQA14G,WQA14
  REAL WQR15,WQA15C,WQA15D,WQA15G,WQA15,WQB15
  REAL WQM16,WQA16D,WQR16,WQR17,WQR18
  REAL TMP19,TEMFAC,DTWQxH,DTWQxH2,WQA19C,WQA19D,WQA19G
  REAL WQA19,WQA19A,WQSUM,WQRea,WQPOC,WQDOC,WQNH3,WQCOD
  REAL WQT20,WQR21,TIMTMP,WQTAMD
  REAL PPCDO,WQA22, WQA22C, WQA22D, WQA22G, WQCDDOC
  REAL WQCDREA, WQCDSUM
  REAL EXPA0,EXPA1 !VARIABLES FOR LIGHT EXTINCTION

  REAL WQGCO2M,WQGCO2C,WQGCO2G,WQGCO2D		!  CO2 Limitation Consts added by AA
  !REAL WQPHDC, WQPHDD, WQPHDG			!  pH Limitation changes - VJ - These are the actual multiplying factors.
  REAL WQPHKH, WQPHKOH, WQPHHCONC		!  For new pH functional form K_H and K_OH hydration constants. WQPHHCONC - Hydrogen ion conc
  REAL CTime, CPos, CPH                 !  Temporary variables needed for handling negative half-sat concentrations in pH calculations 


  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DELKC  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DZCHP
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::TSSS_ABOVE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::CHLS_ABOVE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::POMS_ABOVE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DZC_ABOVE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::WQISC
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::WQISD
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::WQISG
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::WQI0TOP

  ! ***  1) CHC - cyanobacteria
  ! ***  2) CHD - diatom algae
  ! ***  3) CHG - green algae
  ! ***  4) ROC - refractory particulate organic carbon
  ! ***  5) LOC - labile particulate organic carbon
  ! ***  6) DOC - dissolved organic carbon
  ! ***  7) ROP - refractory particulate organic phosphorus
  ! ***  8) LOP - labile particulate organic phosphorus
  ! ***  9) DOP - dissolved organic phosphorus
  ! *** 10) P4D - total phosphate
  ! *** 11) RON - refractory particulate organic nitrogen 23) macroalgae
  ! *** 12) LON - labile particulate organic nitrogen
  ! *** 13) DON - dissolved organic nitrogen
  ! *** 14) NHX - ammonia nitrogen
  ! *** 15) NOX - nitrate nitrogen
  ! *** 16) SUU - particulate biogenic silica
  ! *** 17) SAA - dissolved available silica
  ! *** 18) COD - chemical oxygen demand
  ! *** 19) DOX - dissolved oxygen
  ! *** 20) TAM - total active metal
  ! *** 21) FCB - fecal coliform bacteria
  ! *** 22) CO2 - dissolved carbon dioxide
  ! *** 23) macroalgae

  ! *** DTWQ - Water quality time step, which is in units of days
  ! *** DTWQO2 = DTWQ*0.5 

  ! *** WQCHL   = Chlorophyll a (ug/l)
  ! *** WQCHLC  = carbon-to-chlorophyll ratio for cyanobacteria (mg C / ug Chl)
  ! *** WQCHLD  = carbon-to-chlorophyll ratio for algae diatoms (mg C / ug Chl)
  ! *** WQCHLG  = carbon-to-chlorophyll ratio for algae greens (mg C / ug Chl)
  ! *** WQKECHL = Light Extinction Coeff for CHLa (1/m per mg/l)
  ! *** WQKETSS = Light Extinction Coeff for TSS (1/m per mg/l)
  ! *** WQKEPOM = Light Extinction Coeff for POM (1/m per mg/l)
  ! *** WQKETOT(L,K) = Total Light Extinction
  ! *** WQDOPG  = Optimal Depth for Growth - Green Algae

  ! *** RNH4WQ(L)  = Ammonia (for Current Layer)
  ! *** RNO3WQ(L)  = Nitrate (for Current Layer)
  ! *** PO4DWQ(L)  = Phosphate (for Current Layer)
  ! *** RNH4NO3(L) = Total Inorganic Nitrogen (for Current Layer)

  ! *** WQKHNG = Nitrogen half-saturation for Algae-Greens (mg/L)
  ! *** WQKHPG = Phosphorus half-saturation for Algae-Greens (mg/L)
  ! *** WQKHCO2G = CO2 half-saturation for Algae-Greens (mg/L)

  ! *** XLIMIG = Rate Limiting Factor - Light
  ! *** XLIMTG = Rate Limiting Factor - Temperature   (Lookup Table: WQTDGG)
  ! *** XLIMNG = Rate Limiting Factor - Nitrogen      (Local-WQGNG)
  ! *** XLIMPG = Rate Limiting Factor - Phosphorus    (Local-WQGPG)
  ! *** XLIMCO2G = Rate Limiting Factor - CO2    (Local-WQGPG)
  ! *** WQF1NG = Rate Limiting Factor, Minimum of N & P

  ! *** WQPMG  = Maximum Growth Rate for Algae-Greens (1/d)
  ! *** WQPG   = Current Growth Rate for Algae-Greens (1/d)
  ! *** WQBMG  = Current Basal Metabolism Rate for Algae-Greens (1/d)
  ! *** WQPRG  = Current Predation Metabolism Rate for Algae-Greens (1/d)

  ! *** WQBMRG   = Basal Metabolism Rate for Algae-Greens (1/d)
  ! *** WQPRRG   = Predation Rate for Algae-Greens (1/d)
  ! *** WQTDRG   = Lookup Table for Temperature Rate Effect - Algae-Greens

  ! *** WQPC   = Final Net Growth Rate - Cyanobacteria
  ! *** WQPD   = Final Net Growth Rate - Diatoms Algae  
  ! *** WQPG   = Final Net Growth Rate - Green Algae
  ! *** WQPM   = Final Net Growth Rate - Macroalgae

  ! *** WQPNC = Preference for ammonium uptake - Cyanobacteria
  ! *** WQPND = Preference for ammonium uptake - Diatoms Algae
  ! *** WQPNG = Preference for ammonium uptake - Green Algae

  ! *** WQOBTOT  = Total Algal Biomass (mg/l)
  ! *** WQKRC    = Minimum Dissolution Rate of Refractory POC (1/day)
  ! *** WQKLC    = Minimum Dissolution Rate of Labile POC (1/day)
  ! *** WQKLCALG = Constant Refractory POC Dissolution Rate
  ! *** WQTDHDR  = Lookup Table for Temperature Rate Effect for Hydrolysis
  ! *** WQKRPC   = Current Dissolution Rate for POC

  ! *** WQI0   = SOLAR RADIATION for Current Time
  ! *** WQI1   = SOLAR RADIATION ON PREVIOUS DAY  
  ! *** WQI2   = SOLAR RADIATION TWO DAYS AGO  
  ! *** WQI3   = SOLAR RADIATION THREE DAYS AGO  

  ! *** WQKHR  = DOC Heterotrophic Respiration Rate

  ! *** WQWSSET = Water quality settling speed L;(L:1) is for top water layer; (L,2) is for lower water layers
  ! *** WQTTM   = Temporary concentration variable

  IF(.NOT.ALLOCATED(DELKC))THEN
    ALLOCATE(DELKC(KCM))  
    ALLOCATE(DZCHP(LCM))  
    ALLOCATE(TSSS_ABOVE(LCM))  
    ALLOCATE(CHLS_ABOVE(LCM))  
    ALLOCATE(POMS_ABOVE(LCM))  
    ALLOCATE(DZC_ABOVE(LCM))  
    ALLOCATE(WQISC(LCM))  
    ALLOCATE(WQISD(LCM))  
    ALLOCATE(WQISG(LCM))  
    ALLOCATE(WQI0TOP(LCM))  

    DO K=1,KC  
      DELKC(K)=0.  
    ENDDO  
    DELKC(KC)=1.  
    DZCHP=0.0
  ENDIF
  
  CNS1=2.718  
  NS=1  

  ! COMPUTE WQCHL,WQTAMP,WQPO4D,WQSAD AT A NEW TIME STEP: WQCHLX=1/WQCHLX  

  ! *** Compute WQCHL (Chlorophyll) Using Algal Biomass & factors
  DO K=1,KC  
    DO L=2,LA  
     WQCHL(L,K) = WQV(L,K,1)*WQCHLC + WQV(L,K,2)*WQCHLD + WQV(L,K,3)*WQCHLG  
    ENDDO  
  ENDDO  

  ! ***************************************************************************
  ! INITIALIZE SOLAR RADIATION AND OPTIMAL LIGHT

  ! *** INITIAL SOLAR RADIATION AT TOP OF SURFACE LAYER
  IF(USESHADE)THEN
    DO L=2,LA
      WQI0TOP(L)=WQI0 * PSHADE(L)
    ENDDO
  ELSE
    DO L=2,LA
      WQI0TOP(L)=WQI0
    ENDDO
  ENDIF
  ! ***  COMPUTE THE CURRENT OPTIMAL LIGHT INTENSITY
  IF(IWQSUN  ==  2)THEN  
    WQAVGIO = WQCIA*WQI1 + WQCIB*WQI2 + WQCIC*WQI3  
  ELSE
    WQAVGIO = WQCIA*WQI0 + WQCIB*WQI1 + WQCIC*WQI2
  ENDIF 
  ! *** CORRECT TO AVERAGE SOLAR RADIATION DURING DAYLIGHT HOURS

  ! ***************************************************************************
  ! *** DZWQ=1/H (for a layer), VOLWQ=1/VOL (m^-3)
  TSSS_ABOVE=0.0 
  CHLS_ABOVE=0.0
  POMS_ABOVE=0.0 
  DZC_ABOVE =0.0
  DO K=KC,1,-1  
    DO L=2,LA  
      TWQ(L)=TEM(L,K)              !layer temperature for WQ calcs
      SWQ(L)=MAX(SAL(L,K), 0.0)    !layer salinity for WQ calcs
      DZCHP(L)=DZC(K)*HP(L)        !layer thickness of a cell in meters
      DZWQ(L) = 1.0 / DZCHP(L)     !inverse layer thickness
      VOLWQ(L) = DZWQ(L) / DXYP(L) !inverse volume of each cell in a layer
      IMWQZT(L)=IWQZMAP(L,K)       !binary map for WQ calcs
    ENDDO  
        
    ! *** ZERO WQWPSL IF FLOWS ARE NEGATIVE.  THESE ARE HANDLED IN CALFQC (PMC)
    IF(IWQPSL /= 2)THEN
      DO NQ=1,NQSIJ  
        IF((QSERCELL(K,NQ)+QSS(K,NQ)) <= 0.0)THEN
          ! *** ZERO THE FLUX
          L=LQS(NQ)  
          DO NW=1,NWQV
            WQWPSL(L,K,NW)=0.0
          ENDDO
        ENDIF
      ENDDO
    ENDIF
            
    WQBCSET=0.0  
    
    !EVALUATING THE RATE OF ALGAE LEAVING THE CELL THROUGH SETTLING OR FLOATING
    DO L=2,LA  
      IF(WQWSC(IMWQZT(L)) < 0) THEN    !VB PERMITS CYANOBACTERIA TO FLOAT AND/OR SETTLE
        IF(K == KC) THEN      
          WQBCSET(L,1) = 0.0      !ALGAE AT THE WATER SURFACE CANT FLOAT INTO THE CELL ABOVE
        ELSE
          WQBCSET(L,1) = -WQWSC(IMWQZT(L))*DZWQ(L)   ! *** CYANOBACTERIA  !NEEDS TO BE A POSITIVE QTY  
        ENDIF
      ELSE
        WQBCSET(L,1) = WQWSC(IMWQZT(L))*DZWQ(L)   ! *** CYANOBACTERIA
      ENDIF

      IF(WQWSD(IMWQZT(L)) < 0) THEN    !VB PERMITS DIATOMS TO FLOAT AND/OR SETTLE
        IF(K == KC) THEN
          WQBDSET(L,1) = 0.0
        ELSE
          WQBDSET(L,1) = -WQWSD(IMWQZT(L))*DZWQ(L)   ! *** Diatoms  
        ENDIF
      ELSE
        WQBDSET(L,1) = WQWSD(IMWQZT(L))*DZWQ(L)   ! *** Diatoms
      ENDIF
      IF(WQWSG(IMWQZT(L)) < 0) THEN    !VB PERMITS GREEN ALGAE TO FLOAT AND/OR SETTLE
        IF(K == KC) THEN
          WQBGSET(L,1) = 0.0
        ELSE
          WQBGSET(L,1) = -WQWSG(IMWQZT(L))*DZWQ(L)   ! *** ALGAE  
        ENDIF
      ELSE
        WQBGSET(L,1) = WQWSG(IMWQZT(L))*DZWQ(L)   ! *** ALGAE
      ENDIF
    ENDDO


    ! *** ZONE SPECIFIC SETTING VELOCITIES, (m/day)   
    DO L=2,LA  
      !WQBCSET(L,1) = WQWSC(IMWQZT(L))*DZWQ(L)   ! *** Cyanobacteria   
      !WQBDSET(L,1) = WQWSD(IMWQZT(L))*DZWQ(L)   ! *** Diatoms
      !WQBGSET(L,1) = WQWSG(IMWQZT(L))*DZWQ(L)   ! *** Green
      WQRPSET(L,1) = WQWSRP(IMWQZT(L))*DZWQ(L)  ! *** Refractory POM 
      WQLPSET(L,1) = WQWSLP(IMWQZT(L))*DZWQ(L)  ! *** Labile POM 
    ENDDO

    ! *** SET SETTLING FOR TAM SORPTION: CURRENT LAYER  
    IF(IWQSRP == 1)THEN  
      DO L=2,LA  
        WQWSSET(L,1) = WQWSS(IMWQZT(L))*DZWQ(L)  
      ENDDO  
    ENDIF  

    DO L=2,LA  
      IF(K /= KC)THEN 
        IMWQZT1(L)=IWQZMAP(L,K+1)  
      ENDIF
      IF(K /= 1) THEN
        IMWQZT2(L)=IWQZMAP(L,K-1)  
      ENDIF
    ENDDO 

    DO L=2,LA
      IF(K /= KC)THEN 
        IF(WQWSC(IMWQZT1(L)) < 0) THEN
          WQBCSET(L,2) = 0.0  
        ELSE
          WQBCSET(L,2) = WQWSC(IMWQZT1(L))*DZWQ(L)  
        ENDIF
        IF(WQWSD(IMWQZT1(L)) < 0) THEN
          WQBDSET(L,2) = 0.0  
        ELSE
          WQBDSET(L,2) = WQWSD(IMWQZT1(L))*DZWQ(L)  
        ENDIF
        IF(WQWSG(IMWQZT1(L)) < 0) THEN
          WQBGSET(L,2) = 0.0  
        ELSE
          WQBGSET(L,2) = WQWSG(IMWQZT1(L))*DZWQ(L)  
        ENDIF
      ENDIF

      IF(K /= 1)THEN
        IF(WQWSC(IMWQZT2(L)) < 0) THEN
          WQBCSET(L,2) = WQBCSET(L,2)-WQWSC(IMWQZT2(L))*DZWQ(L)  
        !  ELSE
        !    WQBCSET(L,2) = WQBCSET(L,2) 0.0  
        ENDIF
        IF(WQWSD(IMWQZT2(L)) < 0) THEN
          WQBDSET(L,2) = WQBDSET(L,2)-WQWSD(IMWQZT2(L))*DZWQ(L) 
        !  ELSE
        !    WQBDSET(L,2) = 0.0  
        ENDIF
        IF(WQWSG(IMWQZT2(L)) < 0) THEN
          WQBGSET(L,2) = WQBGSET(L,2)-WQWSG(IMWQZT2(L))*DZWQ(L)  
        !  ELSE
        !    WQBGSET(L,2) = 0.0  
        ENDIF
      ENDIF
    ENDDO

    IF(K /= KC) THEN
      DO L=2,LA 
      !    WQBCSET(L,2) = WQWSC(IMWQZT1(L))*DZWQ(L)  
      !    WQBDSET(L,2) = WQWSD(IMWQZT1(L))*DZWQ(L)  
      !    WQBGSET(L,2) = WQWSG(IMWQZT1(L))*DZWQ(L)  
        WQRPSET(L,2) = WQWSRP(IMWQZT1(L))*DZWQ(L)  
        WQLPSET(L,2) = WQWSLP(IMWQZT1(L))*DZWQ(L)  
      ENDDO  
    ENDIF

    ! *** SET SETTLING FOR TAM SORPTION: ONE LAYER UP
    IF(IWQSRP == 1.AND.K /= KC)THEN  
      DO L=2,LA  
        WQWSSET(L,2) = WQWSS(IMWQZT1(L))*DZWQ(L)  
      ENDDO  
    ENDIF

    ! *** FIND AN INDEX FOR LOOK-UP TABLE FOR TEMPERATURE DEPENDENCY  
    DO L=2,LA  
      IWQT(L)=NINT((TWQ(L)-WQTDMIN)/WQTDINC)+1  
      IF(IWQT(L) < 1 .OR. IWQT(L) > NWQTD)THEN  
        OPEN(1,FILE='ERROR.LOG',POSITION='APPEND',STATUS='UNKNOWN')  
        WRITE(1,*)' *** ERROR IN WQSKE1:TEMPERATURE LOOKUP TABLE'
        WRITE(1,911) TIMEDAY, L, IL(L), JL(L), K, TWQ(L),TEM(L,K)  
        WRITE(6,600)IL(L),JL(L),K,TWQ(L)  
        IWQT(L)=MAX(IWQT(L),1)  
        IWQT(L)=MIN(IWQT(L),NWQTD)  
        CLOSE(1,STATUS='KEEP')
      ENDIF  
    ENDDO  

    600 FORMAT(' I,J,K,TEM = ',3I5,E13.4)  
    911 FORMAT('ERROR: TIME, L, I, J, K, TWQ(L),TEM(L,K) = ',F10.5, 4I4, 2F10.4,/)  

    !C NOTE: MRM 04/29/99  ADDED ARRAYS TO KEEP TRACK OF  
    !C       NITROGEN, PHOSPHORUS, LIGHT, AND TEMPERATURE LIMITS  
    !C       FOR ALGAE GROWTH FOR CYANOBACTERIA, DIATOMS, GREENS,  
    !C       AND MACROALGAE.  THESE ARE THE ARRAYS:  
    !C        XLIMNX(L,K) = NITROGEN    LIMITATION FOR ALGAE GROUP X  
    !C        XLIMPX(L,K) = PHOSPHORUS  LIMITATION FOR ALGAE GROUP X  
    !C        XLIMIX(L,K) = LIGHT       LIMITATION FOR ALGAE GROUP X  
    !C        XLIMTX(L,K) = TEMPERATURE LIMITATION FOR ALGAE GROUP X  
   
    ! *** BEGIN HORIZONTAL LOOP FOR ALGAE PARMETERS  
    DO L=2,LA
      RNH4WQ(L) = MAX (WQV(L,K,14), 0.0)  ! *** Ammonia
      RNO3WQ(L) = MAX (WQV(L,K,15), 0.0)  ! *** Nitrate
      PO4DWQ(L) = MAX (WQPO4D(L,K), 0.0)  ! *** Phosphate
      CO2WQ(L)  = MAX (WQV(L,K,22), 0.0)  ! *** CO2      !AA Added
      RNH4NO3(L) = RNH4WQ(L) + RNO3WQ(L)  ! *** Total Inorganic Nitrogen
      WQGNC = RNH4NO3(L) / (WQKHNC+RNH4NO3(L)+ 1.E-18)  
      WQGND = RNH4NO3(L) / (WQKHND+RNH4NO3(L)+ 1.E-18)  
      WQGNG = RNH4NO3(L) / (WQKHNG+RNH4NO3(L)+ 1.E-18)  
      WQGPC = PO4DWQ(L) / (WQKHPC+PO4DWQ(L)+ 1.E-18)  
      WQGPD = PO4DWQ(L) / (WQKHPD+PO4DWQ(L)+ 1.E-18)  
      WQGPG = PO4DWQ(L) / (WQKHPG+PO4DWQ(L)+ 1.E-18)  
      WQGCO2C = CO2WQ(L) / (ABS(WQKHCO2C)+CO2WQ(L)+ 1.E-18)    !AA
      WQGCO2D = CO2WQ(L) / (ABS(WQKHCO2D)+CO2WQ(L)+ 1.E-18)    !AA
      WQGCO2G = CO2WQ(L) / (ABS(WQKHCO2G)+CO2WQ(L)+ 1.E-18)    !AA

      XLIMNC(L,K) = XLIMNC(L,K) + WQGNC  
      XLIMND(L,K) = XLIMND(L,K) + WQGND  
      XLIMNG(L,K) = XLIMNG(L,K) + WQGNG  
      XLIMPC(L,K) = XLIMPC(L,K) + WQGPC  
      XLIMPD(L,K) = XLIMPD(L,K) + WQGPD  
      XLIMPG(L,K) = XLIMPG(L,K) + WQGPG
      XLIMCO2C(L,K) = XLIMCO2C(L,K) + WQGCO2C		!AA 
      XLIMCO2D(L,K) = XLIMCO2D(L,K) + WQGCO2D		!AA
      XLIMCO2G(L,K) = XLIMCO2G(L,K) + WQGCO2G		!AA

      IF(IDNOTRVA > 0 .AND. K == 1)THEN				! For macrophytes     
        WQGNM = RNH4NO3(L) / (WQKHNM+RNH4NO3(L) + 1.E-18)  
        WQGPM = PO4DWQ(L) / (WQKHPM+PO4DWQ(L) + 1.E-18) 
        WQGCO2M = CO2WQ(L) / (ABS(WQKHCO2M)+CO2WQ(L) + 1.E-18)	!AA
        IF(ISTRWQ(22) > 0)THEN
          WQF1NM = MIN(WQGNM, WQGPM, WQGCO2M)		!AA  Minimum of the N/P/CO2 Limit: Mats
        ELSE
          WQF1NM = MIN(WQGNM, WQGPM)				!AA  Minimum of the N/P     Limit: Mats
        ENDIF
        XLIMNM(L,K) = XLIMNM(L,K) + WQGNM  
        XLIMPM(L,K) = XLIMPM(L,K) + WQGPM
        XLIMCO2M(L,K) = XLIMCO2M(L,K) + WQGCO2M
      ENDIF

      IF(ISTRWQ(22) > 0)THEN 
        WQF1NC = MIN(WQGNC, WQGPC, WQGCO2C)		!AA  Minimum of the N/P/CO2 Limit: Cyanobacteria
      ELSE
        WQF1NC = MIN(WQGNC, WQGPC)				!AA  Minimum of the N/P     Limit: Cyanobacteria
      ENDIF

      IF(IWQSI == 1 .AND. ISTRWQ(22) > 0)THEN
        SADWQ = MAX (WQSAD(L,K), 0.0)
        WQGSD = SADWQ / (WQKHS+SADWQ+ 1.E-18)  
        WQF1ND = MIN(WQGND, WQGPD, WQGSD, WQGCO2D)	!AA  Minimum of the N/P/S/CO2 Limit: Diatoms
      ENDIF
      IF(IWQSI == 0 .AND. ISTRWQ(22) > 0)THEN  
        WQF1ND = MIN(WQGND, WQGPD, WQGCO2D)			!AA  Minimum of the N/P/CO2   Limit: Diatoms
      ENDIF
      IF(IWQSI == 0 .AND. ISTRWQ(22) == 0)THEN  
        WQF1ND = MIN(WQGND, WQGPD)				    !AA  Minimum of the N/P       Limit: Diatoms
      ENDIF 

      IF(ISTRWQ(22) > 0)THEN  
        WQF1NG = MIN(WQGNG, WQGPG, WQGCO2G)			!AA  Minimum of the N/P/CO2   Limit: Greens
      ELSE
        WQF1NG = MIN(WQGNG, WQGPG)				    !AA  Minimum of the N/P       Limit: Greens
      ENDIF

      IF(IDNOTRVA > 0)THEN  
        PO4DWQ(L) = MAX (WQPO4D(L,K), 0.0)  
      ENDIF  
      ! *** IN C&C, F2IC=F2IC/FCYAN, FACTOR TO ALLOW CYANOBACTERIA MAT FORMATION  
      ! *** LIGHT EXTINCTION (THIS WILL ALWAYS BE TRUE EXCEPT FOR IWQSUN=2)
      IF(WQI0 > 0.1)THEN
        DZC_ABOVE(L)=DZC_ABOVE(L)+DZC(K)
        IF(ISTRAN(6) > 0.OR.ISTRAN(7) > 0)THEN
          TSSS_ABOVE(L)=TSSS_ABOVE(L)+(SEDT(L,K)+SNDT(L,K))*DZC(K)
        ENDIF
        POMS_ABOVE(L)=POMS_ABOVE(L)+(WQV(L,K,4)+WQV(L,K,5))*DZC(K)
        ! *** COMPUTE TOTAL EXTINCTION COEFFICIENT
        WQKESS=WQKEB(IMWQZT(L))
        WQKESS=WQKESS+WQKETSS*(SEDT(L,K)+SNDT(L,K))  *DZCHP(L) !SCJ this seems to be the correct way to do this
        WQKESS=WQKESS+WQKEPOM*(WQV(L,K,4)+WQV(L,K,5))*DZCHP(L) !SCJ this seems to be the correct way to do this
        IF(WQKECHL < 0.0)THEN
          ! *** Compute Extinction Factor as a fn(Chla)
          XMRM = 0.054*WQCHL(L,K)**0.6667 + 0.0088*WQCHL(L,K)
        ELSE
          XMRM = WQKECHL*WQCHL(L,K)
          !XMRM = WQKECHL*WQV(L,K,3) !Fix for the LDRD greenhouse only
        ENDIF  
        WQKESS = WQKESS+XMRM  
        ! *** OPTIMAL LIGHT INTENSITY AT OPTIMAL DEPTH
        IF(K == KC)THEN
          WQISC(L) = MAX( WQAVGIO*EXP(-WQKESS*WQDOPC),WQISMIN)  
          WQISD(L) = MAX( WQAVGIO*EXP(-WQKESS*WQDOPD),WQISMIN)  
          WQISG(L) = MAX( WQAVGIO*EXP(-WQKESS*WQDOPG),WQISMIN) 
        ENDIF
        !WQTT1 = (CNS1 * WQFD * DZWQ(L)) / WQKETOT(L,K)	!Not needed - erroneous formulation WQTT1=(2.718*FD)/(Kess*dz)
        ! *** CURRENT LIGHT GROWTH LIMITING FACTOR

        IF(K == KC)THEN		!AA Initialize free surface intensity
          WQITOP(L,K) = WQI0TOP(L)
        ELSE				!AA Calculate surface intensity as a function of the layer above's surface intensity
          WQITOP(L,K) = WQITOP(L,K+1)*EXP(-WQKESS*DZCHP(L))
        ENDIF !SEE DiTORO ET AL (1971, EQNS. (11)&(12)) 
        EXPA0=EXP(-WQITOP(L,K)/WQISC(L))
        EXPA1=EXP(-WQITOP(L,K)/WQISC(L)*EXP(-WQKESS*DZCHP(L)))
        XLIMIC(L,K)=EXP(1.0)*WQFD/(DZCHP(L)*WQKESS)*(EXPA1-EXPA0)
        EXPA0=EXP(-WQITOP(L,K)/WQISD(L))
        EXPA1=EXP(-WQITOP(L,K)/WQISD(L)*EXP(-WQKESS*DZCHP(L)))
        XLIMID(L,K)=EXP(1.0)*WQFD/(DZCHP(L)*WQKESS)*(EXPA1-EXPA0)
        EXPA0=EXP(-WQITOP(L,K)/WQISG(L))
        EXPA1=EXP(-WQITOP(L,K)/WQISG(L)*EXP(-WQKESS*DZCHP(L)))
        XLIMIG(L,K)=EXP(1.0)*WQFD/(DZCHP(L)*WQKESS)*(EXPA1-EXPA0)

        WQF2IC = XLIMIC(L,K)
        WQF2ID = XLIMID(L,K)
        WQF2IG = XLIMIG(L,K)
      ELSE
        WQF2IC=0.0
        WQF2ID=0.0
        WQF2IG=0.0
      ENDIF
      ! *** UPDATE SOLAR RADIATION AT BOTTOM OF THIS LAYER  !AA
      ! *** MACROALGAE SUBMODEL - COMPUTE IF HAVE SOLAR RADIATION
      IF(IDNOTRVA > 0 .AND. K == 1.AND.WQI0TOP(L) > 1.0E-18)THEN
        IZ=IWQZMAP(L,K)  
        WQFDI0 = - WQI0TOP(L) / (WQFD + 1.E-18)
        WQISM = MAX( WQAVGIO*EXP(-WQKESS*WQDOPM(IZ)), WQISMIN )  
        WQFDM = WQFDI0 / (WQISM + 1.E-18)  
        WQF2IM = WQTT1 * (EXP(WQFDM*WQTTB) - EXP(WQFDM*WQTTT))  
        !WQF2IM = WQF2IM * PSHADE(L)  
        UMRM = 0.5*( U(L,K) + U(L+1   ,K) )  
        VMRM = 0.5*( V(L,K) + V(LNC(L),K) )  
        WQVEL=SQRT(UMRM*UMRM + VMRM*VMRM)  
        WQLVF=1.0  
 
        !C OPTION 1 FOR VELOCITY LIMITATION ASSUMES MACROALGAE GROWTH  
        !C IS LIMITED AT LOW VELOCITIES DUE TO REDUCED AVAILABILITY OF  
        !C NUTRIENTS REACHING THE ALGAE BIOMASS.  USES A MICHAELIS-MENTON  
        !C TYPE OF EQUATION.  
        IF(IWQVLIM == 1)THEN
          IF(WQVEL > WQKMVMIN(L))THEN
            WQLVF = WQVEL / (WQKMV(L) + WQVEL)  
          ELSE  
            WQLVF = WQKMVMIN(L) / (WQKMV(L) + WQKMVMIN(L))  
          ENDIF  
        ENDIF  
 
        !C OPTION 2 FOR VELOCITY LIMITATION APPLIES A FIVE-PARAMETER LOGISTIC  
        !C FUNCTION THAT CAN BE ADJUSTED TO LIMIT MACROALGAE GROWTH FOR  
        !C EITHER LOW OR HIGH (SCOUR) VELOCITIES.  IN STREAMS WITH LOW NUTRIENTS,  
        !C THE LOW VELOCITY WILL LIKELY BE LIMITING SINCE AMPLE NUTRIENTS MAY  
        !C NOT REACH THE ALGAE BIOMASS DUE TO REDUCED FLOW.  IN STREAMS WITH  
        !C ABUNDANT NUTRIENTS, LOW VELOCITIES WILL NOT LIMIT MACROALGAE GROWTH,  
        !C INSTEAD, HIGH VELOCITIES WILL LIKELY SCOUR THE MACROALGAE AND DETACH  
        !C IT FROM THE SUBSTRATE.  
        IF(IWQVLIM  == 2)THEN
          XNUMER = WQKMVA(L) - WQKMVD(L)  
          XDENOM = 1.0 + (WQVEL/WQKMVC(L))**WQKMVB(L)  
          WQLVF = WQKMVD(L) + ( XNUMER / (XDENOM**WQKMVE(L)) )  
        ENDIF  
 
        ! *** USE THE MORE SEVERELY LIMITING OF VELOCITY OR NUTRIENT FACTORS:  
        WQF1NM = MIN(WQLVF, WQF1NM)  
 
        !C FIRST CONVERT FROM MACROALGAE FROM A CONCENTRATION (MG C/M3)  
        !C TO A DENSITY (MG C/M2).  
        XMRM = WQV(L,K,IDNOTRVA)*DZCHP(L)  
        WQLDF = WQKBP(L) / (WQKBP(L) + XMRM)  
        WQPM(L)= WQPMM(IMWQZT(L))*WQF1NM*WQF2IM*WQTDGM(IWQT(L))*WQLDF
        XLIMVM(L,K) = XLIMVM(L,K) + WQLVF  
        XLIMDM(L,K) = XLIMDM(L,K) + WQLDF  
        XLIMIM(L,K) = XLIMIM(L,K) + WQF2IM  
        XLIMTM(L,K) = XLIMTM(L,K) + WQTDGM(IWQT(L))  
      ENDIF  
      XLIMTC(L,K) = XLIMTC(L,K) + WQTDGC(IWQT(L))  
      XLIMTD(L,K) = XLIMTD(L,K) + WQTDGD(IWQT(L))  
      XLIMTG(L,K) = XLIMTG(L,K) + WQTDGG(IWQT(L))

!  pH Limitation changes - VJ - Calculate the pH limiting factors if ISTRWQ(22)>1
      IF(ISTRWQ(22)  >  1)THEN !If CO2 is tracked, and ISTRWQ(22)>1 (e.g., 2) then there is also pH limitation
!  Calculate pH from aser.inp if WQKHCO2C is negative. Find current time, look into TASER array and find the two times between which current time falls.
        IF(WQKHCO2C<0.0)THEN
!  Then do a linear interpolation using the pH values at the above two times and calculate the pH at the current time.
          CTime = NITERAT*dt/86400.  ! Current time is calculated from global variable NITERAT couting the number of iterations.
          CPos = 0                   ! Index
          CPH = 0.0                  ! Interpolated pH
          DO CITER=1,MASER(NASER)    ! Iterating through the 1st time column in aser.inp and finding the position where current simulation time lies.
            IF (TASER(CITER,NASER) > CTime)THEN
              CPos = CITER
              EXIT 
            ENDIF 
          ENDDO
     ! Linear interpolation (PHVAL read from aser.inp)
          CPH = PHVAL(CPos-1,1) + (PHVAL(CPos,1) - PHVAL(CPos-1,1)) / (TASER(CPos,1) - TASER(CPos-1,1)) * (CTime - TASER(CPos-1,1))
          WQGPH(L,K)=CPH
        ELSE !Otherwise it is calculated from CO2 concentrations
  !  Calculate local pH as function of CO2
	      WQGPH(L,K) = -0.5*LOG10(1.008E-14 + (1.7E-3*4.45E-7*WQV(L,K,22)))
        ENDIF
        WQPHHCONC = 10.0**(-WQGPH(L,K)) ! Hydrogen ion concentration as function of pH
        WQPHKH = 1E-7*( -5.*(TWQ(L)**3) + 300.*(TWQ(L)**2) - 5000.*TWQ(L) + 3E4 ) !K_H as function of temp
        WQPHKOH = 1E-13*( 8.*(TWQ(L)**2) - 500.*TWQ(L) + 8000) ! K_OH as function of temp
        WQPHDC = WQPHHCONC / ( WQPHHCONC + WQPHKOH + ((WQPHHCONC**2)/WQPHKH) )
        WQPHDD = WQPHHCONC / ( WQPHHCONC + WQPHKOH + ((WQPHHCONC**2)/WQPHKH) )
        WQPHDG = WQPHHCONC / ( WQPHHCONC + WQPHKOH + ((WQPHHCONC**2)/WQPHKH) )
	    XLIMPHC(L,K) = XLIMPHC(L,K) + WQGPH(L,K)
	    XLIMPHD(L,K) = XLIMPHD(L,K) + WQGPH(L,K)
	    XLIMPHG(L,K) = XLIMPHG(L,K) + WQGPH(L,K)
	  ELSE !Even if CO2 is turned on, if ISTRWQ(22)=1, then there is no pH limitation, only CO2 limitation
	    WQPHDC=1.0
	    WQPHDD=1.0
	    WQPHDG=1.0
	  ENDIF

      ! *** Compute the Growth Rate based on Maximums & Limiting Factors
      IF(IWQSTOX == 1)THEN  
        WQF4SC = WQSTOX / (WQSTOX + SWQ(L)*SWQ(L)+1.E-12)  
        WQPC(L)=WQPMC(IMWQZT(L))*WQF1NC*WQF2IC*WQTDGC(IWQT(L))*WQF4SC*WQPHDC  
      ELSE  
        WQPC(L) = WQPMC(IMWQZT(L))*WQF1NC*WQF2IC*WQTDGC(IWQT(L))*WQPHDC  
      ENDIF  
      WQPD(L) = WQPMD(IMWQZT(L))*WQF1ND*WQF2ID*WQTDGD(IWQT(L))*WQPHDD
      WQPG(L) = WQPMG(IMWQZT(L))*WQF1NG*WQF2IG*WQTDGG(IWQT(L))*WQPHDG

! ALGAL BASAL METABOLISM & PREDATION  

      WQBMC(L) = WQBMRC(IMWQZT(L)) * WQTDRC(IWQT(L))  
      WQPRC(L) = WQPRRC(IMWQZT(L)) * WQTDRC(IWQT(L))  

! THE VARIABLE WQTDGP ADJUSTS PREDATION AND BASAL METABOLISM BASED ON A  
! LOWER/UPPER OPTIMUM TEMPERATURE FUNCTION.  THIS WILL ALLOW DIATOMS TO  
! BLOOM IN WINTER IF WQTDGP IS CLOSE TO ZERO.  

      WQBMD(L)=WQBMRD(IMWQZT(L))*WQTDRD(IWQT(L))*WQTDGP(IWQT(L))  
      WQPRD(L)=WQPRRD(IMWQZT(L))*WQTDRD(IWQT(L))*WQTDGP(IWQT(L))  
      WQBMG(L) = WQBMRG(IMWQZT(L)) * WQTDRG(IWQT(L))  
      WQPRG(L) = WQPRRG(IMWQZT(L)) * WQTDRG(IWQT(L))  
      IF(IDNOTRVA > 0.AND.K == 1)THEN
        WQBMM(L) = WQBMRM(IMWQZT(L)) * WQTDRM(IWQT(L))  
        WQPRM(L) = WQPRRM(IMWQZT(L)) * WQTDRM(IWQT(L))  
      ENDIF  
    ENDDO  ! *** END HORIZONTAL LOOP FOR ALGAE PARMETERS  

    ! *** 
    XMRM = 0.0  
    DO L=2,LA  
      IZ=IWQZMAP(L,K)  
      WQOBTOT(L) = WQV(L,K,1)+WQV(L,K,2)+WQV(L,K,3)  
      WQKRPC(L) = (WQKRC + WQKRCALG*WQOBTOT(L)) * WQTDHDR(IWQT(L))  
      WQKLPC(L) = (WQKLC + WQKLCALG*WQOBTOT(L)) * WQTDHDR(IWQT(L))  
      IF(IDNOTRVA > 0 .AND. K == 1)THEN
        XMRM = WQKDCALM(IZ) * WQV(L,K,IDNOTRVA)  
      ENDIF  

! M. MORTON 08/28/99: ADDED SPATIALLY VARIABLE DOC HYDROLYSIS RATE WQKDC  
!    TO ACHIEVE BETTER CONTROL IN SYSTEMS WITH A COMBINATION OF FRESHWAT  
!    STREAMS AND TIDAL RIVERS WITH DIFFERENT CHARACTERISTICS.  

      WQKDOC=(WQKDC(IZ)+WQKDCALG*WQOBTOT(L) + XMRM)*WQTDMNL(IWQT(L))  
      O2WQ(L) = MAX(WQV(L,K,19), 0.0)  
      WQTT1 = WQKDOC / (WQKHORDO + O2WQ(L)+ 1.E-18)  
      WQKHR(L) = WQTT1 * O2WQ(L)  
      WQDENIT(L)=WQTT1*WQAANOX*RNO3WQ(L)/(WQKHDNN+RNO3WQ(L)+ 1.E-18)  
    ENDDO  

! ***********************************
! 7-10 PHOSPHORUS  

    ! *** HYDROLYSIS  
    DO L=2,LA  
      WQAPC(L)=1.0/(WQCP1PRM+WQCP2PRM*EXP(-WQCP3PRM*PO4DWQ(L)))  
      WQKHP = (WQKHPC+WQKHPD+WQKHPG) / 3.0  
      WQTT1 = WQKHP / (WQKHP+PO4DWQ(L)+ 1.E-18) * WQOBTOT(L)  
      WQKRPP(L) = (WQKRP + WQKRPALG*WQTT1) * WQTDHDR(IWQT(L))  ! *** RPOP--> PO4
      WQKLPP(L) = (WQKLP + WQKLPALG*WQTT1) * WQTDHDR(IWQT(L))  ! *** LPOP--> DOP
      WQKDOP(L) = (WQKDP + WQKDPALG*WQTT1) * WQTDMNL(IWQT(L))  ! *** DOP --> PO4
    ENDDO
    ! *** PHOSPHATE SETTLING   
    DO L=2,LA  
      IF(IWQSRP == 1)THEN
        WQTTM = WQKPO4P*WQTAMP(L,K)  
        WQH10(L) = - WQWSSET(L,1) * WQTTM / (1.0+WQTTM)  
        IF(K /= KC)THEN
          WQTTM = WQKPO4P*WQTAMP(L,K+1)  
          WQT10(L) = WQWSSET(L,2) * WQTTM / (1.0+WQTTM)  
        ENDIF  
      ELSE IF(IWQSRP == 2)THEN
        WQTTS = WQKPO4P*SEDT(L,K)  
        WQH10(L) = - WSEDO(NS) * WQTTS * DZWQ(L) / (1.0+WQTTS)  
        IF(K /= KC)THEN
          WQTTS = WQKPO4P*SEDT(L,K)  
          WQT10(L) = WSEDO(NS) * WQTTS * DZWQ(L) / (1.0+WQTTS)  
        ENDIF  
      ELSE  
        WQH10(L) = 0.0  
        WQT10(L) = 0.0  
      ENDIF 
      WQH10(L) = WQH10(L)*DTWQO2 
    ENDDO  

! ***********************************
! 11-15 NITROGEN  

    ! *** HYDROLYSIS  
    DO L=2,LA  
      WQKHN = (WQKHNC+WQKHND+WQKHNG) / 3.0  
      WQTT1 = WQKHN / (WQKHN+RNH4NO3(L)+ 1.E-18) * WQOBTOT(L)  
      WQKRPN(L) = (WQKRN + WQKRNALG*WQTT1) * WQTDHDR(IWQT(L))  ! *** RPON-->NH3
      WQKLPN(L) = (WQKLN + WQKLNALG*WQTT1) * WQTDHDR(IWQT(L))  ! *** LON -->DON
      WQKDON(L) = (WQKDN + WQKDNALG*WQTT1) * WQTDMNL(IWQT(L))  ! *** DON -->NH3
    ENDDO  
    DO L=2,LA  
      IF(RNH4NO3(L) == 0.0)THEN
        WQPNC(L)=0.0  
        WQPND(L)=0.0  
        WQPNG(L)=0.0  
        WQPNM(L)=0.0  
      ELSE  
        WQTTC = RNH4WQ(L)/(WQKHNC+RNO3WQ(L)+ 1.E-18)  
        WQTTD = RNH4WQ(L)/(WQKHND+RNO3WQ(L)+ 1.E-18)  
        WQTTG = RNH4WQ(L)/(WQKHNG+RNO3WQ(L)+ 1.E-18)  
        WQTTM = RNH4WQ(L)/(WQKHNM+RNO3WQ(L)+ 1.E-18)  
        WQPNC(L) = (RNO3WQ(L)/(WQKHNC+RNH4WQ(L)+ 1.E-18) + WQKHNC/(RNH4NO3(L)+ 1.E-18)) * WQTTC  
        WQPND(L) = (RNO3WQ(L)/(WQKHND+RNH4WQ(L)+ 1.E-18) + WQKHND/(RNH4NO3(L)+ 1.E-18)) * WQTTD
        WQPNG(L) = (RNO3WQ(L)/(WQKHNG+RNH4WQ(L)+ 1.E-18) + WQKHNG/(RNH4NO3(L)+ 1.E-18)) * WQTTG
        WQPNM(L) = (RNO3WQ(L)/(WQKHNM+RNH4WQ(L)+ 1.E-18) + WQKHNM/(RNH4NO3(L)+ 1.E-18)) * WQTTM
      ENDIF  
      WQNIT(L) = O2WQ(L) * WQTDNIT(IWQT(L)) / ( (WQKHNDO+O2WQ(L)) * (WQKHNN+RNH4WQ(L)) + 1.E-18)
    ENDDO  
    IF(IWQSI == 1)THEN
      DO L=2,LA  
        IF(IWQSRP == 1)THEN
          WQTTM = WQKSAP*WQTAMP(L,K)  
          WQN17(L) = - WQWSSET(L,1) * WQTTM / (1.0+WQTTM)  
          IF(K /= KC)THEN
            WQTTM = WQKSAP*WQTAMP(L,K+1)  
            WQT17(L) = WQWSSET(L,2) * WQTTM / (1.0+WQTTM)  
          ENDIF  
        ELSE IF(IWQSRP == 2)THEN
          WQTTS = WQKSAP*SEDT(L,K)  
          WQN17(L) = - WSEDO(NS) * WQTTS * DZWQ(L) / (1.0+WQTTS)  
          IF(K /= KC)THEN
            WQTTS = WQKSAP*SEDT(L,K+1)  
            WQT17(L) = WSEDO(NS) * WQTTS * DZWQ(L) / (1.0+WQTTS)  
          ENDIF  
        ELSE  
          WQN17(L) = 0.0  
          WQT17(L) = 0.0  
        ENDIF  
      ENDDO  
      WQN17(L) = WQN17(L)*DTWQO2 
    ENDIF  

    ! ***********************************
    ! *** DISSOLVED OXYGEN
    PPCDO=-3.45	!PARTIAL PRES OF CO2 IN 10^ppcdo ATM; TEMPORARILY DECLARED HERE. SHLD BE READ IN FROM INPUT FILE
    DO L=2,LA  
      IZ=IWQZMAP(L,K)  
      WQO18(L)= -DTWQO2*WQKCOD(IWQT(L),IZ)*O2WQ(L) / (WQKHCOD(IZ) + O2WQ(L) + 1.E-18)

! *** DO Saturation, Modified by SCJ, see Garcia and Gordon, Limnology and Oceanography 37(6), 1992, Eqn. 8 and Table 1
      TVAL1=log((298.15-TWQ(L))/(273.15+TWQ(L)))
      TVAL2=TVAL1*TVAL1
      TVAL3=TVAL1*TVAL2
      TVAL4=TVAL1*TVAL3
      TVAL5=TVAL1*TVAL4
      RLNSAT1=5.80818+3.20684*TVAL1+4.11890*TVAL2+4.93845*TVAL3+1.01567*TVAL4+1.41575*TVAL5
      RLNSAT2=SWQ(L)*(-7.01211E-3-7.25958E-3*TVAL1-7.93334E-3*TVAL2-5.54491E-3*TVAL3)-1.32412E-7*SWQ(L)*SWQ(L)
      WQDOS(L) = EXP(RLNSAT1+RLNSAT2)*32E-3 !32E-3 approximately converts micromol/L to mg/L or g/m^3
      XDOSAT(L,K) = XDOSAT(L,K) + WQDOS(L)*DTWQ*DZCHP(L)
     
!************* CO2 parameters
!VB COMPUTING THE pK FOR SAT CONC OF CO2; K - HENRY'S CONST
     CDOSATIDX(L) = -2385.73/(TWQ(L) + 273.15) -  0.0152642 * (TWQ(L) + 273.15) + 14.0184
!          K * MOL WT OF CO2 * PARTAL PRES OF CO2 IN ATM
     WQCDOS(L) = 10.**(-CDOSATIDX(L)+PPCDO) * (44.* 1000.) !VB EVALUATING CONC OF CO2 IN G/M^3 
!************* CO2 parameters
      ! *** Compute Reaeration
      IF(K == KC)THEN
       

        ! DO NOT ALLOW WIND SPEEDS ABOVE 11 M/SEC IN THE FOLLOWING EQUATION
        WINDREA = MIN(WINDST(L),11.0)
        WQWREA=0.728*SQRT(WINDREA)+(0.0372*WINDREA-0.317)*WINDREA  
 
        IF(IWQKA(IZ) == 0)THEN
          WQVREA = WQKRO(IZ)  
          WQWREA = 0.0  
        ELSEIF(IWQKA(IZ) == 1)THEN  
          ! *** Constand plus Wind
          WQVREA = WQKRO(IZ)  
        ELSEIF(IWQKA(IZ) == 2)THEN
          UMRM = 0.5*(U(L,K)+U(L+1,K))  
          VMRM = 0.5*(V(L,K)+V(LNC(L),K))  
          XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)  
          ! *** WQKRO = 3.933 TYPICALLY  
          WQVREA = WQKRO(IZ) * XMRM**0.5 / HP(L)**0.5  
        ELSEIF(IWQKA(IZ) == 3)THEN
          UMRM = MAX(U(L,K), U(L+1,K))  
          VMRM = MAX(V(L,K), V(LNC(L),K))  
          XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)  
          ! *** WQKRO = 5.32 TYPICALLY  
          WQVREA = WQKRO(IZ) * XMRM**0.67 / HP(L)**1.85  
        ELSEIF(IWQKA(IZ) == 4)THEN
          ! *** MODIFIED OWENS AND GIBBS REAERATION EQUATION:  
          ! *** NOTE: NORMALIZED TO A DEPTH OF 1.0 FT, I.E., THIS EQUATION GIVES THE  
          ! ***       SAME REAERATION AS OWENS & GIBBS AT 1.0 FT DEPTH; AT HIGHER  
          ! ***       DEPTHS IT GIVES LARGER REAERATION THAN OWENS & GIBBS.  
          UMRM = MAX(U(L,K), U(L+1,K))  
          VMRM = MAX(V(L,K), V(LNC(L),K))  
          XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)  
          YMRM = HP(L)*3.0*(1.0 - HP(L)/(HP(L)+0.1524))  
          ! *** WQKRO = 5.32 TYPICALLY  
          WQVREA = WQKRO(IZ) * XMRM**0.67 / YMRM**1.85  
        ELSEIF(IWQKA(IZ) == 5)THEN
          UMRM = MAX(U(L,K), U(L+1,K))  
          VMRM = MAX(V(L,K), V(LNC(L),K))  
          XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)  
          WQVREA = 3.7*XMRM  
        ENDIF  

        ! *** NOW COMBINE REAERATION DUE TO WATER VELOCITY AND WIND STRESS
        WQVREA = WQVREA * REAC(IZ)  ! *** Diffusive Flux
        WQWREA = WQWREA * REAC(IZ)  ! *** Wind Component
        WQP19(L) = - (WQVREA + WQWREA) * DZWQ(L)* WQTDKR(IWQT(L),IZ)  
        WQKRDOS(L) = -WQP19(L)*WQDOS(L)
        WQP22(L) = WQP19(L)*((32./44.)**0.25) 		!VB Kr FOR CO2 ANALOGOUS TO WQP19 ; 44 = MOL WT OF CO2
        WQKRCDOS(L) = -WQP22(L) * WQCDOS(L)		!VB EVALUATING Kr*SAT CONC OF CO2
      ELSE  
        WQP19(L) = 0.0  
        WQP22(L) = 0.0							!VB Kr FOR CO2 IS ZERO FOR CELLS NOT AT THE SURFACE
      ENDIF  
    ENDDO  
    ! ***********************************
    IF(IWQSRP == 1)THEN
      DO L=2,LA  
        WQR20(L) = WQWPSL(L,K,20)*VOLWQ(L) + (WQV(L,K,20) - WQTAMP(L,K)) * WQWSSET(L,1)
      ENDDO

      IF(K == 1)THEN
        DO L=2,LA  
          IF(LMASKDRY(L))THEN
            WQR20(L) = WQR20(L) + WQTDTAM(IWQT(L))*DZWQ(L)/(WQKHBMF+O2WQ(L)+ 1.E-18)
          ENDIF
        ENDDO
      ENDIF

      IF(K /= KC)THEN
        DO L=2,LA  
          WQR20(L) = WQR20(L) + (WQV(L,K+1,20) - WQTAMP(L,K+1)) * WQWSSET(L,2)
        ENDDO 
      ELSE     ! K == KC
        WQR20(L)=WQR20(L)+(WQWDSL(L,KC,20)+WQATML(L,KC,20))*VOLWQ(L)
      ENDIF 
    ENDIF  
    ! ***********************************
    ! *** MACROPHYTE
    IF(IDNOTRVA > 0.AND.K == 1)THEN  
      DO L=2,LA  
        WQA1C = (WQPM(L) - WQBMM(L) - WQPRM(L)-WQWSM*DZWQ(L))*DTWQO2  
        WQVA1C = 1.0 / (1.0 - WQA1C)  
        WQV(L,K,IDNOTRVA)=(WQV(L,K,IDNOTRVA)+WQA1C*WQV(L,K,IDNOTRVA))*WQVA1C*SMAC(L)
        WQV(L,K,IDNOTRVA) = MAX(WQV(L,K,IDNOTRVA),WQMCMIN)*SMAC(L)  
        WQO(L,IDNOTRVA) = WQVO(L,K,IDNOTRVA)+WQV(L,K,IDNOTRVA)  
      ENDDO  
    ENDIF  
!******************************************************************************

! *** NOW COMPUTE KINETICS FOR EACH CONSTITUENT

! ****  PARAM 01  CHC - cyanobacteria

    IF(ISTRWQ(1) == 1)THEN  
      DO L=2,LA  
        ! *** GROWTH BASAL_METAB PREDATION SETTLING  TIME STEP  
        WQA1C=(WQPC(L)-WQBMC(L)-WQPRC(L)-WQBCSET(L,1))*DTWQO2 !production per unit time multiplied by half time step
        WQKK(L) = 1.0 / (1.0 - WQA1C)  !a number?

        ! ***   PT_SRC_LOADS    VOLUME  
        WQR1C = WQWPSL(L,K,1) * VOLWQ(L)  !point source load rate multiplied by inverse cell volume  g/m^3/t
        WQRR(L) = WQV(L,K,1) + DTWQ*WQR1C + WQA1C*WQV(L,K,1)   !transported biomass conc. (CALWQC) + point source load rate X time step + growth rate X previous biomass conc.
      ENDDO  

      IF(K /= KC)THEN  
        ! *** Add in settling from above
        DO L=2,LA  
          WQRR(L) = WQRR(L) + DTWQO2*WQBCSET(L,2)*WQO(L,1) !biomass conc. + DtX(1/t)* biomass conc.
        ENDDO
      ELSE  ! K == KC
        DO L=2,LA  
          ! ***    ATM DRY DEP   ATM WET DEP    VOLUME  
          WQR1C = (WQWDSL(L,KC,1)+WQATML(L,KC,1))*VOLWQ(L)  !atmospheric loading mass per time / cell volume
          WQRR(L) = WQRR(L) + DTWQ*WQR1C  !biomass conc. + Dt*loading rate per unit volume
        ENDDO
      ENDIF  
      DO L=2,LA  
        WQV(L,K,1)=SCB(L)*(WQRR(L)*WQKK(L))+(1.0-SCB(L))*WQV(L,K,1)  !boundary condition implementation
        WQO(L,1)= WQVO(L,K,1)+WQV(L,K,1) !depth totaled biomass conc = old biomass conc in cell + biomass conc from this iteration
      ENDDO  
    ENDIF  
! ****  PARAM 02  CHD - diatom algae

    IF(ISTRWQ(2) == 1)THEN  
      DO L=2,LA  
        ! *** GROWTH BASAL_METAB PREDATION SETTLING  TIME STEP  
        WQA2D=(WQPD(L)-WQBMD(L)-WQPRD(L)-WQBDSET(L,1))*DTWQO2
        WQKK(L) = 1.0 / (1.0 - WQA2D)  

        ! ***   PT_SRC_LOADS    VOLUME
        WQR2D = WQWPSL(L,K,2) * VOLWQ(L)  
        WQRR(L) = WQV(L,K,2) + DTWQ*WQR2D + WQA2D*WQV(L,K,2)  
      ENDDO  

      IF(K /= KC)THEN  
        ! *** Add in settling from above
        DO L=2,LA  
          WQRR(L) = WQRR(L) + DTWQO2*WQBDSET(L,2)*WQO(L,2)
        ENDDO 
      ELSE  ! K == KC
        DO L=2,LA  
          ! ***    ATM DRY DEP   ATM WET DEP    VOLUME
          WQR2D = (WQWDSL(L,KC,2)+WQATML(L,KC,2))*VOLWQ(L)  
          WQRR(L) = WQRR(L) + DTWQ*WQR2D  
        ENDDO
      ENDIF  
      DO L=2,LA  
        WQV(L,K,2)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,2)  
        WQO(L,2)=WQVO(L,K,2)+WQV(L,K,2)
      ENDDO  
    ENDIF  

! ****  PARAM 03  CHG - green algae

    IF(ISTRWQ(3) == 1)THEN
      DO L=2,LA  
        ! *** GROWTH BASAL_METAB PREDATION SETTLING  TIME STEP  
        WQA3G=(WQPG(L)-WQBMG(L)-WQPRG(L)-WQBGSET(L,1))*DTWQO2
        WQKK(L) = 1.0 / (1.0 - WQA3G)  

        ! ***   PT_SRC_LOADS    VOLUME  
        WQR3G = WQWPSL(L,K,3) * VOLWQ(L) 
        ! ***                   External      Internal 
        WQRR(L) = WQV(L,K,3) + DTWQ*WQR3G + WQA3G*WQV(L,K,3)  
      ENDDO  
      IF(K /= KC)THEN
        ! *** Add the Algae settled in from the cell above  
        DO L=2,LA  
          WQRR(L) = WQRR(L) + DTWQO2*WQBGSET(L,2)*WQO(L,3)
        ENDDO  
      ELSE  ! K == KC
        DO L=2,LA  
          ! ***    ATM DRY DEP   ATM WET DEP     VOLUME
          WQR3G = (WQWDSL(L,KC,3)+WQATML(L,KC,3))*VOLWQ(L)  
          WQRR(L) = WQRR(L) + DTWQ*WQR3G  
        ENDDO
      ENDIF  
      DO L=2,LA  
        WQV(L,K,3)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,3)  
        WQO(L,3)=WQVO(L,K,3)+WQV(L,K,3)
      ENDDO  
    ENDIF  

! ****  PARAM 04  ROC - refractory particulate organic carbon

    IF(ISTRWQ(4) == 1)THEN  
      DO L=2,LA  
        ! ***     HYDROLYSIS    SETTLING  
        WQB4 = -( WQKRPC(L) + WQRPSET(L,1))*DTWQO2
        WQKK(L) = 1.0 / (1.0 - WQB4)  
        ! ***  ALGAE PREDATION SOURCE OF RPOC  
        WQA4 = WQFCRP * (WQPRC(L)*WQO(L,1) + WQPRD(L)*WQO(L,2) + WQPRG(L)*WQO(L,3))
        IF(IDNOTRVA > 0.AND.K == 1)THEN  
          WQA4 = WQA4+WQFCRPM*WQPRM(L)*WQVO(L,K,IDNOTRVA)  
        ENDIF  

        ! ***  PT_SRC_LOADS    VOLUME
        WQR4 = WQWPSL(L,K,4) * VOLWQ(L)  
        WQRR(L) = WQV(L,K,4) + DTWQ*WQR4 + DTWQO2*WQA4 + WQB4*WQV(L,K,4)
      ENDDO

      IF(K /= KC)THEN  
        ! *** Add in settling from above
        DO L=2,LA  
          WQRR(L) = WQRR(L) + DTWQO2*WQRPSET(L,2)*WQO(L,4)  
        ENDDO  
      ELSE  ! K == KC
        DO L=2,LA  
          ! ***    ATM DRY DEP   ATM WET DEP    VOLUME
          WQR4 = (WQWDSL(L,KC,4)+WQATML(L,KC,4))*VOLWQ(L)  
          WQRR(L) = WQRR(L) + DTWQ*WQR4  
        ENDDO
      ENDIF  
      DO L=2,LA  
        WQV(L,K,4)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,4)  
        WQO(L,4)=WQVO(L,K,4)+WQV(L,K,4)
      ENDDO  
    ENDIF  

! ****  PARAM 05  LOC - labile particulate organic carbon

    IF(ISTRWQ(5) == 1)THEN  
      DO L=2,LA  
        ! ***     HYDROLYSIS    SETTLING  
        WQC5 = - (WQKLPC(L)  + WQLPSET(L,1))*DTWQO2 
        WQKK(L) = 1.0 / (1.0 - WQC5)
        !Predation
        WQA5 = WQFCLP * (WQPRC(L)*WQO(L,1) + WQPRD(L)*WQO(L,2) + WQPRG(L)*WQO(L,3))  

        IF(IDNOTRVA > 0.AND.K == 1)THEN  
          WQA5 =WQA5 + WQFCLPM * WQPRM(L)*WQVO(L,K,IDNOTRVA)  
        ENDIF  
        ! ***  PT_SRC_LOADS    VOLUME
        WQR5 = WQWPSL(L,K,5) * VOLWQ(L)

        WQRR(L) = WQV(L,K,5) + DTWQ*WQR5 + DTWQO2*WQA5 + WQC5*WQV(L,K,5)
      ENDDO  
      IF(K /= KC)THEN  
        ! *** Add in settling from above
        DO L=2,LA  
          WQRR(L) = WQRR(L) + DTWQO2*WQLPSET(L,2)*WQO(L,5)  
        ENDDO  
      ELSE  ! K == KC
        DO L=2,LA  
          ! ***    ATM DRY DEP   ATM WET DEP    VOLUME
          WQR5 = (WQWDSL(L,K,5)+WQATML(L,KC,5))*VOLWQ(L)  
          WQRR(L) = WQRR(L) + DTWQ*WQR5  
        ENDDO
      ENDIF  
      DO L=2,LA  
        WQV(L,K,5)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,5)  
        WQO(L,5)=WQVO(L,K,5)+WQV(L,K,5)
      ENDDO  
    ENDIF  

! ****  PARAM 06  DOC - dissolved organic carbon
    IF(ISTRWQ(6) == 1)THEN  
      DO L=2,LA  
        ! ***    RESPIRATION  DENITRIFICATION
        WQD6 = - ( WQKHR(L) +   WQDENIT(L)) *DTWQO2
        WQKK(L) = 1.0 / (1.0 - WQD6)  
        WQA6C=WQFCDC + CFCDCWQ*( WQKHRC/(WQKHRC+O2WQ(L)+ 1.E-18) )  
        WQA6D=WQFCDD + CFCDDWQ*( WQKHRD/(WQKHRD+O2WQ(L)+ 1.E-18) )  
        WQA6G=WQFCDG + CFCDGWQ*( WQKHRG/(WQKHRG+O2WQ(L)+ 1.E-18) )  
        WQA6 = ( WQA6C*WQBMC(L) + WQFCDP*WQPRC(L) )*WQO(L,1) + ( WQA6D*WQBMD(L) + WQFCDP*WQPRD(L) )*WQO(L,2) + ( WQA6G*WQBMG(L) + WQFCDP*WQPRG(L) )*WQO(L,3)
        IF(IDNOTRVA > 0.AND.K == 1)THEN  
          IZ=IWQZMAP(L,K)  
          WQA6M=(WQFCDM+(1.-WQFCDM)*WQKHRM(IZ) / (WQKHRM(IZ) + O2WQ(L) + 1.E-18))*WQBMM(L)
          WQA6 =WQA6+ (WQA6M+ WQFCDPM*WQPRM(L))*WQVO(L,K,IDNOTRVA)  
        ENDIF  
        ! ***  PT_SRC_LOADS    VOLUME
        WQR6 = WQWPSL(L,K,6) * VOLWQ(L)
        WQRR(L) = WQV(L,K,6) + DTWQ*WQR6 + WQD6*WQV(L,K,6) + DTWQO2*(WQA6 +WQKLPC(L)*WQO(L,5))
      ENDDO
      IF(K == KC)THEN
        DO L=2,LA  
          ! ***    ATM DRY DEP   ATM WET DEP    VOLUME
          WQR6 = (WQWDSL(L,K,6)+WQATML(L,KC,6))*VOLWQ(L)  
          WQRR(L) = WQRR(L) + DTWQ*WQR6  
        ENDDO
      ENDIF
      DO L=2,LA  
        WQV(L,K,6)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,6)  
        WQO(L,6)=WQVO(L,K,6)+WQV(L,K,6)
      ENDDO  
    ENDIF  
! ****  PARAM 07  ROP - refractory particulate organic phosphorus
    IF(ISTRWQ(7) == 1)THEN  
      DO L=2,LA  
        WQE7 = - (WQKRPP(L)+WQRPSET(L,1))*DTWQO2 
        WQKK(L) = 1.0 / (1.0 - WQE7)  
        WQA7C = (WQFPRC*WQBMC(L) + WQFPRP*WQPRC(L)) * WQO(L,1)  
        WQA7D = (WQFPRD*WQBMD(L) + WQFPRP*WQPRD(L)) * WQO(L,2)  
        WQA7G = (WQFPRG*WQBMG(L) + WQFPRP*WQPRG(L)) * WQO(L,3)  
        WQA7 = (WQA7C+WQA7D+WQA7G) * WQAPC(L)  
        IF(IDNOTRVA > 0.AND.K == 1)THEN  
          WQA7 = WQA7 + (WQFPRM*WQBMM(L) + WQFPRPM*WQPRM(L)) * WQVO(L,K,IDNOTRVA)* WQAPC(L)*WQAPCM
        ENDIF  
        ! ***  PT_SRC_LOADS    VOLUME
        WQR7 = WQWPSL(L,K,7) * VOLWQ(L)  
        WQRR(L) = WQV(L,K,7) + DTWQ*WQR7 + DTWQO2*WQA7 + WQE7*WQV(L,K,7)
      ENDDO  
      IF(K /= KC)THEN  
        ! *** Add in settling from above
        DO L=2,LA  
          WQRR(L) = WQRR(L) + DTWQO2*WQRPSET(L,2)*WQO(L,7)
        ENDDO  
      ELSE
        DO L=2,LA  
          ! ***    ATM DRY DEP   ATM WET DEP    VOLUME
          WQR7 = (WQWDSL(L,K,7)+WQATML(L,KC,7))*VOLWQ(L)  
          WQRR(L) = WQRR(L) + DTWQ*WQR7  
        ENDDO
      ENDIF  
      DO L=2,LA  
        WQV(L,K,7)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,7)  
        WQO(L,7)=WQVO(L,K,7)+WQV(L,K,7)
      ENDDO  
    ENDIF  
! ****  PARAM 08  LOP - labile particulate organic phosphorus
    IF(ISTRWQ(8) == 1)THEN  
      DO L=2,LA  
        ! ***    HYDROLYSIS  SETTLING
        WQF8 = - (WQKLPP(L)+WQLPSET(L,1))*DTWQO2
        WQKK(L) = 1.0 / (1.0 - WQF8)  
        WQA8C = (WQFPLC*WQBMC(L) + WQFPLP*WQPRC(L)) * WQO(L,1)  
        WQA8D = (WQFPLD*WQBMD(L) + WQFPLP*WQPRD(L)) * WQO(L,2)  
        WQA8G = (WQFPLG*WQBMG(L) + WQFPLP*WQPRG(L)) * WQO(L,3)  
        WQA8 = (WQA8C+WQA8D+WQA8G) * WQAPC(L)  
        IF(IDNOTRVA > 0.AND.K == 1)THEN  
          WQA8 = WQA8 +     (WQFPLM*WQBMM(L) + WQFPLPM*WQPRM(L)) * WQVO(L,K,IDNOTRVA)* WQAPC(L)*WQAPCM
        ENDIF  
        ! ***  PT_SRC_LOADS    VOLUME
        WQR8 = WQWPSL(L,K,8) * VOLWQ(L)  
        WQRR(L) = WQV(L,K,8) + DTWQ*WQR8 + DTWQO2*WQA8 + WQF8*WQV(L,K,8)
      ENDDO  
      IF(K /= KC)THEN  
        ! *** Add in settling from above
        DO L=2,LA  
          WQRR(L) = WQRR(L) + DTWQO2*WQLPSET(L,2)*WQO(L,8)
        ENDDO  
      ELSE
        DO L=2,LA  
          ! ***    ATM DRY DEP   ATM WET DEP    VOLUME
          WQR8 = (WQWDSL(L,K,8)+WQATML(L,KC,8))*VOLWQ(L)  
          WQRR(L) = WQRR(L) + DTWQ*WQR8  
        ENDDO
      ENDIF  
      DO L=2,LA  
        WQV(L,K,8)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,8)  
        WQO(L,8)=WQVO(L,K,8)+WQV(L,K,8)
      ENDDO  
    ENDIF  
! ****  PARAM 09  DOP - dissolved organic phosphorus
    IF(ISTRWQ(9) == 1)THEN  
      DO L=2,LA
        WQF9 = - DTWQO2*WQKDOP(L)  
        WQKK(L) = 1.0 / (1.0 - WQF9)
        WQA9C = (WQFPDC*WQBMC(L) + WQFPDP*WQPRC(L)) * WQO(L,1)  
        WQA9D = (WQFPDD*WQBMD(L) + WQFPDP*WQPRD(L)) * WQO(L,2)  
        WQA9G = (WQFPDG*WQBMG(L) + WQFPDP*WQPRG(L)) * WQO(L,3)  
        WQA9 = (WQA9C+WQA9D+WQA9G) * WQAPC(L)  
        IF(IDNOTRVA > 0.AND.K == 1)THEN  
          WQA9 = WQA9 + (WQFPDM*WQBMM(L) + WQFPDPM*WQPRM(L)) * WQVO(L,K,IDNOTRVA) * WQAPC(L)*WQAPCM
        ENDIF  
        ! ***  PT_SRC_LOADS    VOLUME
        WQR9 = WQWPSL(L,K,9) * VOLWQ(L)  
        WQRR(L) = WQV(L,K,9) + DTWQ*WQR9 + WQF9*WQV(L,K,9) + DTWQO2*(WQA9 + WQKLPP(L)*WQO(L,8) )
      ENDDO
      IF(K == KC)THEN
        DO L=2,LA  
          ! ***    ATM DRY DEP    ATM WET DEP    VOLUME
          WQR9 = (WQWDSL(L,KC,9)+WQATML(L,KC,9))*VOLWQ(L)  
          WQRR(L) = WQRR(L) + DTWQ*WQR9  
        ENDDO
      ENDIF
      DO L=2,LA  
        WQV(L,K,9)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,9)  
      ENDDO  
    ENDIF  

! ****  PARAM 10  P4D - total phosphate
    IF(ISTRWQ(10) == 1)THEN  
      DO L=2,LA  
        WQA10C=(WQFPIC*WQBMC(L)+WQFPIP*WQPRC(L)-WQPC(L))*WQO(L,1)  
        WQA10D=(WQFPID*WQBMD(L)+WQFPIP*WQPRD(L)-WQPD(L))*WQO(L,2)  
        WQA10G=(WQFPIG*WQBMG(L)+WQFPIP*WQPRG(L)-WQPG(L))*WQO(L,3)  
        WQKK(L) = (WQA10C+WQA10D+WQA10G) * WQAPC(L)  
        IF(IDNOTRVA > 0.AND.K == 1)THEN  
          WQKK(L) =WQKK(L)+(WQFPIM*WQBMM(L)+WQFPIP*WQPRM(L)-WQPM(L)) *WQVO(L,K,IDNOTRVA) * WQAPC(L)*WQAPCM
        ENDIF  
        ! ***    PT_SRC_LOADS    VOLUME
        WQRR(L) = WQWPSL(L,K,10) * VOLWQ(L)  
      ENDDO  
      IF(K == 1)THEN  
        DO L=2,LA 
          IF(LMASKDRY(L))THEN 
            WQRR(L) = WQRR(L) + WQBFPO4D(L)*DZWQ(L) ! *** Add in Benthic Flux
          ENDIF
        ENDDO  
      ENDIF
      IF(K == KC)THEN
        DO L=2,LA  
          ! ***      ATM DRY DEP    ATM WET DEP     VOLUME
          WQR10 = (WQWDSL(L,KC,10)+WQATML(L,KC,10))*VOLWQ(L)  
          WQRR(L) = WQRR(L) + WQR10  
        ENDDO
      ENDIF  
      DO L=2,LA  
         WQRR(L) = WQV(L,K,10) + DTWQ*WQRR(L) + WQH10(L)*WQV(L,K,10) + DTWQO2*(WQKK(L) + WQKRPP(L)*WQO(L,7) + WQKDOP(L)*WQO(L,9))
      ENDDO  
      IF(K /= KC)THEN  
        ! *** Add in settling from above
        DO L=2,LA  
          WQRR(L) = WQRR(L) + DTWQO2*WQT10(L)*WQO(L,10)
        ENDDO  
      ENDIF
      WQV10LAST=WQV(3,1,10)  
      DO L=2,LA  
        WQKKL = 1.0 / (1.0 - WQH10(L)) 
        WQV(L,K,10)=SCB(L)*(WQRR(L)*WQKKL)+(1.-SCB(L))*WQV(L,K,10)  
        WQO(L,10)=WQVO(L,K,10)+WQV(L,K,10)
      ENDDO  
    ENDIF  

! ****  PARAM 11  RON - refractory particulate organic nitrogen
    IF(ISTRWQ(11) == 1)THEN  
      DO L=2,LA  
        ! ***     HYDROLYSIS     SETTLING
        WQI11 = - (WQKRPN(L) + WQRPSET(L,1))*DTWQO2
        WQKK(L) = 1.0 / (1.0 - WQI11)  
        WQA11C=(WQFNRC*WQBMC(L)+WQFNRP*WQPRC(L))*WQANCC*WQO(L,1)  
        WQA11D=(WQFNRD*WQBMD(L)+WQFNRP*WQPRD(L))*WQANCD*WQO(L,2)  
        WQA11G=(WQFNRG*WQBMG(L)+WQFNRP*WQPRG(L))*WQANCG*WQO(L,3)  
        WQA11 = WQA11C+WQA11D+WQA11G  
        IF(IDNOTRVA > 0.AND.K == 1)THEN  
          WQA11 =WQA11 + (WQFNRM*WQBMM(L)+WQFNRPM*WQPRM(L))*WQANCM*WQVO(L,K,IDNOTRVA)
        ENDIF  

        ! ***    PT_SRC_LOADS    VOLUME
        WQR11 = WQWPSL(L,K,11) * VOLWQ(L)
        WQRR(L) = WQV(L,K,11) + DTWQ*WQR11 + DTWQO2*WQA11 + WQI11*WQV(L,K,11)
      ENDDO 
      IF(K /= KC)THEN  
        ! *** Add in settling from above
        DO L=2,LA  
          WQRR(L) = WQRR(L) + DTWQO2*WQRPSET(L,2)*WQO(L,11)
        ENDDO  
      ELSE   ! K == KC
        DO L=2,LA  
          ! ***     ATM DRY DEP     ATM WET DEP    VOLUME
          WQR11 = (WQWDSL(L,KC,11)+WQATML(L,KC,11))*VOLWQ(L)  
          WQRR(L) = WQRR(L) + DTWQ*WQR11  
        ENDDO
      ENDIF
      DO L=2,LA  
        WQV(L,K,11)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,11)  
        WQO(L,11)=WQVO(L,K,11)+WQV(L,K,11)
      ENDDO  
    ENDIF  
! ****  PARAM 12  LON - labile particulate organic nitrogen
    IF(ISTRWQ(12) == 1)THEN  
      DO L=2,LA  
        ! ***     HYDROLYSIS     SETTLING
        WQJ12 = - (WQKLPN(L)+WQLPSET(L,1))*DTWQO2
        WQKK(L) = 1.0 / (1.0 - WQJ12)  
        WQA12C=(WQFNLC*WQBMC(L)+WQFNLP*WQPRC(L))*WQANCC*WQO(L,1)  
        WQA12D=(WQFNLD*WQBMD(L)+WQFNLP*WQPRD(L))*WQANCD*WQO(L,2)  
        WQA12G=(WQFNLG*WQBMG(L)+WQFNLP*WQPRG(L))*WQANCG*WQO(L,3)  
        WQA12 = WQA12C+WQA12D+WQA12G  
        IF(IDNOTRVA > 0.AND.K == 1)THEN  
          WQA12 =WQA12 +(WQFNLM*WQBMM(L)+WQFNLPM*WQPRM(L)) *WQANCM*WQVO(L,K,IDNOTRVA)
        ENDIF  

        ! ***    PT_SRC_LOADS    VOLUME
        WQR12 = WQWPSL(L,K,12) * VOLWQ(L)  
        WQRR(L) = WQV(L,K,12) + DTWQ*WQR12 + DTWQO2*WQA12 + WQJ12*WQV(L,K,12)
      ENDDO  
      IF(K /= KC)THEN  
        ! *** Add in settling from above
        DO L=2,LA  
          WQRR(L) = WQRR(L) + DTWQO2*WQLPSET(L,2)*WQO(L,12)
        ENDDO  
      ELSE   ! K == KC
        DO L=2,LA  
          ! ***     ATM DRY DEP     ATM WET DEP    VOLUME
          WQR12 = (WQWDSL(L,KC,12)+WQATML(L,KC,12))*VOLWQ(L)  
          WQRR(L) = WQRR(L) + DTWQ*WQR12  
        ENDDO
      ENDIF  
      DO L=2,LA  
        WQV(L,K,12)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,12)  
        WQO(L,12)=WQVO(L,K,12)+WQV(L,K,12)
      ENDDO  
    ENDIF  
! ****  PARAM 13  DON - dissolved organic nitrogen
    IF(ISTRWQ(13) == 1)THEN  
      DO L=2,LA
        WQF13 = - DTWQO2*WQKDON(L)  
        WQKK(L) = 1.0 / (1.0 - WQF13)
        WQA13C=(WQFNDC*WQBMC(L)+WQFNDP*WQPRC(L))*WQANCC*WQO(L,1)  
        WQA13D=(WQFNDD*WQBMD(L)+WQFNDP*WQPRD(L))*WQANCD*WQO(L,2)  
        WQA13G=(WQFNDG*WQBMG(L)+WQFNDP*WQPRG(L))*WQANCG*WQO(L,3)  
        WQA13 = WQA13C+WQA13D+WQA13G  
        IF(IDNOTRVA > 0.AND.K == 1)THEN  
          WQA13 =WQA13 + (WQFNDM*WQBMM(L)+WQFNDPM*WQPRM(L)) *WQANCM*WQVO(L,K,IDNOTRVA)
        ENDIF  
        ! ***    PT_SRC_LOADS    VOLUME
        WQR13 = WQWPSL(L,K,13) * VOLWQ(L)

        WQRR(L) = WQV(L,K,13) + DTWQ*WQR13 + WQF13*WQV(L,K,13) + DTWQO2*( WQA13 + WQKLPN(L)*WQO(L,12))
      ENDDO
      IF(K == KC)THEN
        DO L=2,LA  
          ! ***     ATM DRY DEP     ATM WET DEP    VOLUME
          WQR13 = (WQWDSL(L,KC,13)+WQATML(L,KC,13))*VOLWQ(L)  
          WQRR(L) = WQRR(L) + DTWQ*WQR13  
        ENDDO
      ENDIF  
      DO L=2,LA  
        WQV(L,K,13)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,13)  
        WQO(L,13)=WQVO(L,K,13)+WQV(L,K,13)
      ENDDO  
    ENDIF  
! ****  PARAM 14  NHX - ammonia nitrogen
    IF(ISTRWQ(14) == 1)THEN  
      DO L=2,LA  
        ! ***      PT_SRC_LOADS    VOLUME
        WQRR(L) = WQWPSL(L,K,14) * VOLWQ(L)  
      ENDDO  
      IF(K == 1)THEN  
        DO L=2,LA  
          IF(LMASKDRY(L))THEN 
            WQRR(L) = WQRR(L) + WQBFNH4(L)*DZWQ(L)   ! *** Add in Benthic Flux
          ENDIF
        ENDDO  
      ENDIF  
      IF(K == KC)THEN
        DO L=2,LA  
          ! ***     ATM DRY DEP     ATM WET DEP    VOLUME
          WQR14 = (WQWDSL(L,KC,14)+WQATML(L,KC,14))*VOLWQ(L)  
          WQRR(L) = WQRR(L) + WQR14  
        ENDDO
      ENDIF  
      DO L=2,LA
        WQF14 = - DTWQO2*WQNIT(L)  
        WQKK(L) = 1.0 / (1.0 - WQF14) 
        WQA14C=WQFNIC*WQBMC(L)+WQFNIP*WQPRC(L)-WQPNC(L)*WQPC(L)  
        WQA14D=WQFNID*WQBMD(L)+WQFNIP*WQPRD(L)-WQPND(L)*WQPD(L)  
        WQA14G=WQFNIG*WQBMG(L)+WQFNIP*WQPRG(L)-WQPNG(L)*WQPG(L)  
        WQA14 = WQA14C*WQANCC*WQO(L,1) + WQA14D*WQANCD*WQO(L,2) + WQA14G*WQANCG*WQO(L,3)
        IF(IDNOTRVA > 0.AND.K == 1)THEN  
          WQA14 = WQA14 + (WQFNIM*WQBMM(L)+WQFNIPM*WQPRM(L) - WQPNM(L)*WQPM(L))*WQANCM*WQVO(L,K,IDNOTRVA)
        ENDIF  
        WQRR(L) = WQV(L,K,14) + DTWQ*WQRR(L) + WQF14*WQV(L,K,14) + DTWQO2*( WQA14 + WQKRPN(L)*WQO(L,11) + WQKDON(L)*WQO(L,13) )
        WQV(L,K,14)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,14)  
        WQO(L,14)=WQVO(L,K,14)+WQV(L,K,14)
      ENDDO  
    ENDIF  

! ****  PARAM 15  NOX - nitrate nitrogen
    IF(ISTRWQ(15) == 1)THEN  
      DO L=2,LA  
        ! ***      PT_SRC_LOADS    VOLUME
        WQRR(L) = WQWPSL(L,K,15) * VOLWQ(L)  
      ENDDO  
      IF(K == 1)THEN  
        DO L=2,LA  
          IF(LMASKDRY(L))THEN 
            WQRR(L) = WQRR(L) + WQBFNO3(L)*DZWQ(L)   ! *** Add in Benthic Flux
          ENDIF
        ENDDO  
      ENDIF  
      IF(K == KC)THEN
        DO L=2,LA  
          ! ***     ATM DRY DEP     ATM WET DEP    VOLUME  
          WQR15 = (WQWDSL(L,KC,15)+WQATML(L,KC,15))*VOLWQ(L)  
          WQRR(L) = WQRR(L) + WQR15  
        ENDDO
      ENDIF  
      DO L=2,LA  
        WQA15C = (WQPNC(L)-1.0)*WQPC(L) * WQANCC * WQO(L,1)  
        WQA15D = (WQPND(L)-1.0)*WQPD(L) * WQANCD * WQO(L,2)  
        WQA15G = (WQPNG(L)-1.0)*WQPG(L) * WQANCG * WQO(L,3)  
        WQA15 = WQA15C+WQA15D+WQA15G  
        IF(IDNOTRVA > 0.AND.K == 1)THEN  
          WQA15 =WQA15 + (WQPNM(L)-1.0)*WQPM(L)*WQANCM*WQVO(L,K,IDNOTRVA)
        ENDIF  
        WQB15 = WQV(L,K,15) + DTWQ*WQRR(L) + DTWQO2*( WQA15 - WQANDC*WQDENIT(L)*WQO(L,6) + WQNIT(L)*WQO(L,14))
        WQV(L,K,15)=SCB(L)*WQB15 + (1.-SCB(L))*WQV(L,K,15)  
        WQO(L,15)=WQVO(L,K,15)+WQV(L,K,15)
      ENDDO  
    ENDIF  
! ****  PARAM 16  SUU - particulate biogenic silica
    IF(ISTRWQ(16) == 1)THEN  
      IF(IWQSI == 1)THEN  
        DO L=2,LA  
          WQM16 = - (WQKSUA(IWQT(L)) + WQBDSET(L,1)) * DTWQO2
          WQKK(L) = 1.0 / (1.0 - WQM16)  
          WQA16D = (WQFSPD*WQBMD(L) + WQFSPP*WQPRD(L)) * WQASCD * WQO(L,2)  
          ! ***    PT_SRC_LOADS    VOLUME
          WQR16 = WQWPSL(L,K,16) * VOLWQ(L)
          WQRR(L) = WQV(L,K,16) + DTWQ*WQR16 + DTWQO2*WQA16D + WQM16*WQV(L,K,16)
        ENDDO  
        IF(K /= KC)THEN  
        ! *** Add in settling from above
          DO L=2,LA  
            WQRR(L) = WQRR(L) + DTWQO2*WQBDSET(L,2)*WQO(L,16)     ! *** PMC
          ENDDO  
        ELSE
          DO L=2,LA  
            ! ***     ATM DRY DEP     ATM WET DEP    VOLUME
            WQR16 = (WQWDSL(L,KC,16)+WQATML(L,KC,16))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR16  
          ENDDO
        ENDIF  
        DO L=2,LA  
          WQV(L,K,16)=SCB(L)*( WQRR(L)*WQKK(L) ) + (1.-SCB(L))*WQV(L,K,16)
          WQO(L,16)=WQVO(L,K,16)+WQV(L,K,16)
        ENDDO  
      ENDIF  
    ENDIF  
! ****  PARAM 17  SAA - dissolved available silica
    IF(ISTRWQ(17) == 1)THEN  
      IF(IWQSI == 1)THEN  
        DO L=2,LA  
          WQKK(L) = (WQFSID*WQBMD(L) + WQFSIP*WQPRD(L) - WQPD(L)) * WQASCD * WQO(L,2)
          ! ***      PT_SRC_LOADS    VOLUME  
          WQRR(L) = WQWPSL(L,K,17) * VOLWQ(L)  
          ! ***                   ATM WET DEP      VOLUME  
          WQRR(L) = WQRR(L) + DELKC(K)*(WQATML(L,KC,17)*VOLWQ(L))  
        ENDDO  
        IF(K == 1)THEN  
          DO L=2,LA
            IF(LMASKDRY(L))THEN 
              WQRR(L) = WQRR(L) + WQBFSAD(L)*DZWQ(L)   ! *** Add in Benthic Flux
            ENDIF
          ENDDO  
        ENDIF  
        DO L=2,LA  
          WQRR(L) = WQV(L,K,17) + DTWQ*WQRR(L) +WQN17(L)*WQV(L,K,17) + DTWQO2*( WQKK(L) + WQKSUA(IWQT(L))*WQO(L,16))
        ENDDO  
        IF(K /= KC)THEN  
          DO L=2,LA  
            WQRR(L) = WQRR(L) + DTWQO2*WQT17(L)*WQO(L,17)  
          ENDDO  
        ELSE
          DO L=2,LA  
            ! ***     ATM DRY DEP     ATM WET DEP    VOLUME
            WQR17 = (WQWDSL(L,KC,17)+WQATML(L,KC,17))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR17  
          ENDDO
        ENDIF  
        DO L=2,LA  
          WQKK(L) = 1.0 / (1.0 - WQN17(L))
          WQV(L,K,17)=SCB(L)*( WQRR(L)*WQKK(L) ) + (1.-SCB(L))*WQV(L,K,17)
          WQO(L,17)=WQVO(L,K,17)+WQV(L,K,17)
        ENDDO  
      ENDIF  
    ENDIF  
! ****  PARAM 18  COD - chemical oxygen demand
    IF(ISTRWQ(18) == 1)THEN  
      DO L=2,LA  
        WQKK(L) = 1.0 / (1.0 - WQO18(L))  
          ! ***    PT_SRC_LOADS    VOLUME
        WQRR(L) = WQWPSL(L,K,18) * VOLWQ(L)  
      ENDDO  
      IF(K == 1)THEN  
        DO L=2,LA  
          IF(LMASKDRY(L))THEN 
            WQRR(L) = WQRR(L) + WQBFCOD(L)*DZWQ(L)   ! *** Add in Benthic Flux
          ENDIF
        ENDDO  
      ENDIF  
      IF(K == KC)THEN
        DO L=2,LA  
          ! ***     ATM DRY DEP     ATM WET DEP    VOLUME
          WQR18 = (WQWDSL(L,KC,18)+WQATML(L,KC,18))*VOLWQ(L)  
          WQRR(L) = WQRR(L) + WQR18  
        ENDDO
      ENDIF  
      DO L=2,LA  
        WQRR(L)=WQV(L,K,18)+DTWQ*WQRR(L)+WQO18(L)*WQV(L,K,18)  
        WQV(L,K,18)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,18)  
        WQO(L,18) = WQVO(L,K,18)+WQV(L,K,18)  
      ENDDO  
    ENDIF
! ****  PARAM 19  DOX - dissolved oxygen
    ! ***  1) CHC - cyanobacteria
    ! ***  2) CHD - diatom algae
    ! ***  3) CHG - green algae
    ! ***  4) ROC - refractory particulate organic carbon
    ! ***  5) LOC - labile particulate organic carbon
    ! ***  6) DOC - dissolved organic carbon
    ! ***  7) ROP - refractory particulate organic phosphorus
    ! ***  8) LOP - labile particulate organic phosphorus
    ! ***  9) DOP - dissolved organic phosphorus
    ! *** 10) P4D - total phosphate
    ! *** 11) RON - refractory particulate organic nitrogen 
    ! *** 12) LON - labile particulate organic nitrogen
    ! *** 13) DON - dissolved organic nitrogen
    ! *** 14) NHX - ammonia nitrogen
    ! *** 15) NOX - nitrate nitrogen
    ! *** 16) SUU - particulate biogenic silica
    ! *** 17) SAA - dissolved available silica
    ! *** 18) COD - chemical oxygen demand
    ! *** 19) DOX - dissolved oxygen
    ! *** 20) TAM - total active metal
    ! *** 21) FCB - fecal coliform bacteria
    ! *** 22) CO2 - dissolved carbon dioxide
    ! *** 23) MAC - macroalgae

    !C 04/29/99 MRM:  
    !C THE FOLLOWING ARRAYS WERE ADDED TO KEEP TRACK OF THE VARIOUS COMPONENT  
    !C OF DISSOLVED OXYGEN.  THE INSTANTANEOUS VALUES FOR EACH COMPONENT ARE  
    !C SUMMED IN THE ARRAYS AND THEN DUMPED TO THE WQDOCOMP.BIN FILE AT THE  
    !C SAME TIME INTERVAL AS FOR THE WQWCAVG.BIN FILES (I.E., IWQTSDT INTERVA  
    !C USUALLY DAILY AVERAGES).  THE ARRAY DESCRIPTIONS ARE:  
    !C  XDOPSL(L,L) = D.O. COMPONENT FOR POINT SOURCE LOADS
    !C  XDOSOD(L,K) = D.O. COMPONENT FOR SEDIMENT OXYGEN DEMAND  
    !C  XDOKAR(L,K) = D.O. COMPONENT FOR REAERATION  
    !C  XDODOC(L,K) = D.O. COMPONENT FOR DISS. ORG. CARBON DECAY  
    !C  XDONIT(L,K) = D.O. COMPONENT FOR AMMONIA NITRIFICATION  
    !C  XDOCOD(L,K) = D.O. COMPONENT FOR CHEM. OXY. DEMAND OXIDATION  
    !C  XDOPPB(L,K) = D.O. COMPONENT FOR PHOTOSYNTHESIS OF TOTAL CHLOROPHYLL  
    !C  XDORRB(L,K) = D.O. COMPONENT FOR RESPIRATION OF TOTAL CHLOROPHYLL  
    !C  XDOPPM(L,K) = D.O. COMPONENT FOR PHOTOSYNTHESIS OF MACROALGAE  
    !C  XDORRM(L,K) = D.O. COMPONENT FOR RESPIRATION OF MACROALGAE  
    !C  XDOALL(L,K) = SUM OF THE ABOVE 10 D.O. COMPONENTS  
    !C  NLIM = COUNTER FOR NUMBER OF ITEMS SUMMED IN EACH ARRAY SLOT  

    IF(ISTRWQ(19) == 1)THEN 
      ! *** Point Source Loading 
      DO L=2,LA  
        WQRR(L) = WQWPSL(L,K,19) * VOLWQ(L)  
      ENDDO
      
      IF(ISDIURDO > 0)THEN
        DO L=2,LA  
          TMP19=WQRR(L)*DTWQ*DZCHP(L)
          XDOPSL(L,K) = XDOPSL(L,K) + TMP19
          XDOALL(L,K) = XDOALL(L,K) + TMP19  
        ENDDO  
      ENDIF
      
      ! *** Handle Surface Processes
      IF(K == KC)THEN  
        DO L=2,LA  
          WQKK(L) = 1.0 / (1.0 - DTWQO2*WQP19(L)) 
          ! ***             ATM DRY DEP    ATM WET DEP    VOLUME  
          WQRR(L)=WQRR(L)+(WQWDSL(L,KC,19)+WQATML(L,KC,19))*VOLWQ(L)
          ! *** Reaeration - Atm to water flux (Offset later in O2 Mass Balance by WQRea term)
          WQRR(L) = WQRR(L) + WQKRDOS(L)
        ENDDO
      ELSE  
        DO L=2,LA  
          WQKK(L) = 1.0
        ENDDO
      ENDIF 
      ! *** Bottom Processes 
      IF(K == 1)THEN  
        DO L=2,LA  
          IF(LMASKDRY(L))THEN 
            TEMFAC=1.065**(TEM(L,1)-20.)  
            WQRR(L) = WQRR(L) + TEMFAC*WQBFO2(L)*DZWQ(L)   ! *** Add in Benthic Flux
          ENDIF
        ENDDO
          
        IF(ISDIURDO > 0)THEN
          DO L=2,LA  
            IF(LMASKDRY(L))THEN 
              TEMFAC=1.065**(TEM(L,1)-20.)  
              TMP19=TEMFAC*WQBFO2(L)*DTWQ
              XDOSOD(L,K) = XDOSOD(L,K) + TMP19  
              XDOALL(L,K) = XDOALL(L,K) + TMP19
            ENDIF
          ENDDO  
        ENDIF
      ENDIF  
  
      DO L=2,LA  
        DTWQxH = DTWQ*DZCHP(L)  
        DTWQxH2= DTWQO2*DZCHP(L)

        IF(WQI0  <=  0.001)THEN  
          WQTTC = 0.0
          WQTTD = 0.0
          WQTTG = 0.0
        ELSE
          WQTTC = (1.3 - 0.3*WQPNC(L)) * WQPC(L) 
          WQTTD = (1.3 - 0.3*WQPND(L)) * WQPD(L)  
          WQTTG = (1.3 - 0.3*WQPNG(L)) * WQPG(L)
          ! *** PHOTOSYNTHESIS OF TOTAL CHLOROPHYLL  
          TMP19 = WQAOCR*DTWQxH2*(WQTTC*WQO(L,1)+WQTTD*WQO(L,2)+WQTTG*WQO(L,3))
          XDOPPB(L,K) = XDOPPB(L,K) + TMP19
          XDOALL(L,K) = XDOALL(L,K) + TMP19
        ENDIF
        ! *** RESPIRATION OF TOTAL CHLOROPHYLL - CYANOBACTERIA 
        XMRM = CFCDCWQ*O2WQ(L)*WQBMC(L)/(WQKHRC+O2WQ(L)+ 1.E-18)  
        WQA19C = WQTTC - XMRM
        TMP19  = XMRM*WQO(L,1)*WQAOCR * DTWQxH2  
        XDORRB(L,K) = XDORRB(L,K) - TMP19 
        XDOALL(L,K) = XDOALL(L,K) - TMP19
        ! *** RESPIRATION OF TOTAL CHLOROPHYLL - DIATOMS 
        XMRM = CFCDDWQ*O2WQ(L)*WQBMD(L)/(WQKHRD+O2WQ(L)+ 1.E-18)  
        WQA19D = WQTTD - XMRM
        TMP19  = XMRM*WQO(L,2)*WQAOCR * DTWQxH2
        XDORRB(L,K) = XDORRB(L,K) - TMP19  
        XDOALL(L,K) = XDOALL(L,K) - TMP19
        ! *** RESPIRATION OF TOTAL CHLOROPHYLL -  GREENS
        XMRM = CFCDGWQ*O2WQ(L)*WQBMG(L)/(WQKHRG+O2WQ(L)+ 1.E-18)  
        WQA19G = WQTTG - XMRM
        TMP19  = XMRM*WQO(L,3)*WQAOCR * DTWQxH2   
        XDORRB(L,K) = XDORRB(L,K) - TMP19
        XDOALL(L,K) = XDOALL(L,K) - TMP19
        ! *** TOTAL NET RESPIRATION/PHOTOSYNTHESIS
        WQA19=(WQA19C*WQO(L,1) + WQA19D*WQO(L,2) + WQA19G*WQO(L,3)) * WQAOCR
        ! *** MODIFIED BY MRM 05/23/99 TO ALLOW DIFFERENT AOCR CONSTANTS TO BE APPLIED  
        ! ***   TO PHOTOSYNTHESIS AND RESPIRATION TERMS FOR MACROALGAE  
        IF(IDNOTRVA > 0.AND.K == 1)THEN  
          ! *** TRAPEZOIDAL AVERAGE CONCENTRATIONS
          WQO(L,IDNOTRVA)=WQVO(L,K,IDNOTRVA)+WQV(L,K,IDNOTRVA)
          IZ = IWQZMAP(L,K)  
          WQTTM = (1.3 - 0.3*WQPNM(L)) * WQPM(L)  
          XMRM=(1.0-WQFCDM)*O2WQ(L)*WQBMM(L)/(WQKHRM(IZ)+O2WQ(L)+1.E-18)
          WQA19A = WQTTM*WQO(L,IDNOTRVA)*WQAOCRPM - XMRM*WQO(L,IDNOTRVA)*WQAOCRRM  
          WQA19 = WQA19 + WQA19A
          TMP19 = WQTTM*WQO(L,IDNOTRVA)*WQAOCRPM * DTWQxH2     
          XDOPPM(L,K) = XDOPPM(L,K) + TMP19
          XDOALL(L,K) = XDOALL(L,K) + TMP19 
          TMP19 = XMRM*WQO(L,IDNOTRVA)*WQAOCRRM * DTWQxH2
          XDORRM(L,K) = XDORRM(L,K) - TMP19
          XDOALL(L,K) = XDOALL(L,K) - TMP19
        ENDIF  
        ! *** O2 Mass Balance
        ! WQA19                         ! *** Total Net Respiration/Photosynthesis
        WQSUM=DTWQ*WQRR(L)              ! *** Sum of Loadings/Demands
        WQRea=WQP19(L)*WQV(L,K,19)      ! *** Reaeration - Offsetting Flux (WQP19 is negative)
        WQPOC=WQAOCR*WQKRPC(L)*WQO(L,4) ! *** POC
        WQDOC=WQAOCR*WQKHR(L) *WQO(L,6) ! *** DOC
        WQNH3=WQAONT*WQNIT(L) *WQO(L,14)! *** Ammonia
        WQCOD=WQO18(L)*WQO(L,18)        ! *** COD
        WQRR(L) = WQV(L,K,19) + WQSUM  + WQCOD  + DTWQO2*(WQA19 - WQPOC - WQDOC - WQNH3 + WQRea)
        WQV(L,K,19)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,19)  ! *** remove the SCB computations!
        ! *** WQV(L,K,19) After this WQV(L,K,19) can not be < 0.
        WQV(L,K,19) = MAX(WQV(L,K,19), 0.0)  
  
        ! *** COMPUTE AND SAVE D.O. DEFICIT:  
        IF(ISDIURDO > 0)THEN
          WQO(L,19)=WQVO(L,K,19)+WQV(L,K,19)
          XMRM = WQDOS(L) - WQV(L,K,19)  
          XDODEF(L,K) = XDODEF(L,K) + XMRM*DTWQ*DZCHP(L)  
          IF(K == KC)THEN
            TMP19=WQKRDOS(L)*DTWQ*DZCHP(L) + WQP19(L)*WQO(L,19)*DTWQxH2  
            XDOKAR(L,K) = XDOKAR(L,K) + TMP19  
            XDOALL(L,K) = XDOALL(L,K) + TMP19  
          ENDIF

          TMP19=WQAOCR*WQKHR(L)*WQO(L,6)*DTWQxH2  
          XDODOC(L,K)=XDODOC(L,K) - TMP19  
          XDOALL(L,K)=XDOALL(L,K) - TMP19  

          TMP19=WQAONT*WQNIT(L)*WQO(L,14)*DTWQxH2  
          XDONIT(L,K)=XDONIT(L,K) - TMP19
          XDOALL(L,K)=XDOALL(L,K) - TMP19

          TMP19=WQO18(L)*WQO(L,18)*DZCHP(L)
          XDOCOD(L,K)=XDOCOD(L,K) - TMP19  
          XDOALL(L,K)=XDOALL(L,K) - TMP19  

          XDODZ(L,K) = XDODZ(L,K) + DZCHP(L)  
        ENDIF
! O2-reaeration changes
! If DO conc. is greater than saturated DO conc., excess DO
! gets released out, which is calculated as XDOOUT.
        IF(K==KC)THEN
	  IF(WQV(L,KC,19)>WQDOS(L))THEN
	    XDOOUT(L,KC)=(WQV(L,KC,19)-WQDOS(L))*DTWQ*DZCHP(L)
	    WQV(L,KC,19)=WQDOS(L)
	    XDOALL(L,K)=XDOALL(L,K)-XDOOUT(L,KC)
	  ENDIF
	ENDIF

      ENDDO  
      
    ENDIF  ! *** END OF ISTRWQ(19)

    ! ****  PARAM 20  TAM - total active metal
    IF(ISTRWQ(20) == 1)THEN  
      IF(IWQSRP == 1)THEN  
        DO L=2,LA  
          WQT20 = - DTWQ*WQWSSET(L,1)    ! *** DTWQO2
          WQKK(L) = 1.0 / (1.0 - WQT20)  
          WQRR(L)=WQV(L,K,20)+DTWQ*WQR20(L)+WQT20*WQV(L,K,20)  
        ENDDO  
        IF(K /= KC)THEN  
          ! *** Add in settling from above
          DO L=2,LA  
            WQRR(L) = WQRR(L) + DTWQO2*WQWSSET(L,2)*WQO(L,20)
          ENDDO  
        ENDIF  
        DO L=2,LA  
          WQV(L,K,20)=SCB(L)*( WQRR(L)*WQKK(L) ) + (1.-SCB(L))*WQV(L,K,20)  
          WQO(L,20)=WQVO(L,K,20)+WQV(L,K,20)
        ENDDO  
      ENDIF  
    ENDIF  

    ! ****  PARAM 21  FCB - fecal coliform bacteria
    IF(ISTRWQ(21) == 1)THEN  
      IF(IWQFCB == 1)THEN  
        DO L=2,LA  
          WQKK(L) = WQTD2FCB(IWQT(L))  
    
          ! ***      ATM DRY DEP       LOADS          VOLUME  
          WQR21= (WQWDSL(L,K,21)+WQWPSL(L,K,21))*VOLWQ(L)      !VB CHANGED NWQV TO 21
          ! WQR21= (WQWDSL(L,K,NWQV)+WQWPSL(L,K,NWQV))*VOLWQ(L)  
    
          ! ***                   ATM WET DEP      VOLUME  
          WQR21 = WQR21 + DELKC(K)*(WQATML(L,KC,21)*VOLWQ(L))  
          ! WQRR(L) = WQV(L,K,NWQV)*WQTD1FCB(IWQT(L)) + DTWQ*WQR21    !VB CHANGED NWQV TO 21
          WQRR(L) = WQV(L,K,21)*WQTD1FCB(IWQT(L)) + DTWQ*WQR21
          WQV(L,K,21)=SCB(L)*( WQRR(L)*WQKK(L) ) + (1.-SCB(L))*WQV(L,K,21)  
        ENDDO  
      ENDIF  
    ENDIF

    ! ****  PARAM 22 DISSOLVED CARBON DIOXIDE

    !C THE FOLLOWING ARRAYS WERE ADDED TO KEEP TRACK OF THE VARIOUS COMPONENT  
    !C OF DISSOLVED CARBON DIOXIDE.  
    !C THE ARRAY DESCRIPTIONS ARE:  
    !C  XCDOKAR(L,K) = CDO. COMPONENT FOR REAERATION  
    !C  XCDODOC(L,K) = CDO. COMPONENT FOR DISS. ORG. CARBON DECAY  
    !C  XCDOPPB(L,K) = CDO. COMPONENT FOR PHOTOSYNTHESIS OF TOTAL CHLOROPHYLL  
    !C  XCDORRB(L,K) = CDO. COMPONENT FOR RESPIRATION OF TOTAL CHLOROPHYLL  
    !C  XCDOPPM(L,K) = CDO. COMPONENT FOR PHOTOSYNTHESIS OF MACROALGAE  
    !C  XCDORRM(L,K) = CDO. COMPONENT FOR RESPIRATION OF MACROALGAE  
    !C  XCDOALL(L,K) = SUM OF THE ABOVE 6 CDO. COMPONENTS  

!Checking if WQ variable 22 is enabled
    IF(ISTRWQ(22) == 1)THEN
!------------------------------------------------------------- 
!		looping over all horizontal cells
!		WQRR 	- biomass concentration (g / (day*m3))
!		WQWPSL 	- Point source load for each layer (g / day)
!		VOLWQ	- inverse volume of each cell in a layer ( / m3)
!		DTWQ	- WQ time step
!		DZCHP	- layer thickness of a cell in meters = DZC (vertical layer thickness as decimal fraction of water depth) * HP (water depth at cell center)
!		In this loop, the biomass concentration WQRR at each horizontal cell is calculated.
!-------------------------------------------------------------          
      DO L=2,LA  
        WQRR(L) = WQWPSL(L,K,22) * VOLWQ(L)  
      ENDDO
!-------------------------------------------------------------
! 		Handle Surface Processes
!		K = KC implies the surface layer.
!		DTWQO2	- Half time step
!		WQP22	- Production rate per unit time
!		WQKK	- has the form 1 / (1 - production per half time step) 
!		(WQWDSL(L,KC,22)+WQATML(L,KC,22))*VOLWQ(L)	- atmospheric loading mass per time / cell volume            
!-------------------------------------------------------------
      IF(K == KC)THEN  
        DO L=2,LA  
          WQKK(L) = 1.0 / (1.0 - DTWQO2*WQP22(L)) 
          ! ***             ATM DRY DEP    ATM WET DEP    VOLUME  
          WQRR(L)=WQRR(L)+(WQWDSL(L,KC,22)+WQATML(L,KC,22))*VOLWQ(L)
          ! *** Reaeration
          WQRR(L) = WQRR(L) + WQKRCDOS(L)  
        ENDDO
      ELSE  
        DO L=2,LA  
          WQKK(L) = 1.0
        ENDDO
      ENDIF
!-------------------------------------------------------------
!		WQI0	- Solar radiation for current time
!		WQPNC	- Preference for ammonium uptake for cyanobacteria
!		WQPC	- Final net growth rate	for cyanobacteria
!		WQTTC	- Photosynthesis ratio - Eq 2.7 WQ manual and page 47.
!		Here if solar radiation is less than 0.001, then no photosynthesis happens.
!-------------------------------------------------------------
      DO L=2,LA  
        DTWQxH = DTWQ*DZCHP(L)  
        DTWQxH2= DTWQO2*DZCHP(L)
        IF(WQI0  <=  0.001)THEN  
          WQTTC = 0.0
          WQTTD = 0.0
          WQTTG = 0.0
        ELSE
          ! *** PHOTOSYNTHESIS OF TOTAL CHLOROPHYLL  
          WQTTC = (1.3 - 0.3*WQPNC(L)) * WQPC(L) 
          WQTTD = (1.3 - 0.3*WQPND(L)) * WQPD(L)  
          WQTTG = (1.3 - 0.3*WQPNG(L)) * WQPG(L)
        ENDIF

!------------------------------------------------------------- 
!		WQA22C	- Net growth/production of cyanobacteria = photosynthesis - respiration
!		WQBMC	- Current Basal Metabolism Rate for cyanobacteria
!		WQKHRC	- DOC Heterotrophic Respiration Rate for cyanobacteria
!		CFCDCWQ	- Complementary FCD = 1 - FCD, 0 < FCD < 1, Eq 2.33 - 2.37 WQ manual
!		O2WQ	- Dissolved oxygen concentration (g O2 m-3)
!		XMRM	- Respiration
!		WQO		- Depth totaled biomass conc
!		WQV		- Transported biomass conc from current iteration, WQVO - from old iteration
!		WQDENIT	- Denitrification rate
!		SCB		- Boundary cells (0 or 1)
!-------------------------------------------------------------
        ! *** RESPIRATION OF TOTAL CHLOROPHYLL - CYANOBACTERIA 
        XMRM = CFCDCWQ*O2WQ(L)*WQBMC(L)/(WQKHRC+O2WQ(L)+ 1.E-18)  
        WQA22C = WQTTC - XMRM

        ! *** RESPIRATION OF TOTAL CHLOROPHYLL - DIATOMS 
        XMRM = CFCDDWQ*O2WQ(L)*WQBMD(L)/(WQKHRD+O2WQ(L)+ 1.E-18)  
        WQA22D = WQTTD - XMRM

        ! *** RESPIRATION OF TOTAL CHLOROPHYLL -  GREENS
        XMRM = CFCDGWQ*O2WQ(L)*WQBMG(L)/(WQKHRG+O2WQ(L)+ 1.E-18)  
        WQA22G = WQTTG - XMRM
        
        ! *** TOTAL NET RESPIRATION/PHOTOSYNTHESIS
        WQA22=3.67*(WQA22C*WQO(L,1)+WQA22D*WQO(L,2)+WQA22G*WQO(L,3))  !VB 3.67 CONVERTS g CARBON TO g CO2

        ! *** MODIFIED BY MRM 05/23/99 TO ALLOW DIFFERENT AOCR CONSTANTS TO BE APPLIED  
        ! ***   TO PHOTOSYNTHESIS AND RESPIRATION TERMS FOR MACROALGAE  
        ! *** TRAPEZOIDAL AVERAGE CONCENTRATIONS
        ! *** CO2 Mass Balance
        ! WQA22                                         ! *** Total Net Respiration/Photosynthesis
        WQCDSUM=DTWQ*WQRR(L)                            ! *** Sum of Loadings/Demands
        WQCDRea=WQP22(L)*WQV(L,K,22)                    ! *** Reaeration
        WQCDDOC=(WQKHR(L)+WQDENIT(L))*WQO(L,6)*3.67      ! *** DOC FROM HYDROLYSIS AND DENITRIFICATION 3.67 CONVERTS G CARBON TO G CO2    

        WQRR(L) = WQV(L,K,22) + WQCDSUM + DTWQO2*(-WQA22 + WQCDDOC + WQCDRea)
        WQV(L,K,22)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,22)  ! *** remove the SCB computations!
        ! *** WQV(L,K,22) After this WQV(L,K,22) can not be < 0.
        WQV(L,K,22) = MAX(WQV(L,K,22), 0.0)  !this is a bit dangerous because it can cover up problems that arise from too large of a time step
      ENDDO  
    ENDIF  

  ENDDO  ! *** END OF THE KC LOOP
  
  ! ***************************************************************************
  ! ----------------------------------------------------------------  
  !  
  ! INCREMENT COUNTER FOR LIMITATION AND XDOXXX DO COMPONENT ARRAYS:  
  !  
  IF(ISDYNSTP == 0)THEN  
    TIMTMP=DT*FLOAT(N)+TCON*TBEGIN  
    TIMTMP=TIMTMP/TCTMSR  
  ELSE  
    TIMTMP=TIMESEC/TCTMSR  
  ENDIF  
  TIMESUM3 = TIMESUM3 + TIMTMP  
  NLIM = NLIM + 1  

  ! PMC - Moved CHLa Computations to the beginning of the WQ Calculations
  IF(IWQSRP == 1)THEN  
    ! *** Sorption Option: TAM
    DO K=1,KC  
      DO L=2,LA  
        O2WQ(L) = MAX(WQV(L,K,19), 0.0)  
        WQTAMD = MIN( WQTAMDMX*EXP(-WQKDOTAM*O2WQ(L)), WQV(L,K,20) )  
        WQTAMP(L,K) = WQV(L,K,20) - WQTAMD  
        WQPO4D(L,K) = WQV(L,K,10) / (1.0 + WQKPO4P*WQTAMP(L,K))  
        WQSAD(L,K)  = WQV(L,K,17) / (1.0 + WQKSAP*WQTAMP(L,K))  
      ENDDO  
    ENDDO  
  ELSE IF(IWQSRP == 2)THEN
    ! *** Sorption Option: Sediments
    DO K=1,KC  
      DO L=2,LA  
        WQPO4D(L,K) = WQV(L,K,10) / (1.0 + WQKPO4P*SEDT(L,K))  
        WQSAD(L,K)  = WQV(L,K,17) / (1.0 + WQKSAP*SEDT(L,K))  
      ENDDO  
    ENDDO  
  ELSE  
    DO K=1,KC  
      DO L=2,LA  
        WQPO4D(L,K) = WQV(L,K,10)  
        WQSAD(L,K)  = WQV(L,K,17)  
      ENDDO  
    ENDDO  
  ENDIF  

  ! ***************************************************************************
  ! COUPLING TO SEDIMENT MODEL  
  ! EVALUATE DEP. FLUX USING NEW VALUES CAUSE IMPLICIT SCHEME IS USED IN SPM
  IF(IWQBEN == 1)THEN  
    DO L=2,LA  
      IMWQZ = IWQZMAP(L,1)  
      WQDFBC(L) = SCB(L)*WQWSC(IMWQZ)*WQV(L,1,1)  
      WQDFBD(L) = SCB(L)*WQWSD(IMWQZ)*WQV(L,1,2)  
      WQDFBG(L) = SCB(L)*WQWSG(IMWQZ)*WQV(L,1,3)+WQWSM*DZWQ(L)*WQV(L,1,IDNOTRVA)  
      WQDFRC(L) = SCB(L)*WQWSRP(IMWQZ)*WQV(L,1,4)  
      WQDFLC(L) = SCB(L)*WQWSLP(IMWQZ)*WQV(L,1,5)  
      WQDFRP(L) = SCB(L)*WQWSRP(IMWQZ)*WQV(L,1,7)  
      WQDFLP(L) = SCB(L)*WQWSLP(IMWQZ)*WQV(L,1,8)  
      WQDFRN(L) = SCB(L)*WQWSRP(IMWQZ)*WQV(L,1,11)  
      WQDFLN(L) = SCB(L)*WQWSLP(IMWQZ)*WQV(L,1,12)  
      IF(IWQSI == 1) WQDFSI(L) = SCB(L)*WQWSD(IMWQZ)*WQV(L,1,16)  
    ENDDO  
    IF(IWQSRP == 1)THEN  
      DO L=2,LA  
        IMWQZ = IWQZMAP(L,1)  
        WQDFLP(L) = SCB(L)*( WQDFLP(L) + WQWSS(IMWQZ)*( WQV(L,1,10)-WQPO4D(L,1) ) )  
        IF(IWQSI == 1) WQDFSI(L) = SCB(L)*( WQDFSI(L) + WQWSS(IMWQZ)*( WQV(L,1,17)-WQSAD(L,1) ) )  
      ENDDO  
    ELSE IF(IWQSRP == 2)THEN  
      DO L=2,LA  
        WQDFLP(L) = SCB(L)*( WQDFLP(L)+WSEDO(NS)*( WQV(L,1,10)-WQPO4D(L,1) ) )  
        IF(IWQSI == 1) WQDFSI(L) = SCB(L)*( WQDFSI(L) + WSEDO(NS)*( WQV(L,1,17)-WQSAD(L,1) ) )  
      ENDDO  
    ENDIF  
  ENDIF  

  ! ***************************************************************************
  ! DIURNAL DO ANALYSIS  
  IF(NDDOAVG >= 1.AND.DEBUG)THEN  
    OPEN(1,FILE='DIURNDO.OUT',POSITION='APPEND')  
    NDDOCNT=NDDOCNT+1  
    NSTPTMP=NDDOAVG*NTSPTC/2  
    RMULTMP=1./FLOAT(NSTPTMP)  
    DO K=1,KC  
      DO L=2,LA  
        DDOMAX(L,K)=MAX(DDOMAX(L,K),WQV(L,K,19))  
        DDOMIN(L,K)=MIN(DDOMIN(L,K),WQV(L,K,19))  
      ENDDO  
    ENDDO  
    IF(NDDOCNT == NSTPTMP)THEN  
      NDDOCNT=0  
      IF(ISDYNSTP == 0)THEN  
        TIME=DT*FLOAT(N)+TCON*TBEGIN  
        TIME=TIME/TCON  
      ELSE  
        TIME=TIMESEC/TCON  
      ENDIF  
      WRITE(1,1111)N,TIME  
      DO L=2,LA  
        WRITE(1,1112)IL(L),JL(L),(DDOMIN(L,K),K=1,KC),(DDOMAX(L,K),K=1,KC)  
      ENDDO  
      DO K=1,KC  
        DO L=2,LA  
          DDOMAX(L,K)=-1.E6  
          DDOMIN(L,K)=1.E6  
        ENDDO  
      ENDDO  
    ENDIF  
    CLOSE(1)  
  ENDIF  

  ! ***************************************************************************
  ! *** LIGHT EXTINCTION ANALYSIS  
  IF(NDLTAVG >= 1)THEN  
    OPEN(1,FILE='LIGHT.OUT',POSITION='APPEND')  
    NDLTCNT=NDLTCNT+1  
    NSTPTMP=NDLTAVG*NTSPTC/2  
    RMULTMP=1./FLOAT(NSTPTMP)  
    DO K=1,KC  
      DO L=2,LA  
        RLIGHT1=WQKEB(IMWQZT(L))+WQKETSS*SEDT(L,K)  
        XMRM = WQKECHL*WQCHL(L,K)  
        IF(WQKECHL  <  0.0)THEN  
          XMRM = 0.054*WQCHL(L,K)**0.6667 + 0.0088*WQCHL(L,K)  
        ENDIF  
        RLIGHT2 = XMRM  
        RLIGHTT(L,K)=RLIGHTT(L,K)+RLIGHT1  
        RLIGHTC(L,K)=RLIGHTC(L,K)+RLIGHT1+RLIGHT2  
      ENDDO  
    ENDDO  
    IF(NDLTCNT == NSTPTMP)THEN  
      NDLTCNT=0  
      IF(ISDYNSTP == 0)THEN  
        TIME=DT*FLOAT(N)+TCON*TBEGIN  
        TIME=TIME/TCON  
      ELSE  
        TIME=TIMESEC/TCON  
      ENDIF  
      DO K=1,KC  
        DO L=2,LA  
          RLIGHTT(L,K)=RMULTMP*RLIGHTT(L,K)  
          RLIGHTC(L,K)=RMULTMP*RLIGHTC(L,K)  
        ENDDO  
      ENDDO  
      WRITE(1,1111)N,TIME  
      DO L=2,LA  
        WRITE(1,1113)IL(L),JL(L),(RLIGHTT(L,K),K=1,KC),(RLIGHTC(L,K),K=1,KC)  
      ENDDO  
      DO K=1,KC  
        DO L=2,LA  
          RLIGHTT(L,K)=0.  
          RLIGHTC(L,K)=0.  
        ENDDO  
      ENDDO  
    ENDIF  
    CLOSE(1)  
  ENDIF  

  1111 FORMAT(I12,F10.4)  
  1112 FORMAT(2I5,12F7.2)  
  1113 FORMAT(2I5,12E12.4)  
  1414 FORMAT(I12,11E12.4)  

  RETURN  

END  
