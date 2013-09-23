SUBROUTINE CALHDMF
 
  ! *** CALDMF CALCULATES THE HORIZONTAL VISCOSITY AND                                                                    
  ! *** DIFFUSIVE MOMENTUM FLUXES. THE VISCOSITY, AH IS CALCULATED USING                                                  
  ! *** SMAGORINSKY'S SUBGRID SCALE FORMULATION PLUS A CONSTANT AHO                                                       
 
  ! *** ONLY VALID FOR ISHDMF >= 1
  !
  !----------------------------------------------------------------------C  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2011-05           Paul M. Craig      Corrected DSQR equation from /4 to /2
  ! 2011-03           Paul M. Craig      Rewritten to F90 and added OMP
  ! 2008-10           SANG YUK (DS-INTL)   CORRECTED THE DIFFUSIVE MOMENTUM FLUXES COMPUTATION                                     
  ! 2004-11           PAUL M. CRAIG      REWRITTEN AND RESTRUCTURED

  USE GLOBAL
  
  IMPLICIT NONE
  
  INTEGER::LWVMASK(LA) !SCJ added to make compatible with unknown EE upgrades
                                                                        
  INTEGER :: L,LW,K,LF,LL,NQSTMP,IU,JU,KU,IWR,ND,LN,LS
  REAL    :: SLIPCO,DY2DZBR,DX2DZBR,CSDRAG,SLIPFAC,TMPVAL,DSQR,WVFACT
  REAL    :: DTMPH,DTMPX,AHWVX,SXYLN,SXYEE
   
  REAL DYU1L, DXV1L, SXYL
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)::AHEE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)::AHNN
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)::SXY
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)::SXY2CC
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)::SXY2EE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)::SXY2NN
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)::DYU1
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)::DYV1
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)::DXU1
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)::DXV1
  
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:)::ICORDXV  
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:)::ICORDYU  
  
  LWVMASK(:)=1 !SCJ added to make compatible with unknown EE upgrades 
  
  IF(.NOT.ALLOCATED(AHEE))THEN
    ALLOCATE(AHEE(LCM,KCM))
    ALLOCATE(AHNN(LCM,KCM))
    ALLOCATE(SXY(LCM,KCM))
    ALLOCATE(SXY2CC(LCM,KCM))
    ALLOCATE(SXY2EE(LCM,KCM))
    ALLOCATE(SXY2NN(LCM,KCM))
    ALLOCATE(ICORDXV(LCM))  
    ALLOCATE(ICORDYU(LCM))  
    ALLOCATE(DYU1(LCM,KCM))  
    ALLOCATE(DYV1(LCM,KCM))  
    ALLOCATE(DXU1(LCM,KCM))  
    ALLOCATE(DXV1(LCM,KCM))  

    AHEE=0.0
    AHNN=0.0
    SXY=0.0
    SXY2CC=0.0
    SXY2EE=0.0
    SXY2NN=0.0
    ICORDXV=0
    ICORDYU=0
    DYU1=0.
    DYV1=0.
    DXU1=0.
    DXV1=0.
  ENDIF

  AHMAX=AHO

  SLIPCO=1.
  IF(AHD > 0.0)THEN
    SLIPCO=0.5/SQRT(AHD)
  ENDIF
!$OMP PARALLEL DO
  DO ND=1,NDM  
    LF=2+(ND-1)*LDM  
    LL=MIN(LF+LDM-1,LA)

    ! **  CALCULATE TYPE FLAGS                                                                                               
    IF(ISDRY >= 1 .OR. N < 5)THEN
      ! *** ICORDYU
      DO L=LF,LL
        LS=LSC(L)
        IF(SUB(L) < 0.5 .AND. SUB(LS) < 0.5) ICORDYU(L)=0
        IF(SUB(L) > 0.5 .AND. SUB(LS) > 0.5) ICORDYU(L)=1
        IF(SUB(L) < 0.5 .AND. SUB(LS) > 0.5) ICORDYU(L)=2
        IF(SUB(L) > 0.5 .AND. SUB(LS) < 0.5) ICORDYU(L)=3
      ENDDO
      ! *** ICORDXV
      DO L=LF,LL
        LW=L-1
        IF(SVB(L) < 0.5 .AND. SVB(LW) < 0.5) ICORDXV(L)=0
        IF(SVB(L) > 0.5 .AND. SVB(LW) > 0.5)THEN
          ICORDXV(L)=1
          IF(SUB(L) < 0.5) ICORDXV(L)=3
        ENDIF
        IF(SVB(L) < 0.5 .AND. SVB(LW) > 0.5) ICORDXV(L)=2
        IF(SVB(L) > 0.5 .AND. SVB(LW) < 0.5) ICORDXV(L)=3
      ENDDO
    ENDIF
  ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
  DO ND=1,NDM  
    LF=2+(ND-1)*LDM  
    LL=MIN(LF+LDM-1,LA)

    ! **  CALCULATE HORIZONTAL VELOCITY SHEARS                                                                              

    ! *** SXX+SYY DEFINED AT CELL CENTERS AND STORED IN DXU1(L,K)
    DO K=1,KC
      DO L=LF,LL
          LN=LNC(L)   
          ! *** DXU1 = dU/dX, UNITS: 1/S   
          DXU1(L,K)=SUB(L+1)*SUB(L)*(U(L+1,K)-U(L,K))/DXP(L)   !This appears to be the correct way to do things SCJ
   !      DXU1(L,K)=       SUB(L+1)*(U(L+1,K)-U(L,K))/DXP(L)   !Hamrick formulation (this may be required to do the straight west-2-east flow channel) 
          ! *** DYV1 = dV/dY, UNITS: 1/S   
          DYV1(L,K)=SVB(LN )*SVB(L)*(V(LN ,K)-V(L,K))/DYP(L)   !This appears to be the correct way to do things SCJ
   !      DYV1(L,K)=                (V(LN, K)-V(L,K))/DYP(L)   !Hamrick formulation
      ENDDO
    ENDDO

    ! *** DYU1 = dU/dY
    DO K=1,KC
      DO L=LF,LL
        LS=LSC(L)
        IF(ICORDYU(L) == 1)THEN
           DYU1(L,K)=2.*SVB(L)*SVB(LS)*(U(L,K)-U(LS,K))/(DYU(L)+DYU(LS))   !This appears to be the correct way to do things SCJ
 !         DYU1(L,K)=2.*               (U(L,K)-U(LS,K))/(DYU(L)+DYU(LS))   !Hamrick formulation
        ELSE
          DYU1(L,K)=0.
        ENDIF
        IF(ISHDMF == 2)THEN
          ! *** HMD WITH WALL EFFECTS
          IF(ICORDYU(L) == 2)THEN
            DY2DZBR=1.+0.5*DYU(LS)/ZBRWALL
            CSDRAG=0.16/((LOG(DY2DZBR))**2)
            SLIPFAC=SLIPCO*SQRT(CSDRAG)
            DYU1(L,K)=-2.*SLIPFAC*U(LS,K)/DYU(LS)
          ENDIF
          IF(ICORDYU(L) == 3)THEN
            DY2DZBR=1.+0.5*DYU(L)/ZBRWALL
            CSDRAG=0.16/((LOG(DY2DZBR))**2)
            SLIPFAC=SLIPCO*SQRT(CSDRAG)
            DYU1(L,K)=2.*SLIPFAC*U(L,K)/DYU(L)
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    ! *** DXV1 = dV/dX
    DO K=1,KC
      DO L=LF,LL
        LW=L-1
        IF(ICORDXV(L) == 1)THEN
          DXV1(L,K)=2.*SUB(L)*SUB(LW)*(V(L,K)-V(LW,K))/(DXV(L)+DXV(LW))   !This appears to be the correct way to do things SCJ
    !     DXV1(L,K)=2.*               (V(L,K)-V(LW,K))/(DXV(L)+DXV(LW))   !Hamrick formulation
        ELSE
          DXV1(L,K)=0.
        ENDIF
        IF(ISHDMF == 2)THEN
          ! *** WALL EFFECTS
          IF(ICORDXV(L) == 2)THEN
            DX2DZBR=1.+0.5*DXV(LW)/ZBRWALL
            CSDRAG=0.16/((LOG(DX2DZBR))**2)
            SLIPFAC=SLIPCO*SQRT(CSDRAG)
            DXV1(L,K)=-2.*SLIPFAC*V(LW,K)/DXV(LW)
          ENDIF
          IF(ICORDXV(L) == 3)THEN
           DX2DZBR=1.+0.5*DXV(L)/ZBRWALL
           CSDRAG=0.16/((LOG(DX2DZBR))**2)
           SLIPFAC=SLIPCO*SQRT(CSDRAG)
           DXV1(L,K)=2.*SLIPFAC*V(L,K)/DXV(L)
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    IF(NQWR < 1)THEN
      ! *** SXY = dU/dY + dV/dX
      DO K=1,KC
        DO L=LF,LL
          SXY(L,K)=DYU1(L,K)+DXV1(L,K)
          DYU1L=DYU1(L,K)
          DXV1L=DXV1(L,K)
          SXYL=SXY(L,K)
        ENDDO
      ENDDO
    ENDIF
    
  ENDDO  ! *** END OF DOMAIN
!$OMP END PARALLEL DO  
  ! *** WITHDRAWAL/RETURN
  IF(NQWR > 0)THEN
    DO IWR=1,NQWR
      ! *** Handle +/- Flows for Withdrawal/Return Structures
      NQSTMP=NQWRSERQ(IWR)  
      IF( QWRSERT(NQSTMP) >= 0. )THEN
        ! *** Original Withdrawal/Return
        IU=IQWRU(IWR)  
        JU=JQWRU(IWR)  
        KU=KQWRU(IWR)  
      ELSE
        ! *** Reverse Flow Withdrawal/Return
        IU=IQWRD(IWR)  
        JU=JQWRD(IWR)  
        KU=KQWRD(IWR) 
      ENDIF
      DXU1(LIJ(IU,JU),KU)=0.0
      DXV1(LIJ(IU,JU),KU)=0.0
      DYU1(LIJ(IU,JU),KU)=0.0
      DYV1(LIJ(IU,JU),KU)=0.0
    ENDDO
    
    ! *** SXY = dU/dY + dV/dX
    DO ND=1,NDM  
      LF=2+(ND-1)*LDM  
      LL=MIN(LF+LDM-1,LA)
      DO K=1,KC
        DO L=2,LA
          SXY(L,K)=DYU1(L,K)+DXV1(L,K)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
!$OMP PARALLEL DO  
  DO ND=1,NDM  
    LF=2+(ND-1)*LDM  
    LL=MIN(LF+LDM-1,LA)

    IF(AHD > 0.0)THEN
      ! *** CALCULATE SMAGORINSKY HORIZONTAL VISCOSITY
      DO K=1,KC
        DO L=LF,LL
         TMPVAL=AHD*DXP(L)*DYP(L)
         DSQR=DXU1(L,K)*DXU1(L,K)+DYV1(L,K)*DYV1(L,K)+SXY(L,K)*SXY(L,K)*0.5
         AH(L,K)=AHO+TMPVAL*SQRT(DSQR)
       ENDDO
      ENDDO
    ELSEIF(N < 10 .OR. ISWAVE == 2 .OR. ISWAVE==4)THEN
      ! *** ONLY NEED TO ASSIGN INITIALLY
      DO K=1,KC
        DO L=LF,LL
         AH(L,K)=AHO
       ENDDO
      ENDDO
    ENDIF

    ! **  CALCULATE HORIZONTAL SMAG DIFFUSION DUE TO WAVE BREAKING                                                               
    IF(ISWAVE == 2 .OR. ISWAVE==4)THEN
      IF(WVLSH > 0.0 .OR. WVLSX > 0.0)THEN
        IF(ISWAVE == 2 .AND. N < NTSWV)THEN
          TMPVAL=FLOAT(N)/FLOAT(NTSWV)
          WVFACT=0.5-0.5*COS(PI*TMPVAL)
        ELSE
          WVFACT=1.0
        ENDIF

        IF(ISDRY > 0)THEN
          DO K=1,KC
            DO L=LF,LL
              IF(LWVMASK(L))THEN
                IF(LMASKDRY(L))THEN
                  IF (WVDISP(L,K)>0) THEN
                    DTMPH=WVDISP(L,K)**0.3333
                  ELSE
                    DTMPH=0
                  ENDIF
                  TMPVAL=2.*PI/WVFRQL(L)     ! *** WAVE PERIOD
                  AHWVX=WVLSX*TMPVAL*TMPVAL  
                  DTMPX=WVDISP(L,K)/HP(L) 
                  AH(L,K)=AH(L,K)+WVFACT*(WVLSH*DTMPH*HP(L)+AHWVX*DTMPX) 
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ELSE
          DO K=1,KC
            DO L=LF,LL
              IF(LWVMASK(L))THEN
                IF (WVDISP(L,K)>0) THEN
                  DTMPH=WVDISP(L,K)**0.3333
                ELSE
                  DTMPH=0
                ENDIF
                TMPVAL=2.*PI/WVFRQL(L)
                AHWVX=WVLSX*TMPVAL*TMPVAL  
                DTMPX=WVDISP(L,K)/HP(L)   
                AH(L,K)=AH(L,K)+WVFACT*(WVLSH*DTMPH*HP(L)+AHWVX*DTMPX)  
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDIF
    ENDIF
  ENDDO  ! *** END OF DOMAIN
!$OMP END PARALLEL DO
!$OMP PARALLEL DO  
  DO ND=1,NDM  
    LF=2+(ND-1)*LDM  
    LL=MIN(LF+LDM-1,LA)

    ! **  CALCULATE DIFFUSIVE MOMENTUM FLUXES                                                                               
    DO K=1,KC
      DO L=LF,LL
        LN=LNC(L)
        LS=LSC(L)
        FMDUX(L,K)=2.0*SUB(L  )*(DYP(L  )*HP(L  )*AH(L  ,K)*DXU1(L ,K) - DYP(L-1)*HP(L-1)*AH(L-1,K)*DXU1(L-1,K))
        FMDUY(L,K)=    SVB(LN )*(DXU(LN )*HU(LN )*AH(LN ,K)*SXY(LN ,K) - DXU(L  )*HU(L  )*AH(L  ,K)*SXY(L,   K))
        FMDVY(L,K)=2.0*SVB(L  )*(DXP(L  )*HP(L  )*AH(L  ,K)*DYV1(L ,K) - DXP(LS )*HP(LS )*AH(LS ,K)*DYV1(LS ,K))
        FMDVX(L,K)=    SUB(L+1)*(DYV(L+1)*HV(L+1)*AH(L+1,K)*SXY(L+1,K) - DYV(L  )*HV(L  )*AH(L  ,K)*SXY(L   ,K))
      ENDDO
    ENDDO

    ! *** TREAT THE NORTH & WEST WALL SLIPPAGE
    IF(ISHDMF == 2)THEN
      DO K=1,KC
        DO L=LF,LL
          LN=LNC(L)
          IF(SVBO(LN) < 0.5)THEN
            DY2DZBR=1.+0.5*DYU(L)/ZBRWALL
            CSDRAG=0.16/((LOG(DY2DZBR))**2)
            SLIPFAC=SLIPCO*SQRT(CSDRAG)
            SXYLN=-2.*SLIPFAC*U(L,K)/DYU(L)
            FMDUY(L,K)=DXU(L)*HP(L)*AH(L,K)*(SXYLN-SXY(L ,K))
          ENDIF
          IF(SUBO(L+1) < 0.5)THEN
            DX2DZBR=1.+0.5*DXV(L)/ZBRWALL
            CSDRAG=0.16/((LOG(DX2DZBR))**2)
            SLIPFAC=SLIPCO*SQRT(CSDRAG)
            SXYEE=-2.*SLIPFAC*V(L,K)/DXV(L)
            FMDVX(L,K)=DYV(L)*HP(L)*AH(L,K)*(SXYEE-SXY(L,K))
          ENDIF
        ENDDO
      ENDDO
    ENDIF
  ENDDO  ! *** END OF DOMAIN
!$OMP END PARALLEL DO
  
  ! *** ZERO BOUNDARY CELL MOMENTUM DIFFUSION
  DO LL=1,NBCS
    L=LBCS(LL)
    DO K=1,KC
      FMDUX(L,K)=0.0
      FMDUY(L,K)=0.0
      FMDVY(L,K)=0.0
      FMDVX(L,K)=0.0
    ENDDO
  ENDDO

  RETURN

END
