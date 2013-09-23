SUBROUTINE MHKPWRDIS  
!**********************************************************************C  
! **  SUBROUTINE MHKPWRDIS CALCULATES POWER DISSIPATION FROM MARINE HYDROKINETIC  
! **  DEVICES AS A FUNCTION OF WATER VELOCITY  
!
!     04/2010  Bill Arnold and Scott James
!              Sandia National Laboratories
!  
!**********************************************************************C  
  USE GLOBAL
  IMPLICIT NONE
  REAL::THRSTCOEF=0.0,ZTOP=0.0,ZBOTTOM=0.0
  REAL::FS_WF,FS_NF,FS_EF,FS_SF,MAXSPD,SPDN,SPDS,SPDE,SPDW,AVGSPD
  REAL::VELUP,UVEC,VVEC !,UVELUP,VVELUP
  REAL::FMHK,FSUP,UATVFACE,VATUFACE,UATVFACEN,VATUFACEE
  REAL::LAYFRACM(KC),LAYFRACS(KC),NEGLAYFRACM(KC),NEGLAYFRACS(KC),FLOWSPEED(KC)
  REAL::AWEIGHTXW,AWEIGHTXE,AWEIGHTYS,AWEIGHTYN,SUMLAYM,SUMLAYS,SUMNEGLAYM,SUMNEGLAYS
  INTEGER::MHKCOUNT,M,L,K,LW,LE,LS,LN,LNW,LSE,LNE
  LOGICAL::STATUS
  LOGICAL,SAVE::FIRSTTIME=.FALSE.
!**********************************************************************C 
  INQUIRE(FILE='POWERBUG.DAT',EXIST=STATUS)
  IF(STATUS.AND.DEBUG.AND..NOT.FIRSTTIME)THEN
    OPEN(357,FILE='POWERBUG.DAT') !this file is for debugging purposes
    FIRSTTIME=.TRUE.
  ENDIF
  MHKCOUNT=0 !initialize the running count of the number of cells with MHKdevices (will equal TCOUNT)
  FXMHKE(:)=0.0;FYMHKE(:)=0.0;FXSUPE(:)=0.0;FYSUPE(:)=0.0 !initialize x,y-forces from MHK device/support for the external solution
  FXMHK(:,:)=0.0;FYMHK(:,:)=0.0;FXSUP(:,:)=0.0;FYSUP(:,:)=0.0 !initialize x,y-forces from the MHK device/support for the internal solution
  PMHK(:,:)=0.0;PSUP(:,:)=0.0 !initialize power extracted from the MHK device/support from flow
  FLOWSPEED(:)=0.0 !Initialize flow speed
  DO L=2,LA  !loop through model area
    IF(.NOT.LMASKDRY(L).OR.MVEGL(L)<90)CYCLE !if the cell is dry or does not contain an MHK, skip this cell
    MHKCOUNT=MHKCOUNT+1
    LW=L-1      !west cell, I-1,J
    LE=L+1      !east cell, I+1,J
    LS=LSC(L)   !south cell, I,J-1
    LN=LNC(L)   !north cell, I,J+1
    LNW=LNWC(L) !northwest cell, I-1,J+1
    LNE=LNEC(L) !northeast cell, I+1,J+1
    LSE=LSEC(L) !southeast cell, I+1,J-1
    M=MVEGL(L)-90 !This was put in to have MHKs be vegetative inputs > 90; this is the mhktype
    LAYFRACM(:)=0.0;NEGLAYFRACM(:)=0.0 !initialize the layer fraction variable MHK
    LAYFRACS(:)=0.0;NEGLAYFRACS(:)=0.0 !initialize the layer fraction variable support
    DO K=1,KC !MHK device layer filler - which layers does the device occupy and at what fraction
      ZTOP=HP(L)*Z(K)+BELV(L) !layer top elevation
      IF(ZTOP<ZMINMHK(M,L))CYCLE !layer is below device
      ZBOTTOM=HP(L)*Z(K-1)+BELV(L) !layer bottom elevation
      IF(ZBOTTOM>ZMAXMHK(M,L))CYCLE !layer is above device
      IF(ZTOP>=ZMAXMHK(M,L).AND.ZBOTTOM<=ZMINMHK(M,L))THEN !device is wholly contained in this layer (special case)
        LAYFRACM(K)=(ZMAXMHK(M,L)-ZMINMHK(M,L))/(HP(L)*DZC(K)) !calculate fraction of layer that is occupied
        EXIT
      ENDIF
      IF(ZMAXMHK(M,L)>=ZTOP.AND.ZMINMHK(M,L)<=ZBOTTOM)THEN !this layer is fully occupied by the device
        LAYFRACM(K)=1.0
        CYCLE
      ENDIF
      IF(ZBOTTOM<ZMINMHK(M,L).AND.ZMAXMHK(M,L)>=ZTOP)THEN !this layer is partially occupied by the device (bottom)
        LAYFRACM(K)=(ZTOP-ZMINMHK(M,L))/(HP(L)*DZC(K)) !calculate the fraction of layer that is occupied
        CYCLE
      ENDIF
      IF(ZTOP>=ZMAXMHK(M,L).AND.ZMINMHK(M,L)<ZBOTTOM)THEN !this layer is partially occupied by the device (top)
        LAYFRACM(K)=(ZMAXMHK(M,L)-ZBOTTOM)/(HP(L)*DZC(K)) !calculate the fraction of layer that is occupied
        CYCLE
      ENDIF
    ENDDO
    NEGLAYFRACM(:)=LAYFRACM(:)-1.0 !negative of the layer fraction occupied by the MHK turbine
    SUMLAYM=SUM(LAYFRACM(1:KC));SUMNEGLAYM=SUM(NEGLAYFRACM(1:KC)) !Sum of MHK layer fractions
    DO K=1,KC !MHK support layer filler - which layers does the support occupy and at what fraction
      ZTOP=HP(L)*Z(K)+BELV(L) !layer top elevation
      IF(ZTOP<ZMINSUP(M,L))CYCLE !layer is below support
      ZBOTTOM=HP(L)*Z(K-1)+BELV(L) !layer bottom elevation
      IF(ZBOTTOM>ZMAXSUP(M,L))CYCLE !layer is above support
      IF(ZTOP>=ZMAXSUP(M,L).AND.ZBOTTOM<=ZMINSUP(M,L))THEN !support is wholly contained in this layer (special case)
        LAYFRACS(K)=(ZMAXSUP(M,L)-ZMINSUP(M,L))/(HP(L)*DZC(K)) !calculate fraction of layer that is occupied
        EXIT
      ENDIF
      IF(ZMAXSUP(M,L)>=ZTOP.AND.ZMINSUP(M,L)<=ZBOTTOM)THEN !this layer is fully occupied by the support
        LAYFRACS(K)=1.0
        CYCLE
      ENDIF
      IF(ZBOTTOM<ZMINSUP(M,L).AND.ZMAXSUP(M,L)>=ZTOP)THEN !this layer is partially occupied by the support (bottom)
        LAYFRACS(K)=(ZTOP-ZMINSUP(M,L))/(HP(L)*DZC(K)) !calculate the fraction of layer that is occupied
        CYCLE
      ENDIF
      IF(ZTOP>=ZMAXSUP(M,L).AND.ZMINSUP(M,L)<ZBOTTOM)THEN !this layer is partially occupied by the support (top)
        LAYFRACS(K)=(ZMAXSUP(M,L)-ZBOTTOM)/(HP(L)*DZC(K)) !calculate the fraction of layer that is occupied
        CYCLE
      ENDIF
    ENDDO
    NEGLAYFRACS(:)=LAYFRACS(:)-1.0  
    SUMLAYS=SUM(LAYFRACS(1:KC));SUMNEGLAYS=SUM(NEGLAYFRACS(1:KC)) !Sum of support layer fractions
    DO K=1,KC
      IF(ZMAXMHK(M,L)>HP(L)+BELV(L))THEN !Protruding from water?
        PRINT*,'MHK DEVICE PROTRUDING FROM WATER',HP(L),ZMAXMHK(M,L),L,IL(L),JL(L)
        STOP
      ELSEIF(ZMINMHK(M,L)<BELV(L))THEN !Below bed elevation?
        PRINT*,'MKH DEVICE IN SEDIMENT',BELV(L),ZMINMHK(M,L),L,IL(L),JL(L)
        STOP
      ENDIF
      UVEC=0.5*(U(L,K)+U(LE,K))      !I,J cell center u-speed
      VVEC=0.5*(V(L,K)+V(LN,K))      !I,J cell center v-speed
      FLOWSPEED(K)=SQRT(UVEC*UVEC+VVEC*VVEC) !I,J cell center speed
      IF((LAYFRACM(K)==0.0.AND.LAYFRACS(K)==0.0).OR.FLOWSPEED(K)<1.0E-03)CYCLE !no MHK or support or velocity in this layer
      FMHK=0.0;FSUP=0.0 !initialize variables
      IF(UPSTREAM==1)THEN !use the upstream flowspeed to assess power extraction
        UATVFACE= 0.25*(U(L,K)+U(LE,K)+U(LS,K)+U(LSE,K))  !u-velocity at south face (the v-face)
        VATUFACE= 0.25*(V(L,K)+V(LW,K)+V(LN,K)+V(LNW,K))  !v-velocity at west face  (the u-face)
        UATVFACEN=0.25*(U(L,K)+U(LE,K)+U(LN,K)+U(LNE,K))  !u-velocity at north face (the u-north-face)
        VATUFACEE=0.25*(V(L,K)+V(LE,K)+V(LN,K)+V(LNE,K))  !v-velocity at east face  (the v-east-face)
        FS_WF=U(L ,K);FS_EF=U(LE,K) !velocities on the west/east   faces (u velocities into the cell)
        FS_SF=V(L ,K);FS_NF=V(LN,K) !velocities on the south/north faces (v velocities into the cell) !Bug found! Had FS_SF=V(LS,K)
        SPDN=SQRT(UATVFACEN*UATVFACEN+V(LN,K)*V(LN,K)) !speed at north face
        SPDS=SQRT(UATVFACE*UATVFACE+V(L,K)*V(L,K)) !speed at south face
        SPDE=SQRT(U(LE,K)*U(LE,K)+VATUFACEE*VATUFACEE) !speed at east face
        SPDW=SQRT(U(L,K)*U(L,K)+VATUFACE*VATUFACE) !speed at west face
        IF(FS_NF>-0.01)SPDN=0.0 !flow is OUT of north face
        IF(FS_SF< 0.01)SPDS=0.0 !flow is OUT of south face
        IF(FS_WF< 0.01)SPDW=0.0 !flow is OUT of west face
        IF(FS_EF>-0.01)SPDE=0.0 !flow is OUT of east face
        MAXSPD=MAX(SPDN,SPDS,SPDE,SPDW) !identify maximum speed
!        VVELUP = MAX(SPDN,SPDS)
!        IF(V(L,K)<0)VVELUP=-VVELUP 
!        UVELUP = MAX(SPDW,SPDE)
!        IF(U(L,K)<0)UVELUP=-UVELUP
        IF(MAXSPD==SPDN)THEN !what face is it on?
          VELUP=SQRT((0.25*(U(LN,K)+U(LNE,K)+U(LNC(LN),K)+U(LNC(LN)+1,K)))**2+V(LNC(LN),K)**2) 
        ELSEIF(MAXSPD==SPDS)THEN !South
          VELUP=SQRT((0.25*(U(LS,K)+U(LSE,K)+U(LSC(LS),K)+U(LSC(LS)+1,K)))**2+V(LS,K)**2)
        ELSEIF(MAXSPD==SPDE)THEN !East
          VELUP=SQRT(U(LE+1,K)**2+(0.25*(V(LE,K)+V(LE+1,K)+V(LNE,K)+V(MIN(LC,LN+2),K)))**2)
        ELSE !West
          VELUP=SQRT(U(LW  ,K)**2+(0.25*(V(LW,K)+V(LW-1,K)+V(LNW,K)+V(LN-2,K)))**2)
        ENDIF
      ELSEIF(UPSTREAM==0)THEN !use the local cell's flowspeed to assess power extraction
!        UVELUP=UVEC
!        VVELUP=VVEC
        VELUP=FLOWSPEED(K)
      ENDIF
      IF(LAYFRACM(K)>0.0)THEN !MHK device exists in this layer 
        IF(VELUP<VMINCUT(M))THEN !no power generation
          THRSTCOEF=0.0 !no need for these calcs
        ELSEIF(VELUP<VMAXCUT(M))THEN !optimal power generation
          THRSTCOEF=CTMHK(M)
        ELSE !superoptimal flow speed limits power generation to VMAXCUT
          THRSTCOEF=CTMHK(M)
          VELUP=VMAXCUT(M)
        ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        VELUP=U(L-10,K)  !Special case for calibration for flow from west to east
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!FMHK=0.5*ThrustCoef*Area*(U_inf)^2 where U_inf is the upstream velocity, VELUP [m^4/s^2]
        FMHK=0.5*LAYFRACM(K)*THRSTCOEF*VELUP*VELUP*HP(L)*DZC(K)*WIDTHMHK(M) !area is ASSUMED square
!PMHK=FMHK*U where U is the local flowspeed
        PMHK(L,K)=FMHK*FLOWSPEED(K) !ThrustCoef*|u|u^2*area [m^5/s^3] (will yield different power outputs depending on UPSTREAM)
        AWEIGHTXW=DYU(L)*HU(L)/(DYU(L)*HU(L)+DYU(LE)*HU(LE));AWEIGHTXE=1.0-AWEIGHTXW !area-weight for west/east faces
        AWEIGHTYS=DXV(L)*HV(L)/(DXV(L)*HV(L)+DXV(LN)*HV(LN));AWEIGHTYN=1.0-AWEIGHTYS !area-weight for south/north faces
!To get the x and y components, multiply by a velocity vector divided by the local flow speed FXMHK=FMHK*UVEC/FLOWSPEED(K)
        FXMHK(L ,K)=FXMHK(L ,K)+AWEIGHTXW*SUB(L )*FMHK*UVEC/FLOWSPEED(K) !SUB(L)*FMHK(L,K)*(Uvel/q) [m^4/s^2]
        FXMHK(LE,K)=FXMHK(LE,K)+AWEIGHTXE*SUB(LE)*FMHK*UVEC/FLOWSPEED(K) !distribute forces on each U-face of the cell
        FYMHK(L ,K)=FYMHK(L ,K)+AWEIGHTYS*SVB(L )*FMHK*VVEC/FLOWSPEED(K) !y components of "forces" [m^4/s^2]
        FYMHK(LN,K)=FYMHK(LN,K)+AWEIGHTYN*SVB(LN)*FMHK*VVEC/FLOWSPEED(K) !distribute forces on each V-face of the cell
        IF(BSC>0.0)THEN !if variable density, take it into account
          PMHK(L,K)=PMHK(L,K)*(B(L,K)+1.0)*1000.0
        ELSE
          PMHK(L,K)=PMHK(L,K)*1024. !density of seawater is ~1024kg/m^3
        ENDIF
        IF(DEBUG.AND..NOT.UPSTREAM.AND.MOD(N,100)==0)WRITE(357,'(I5,1X,3(I3,1X),5(E12.5,1X))')N,IL(L),JL(L),K,PMHK(L,K),VELUP,UVEC,VVEC,HP(L)
!        IF(DEBUG.AND.     UPSTREAM.AND.MOD(N,100)==0)WRITE(357,'(I5,1X,3(I3,1X),5(E12.5,1X))')N,IL(L),JL(L),K,PMHK(L,K),VELUP,UVEC,VVEC,HP(L)
        IF(DEBUG.AND.M==1.AND.K==3.AND.MOD(N,1)==0)WRITE(357,'(I5,1X,3(I3,1X),5(E12.5,1X))')N,IL(L),JL(L),K,PMHK(L,K),VELUP,FXMHK(L ,K),FYMHK(L ,K),HP(L)
      ENDIF
      IF(LAYFRACS(K)>0.0)THEN !MHK support exists in this layer
        FSUP=0.5*LAYFRACS(K)*CDSUP(M)*FLOWSPEED(K)*FLOWSPEED(K)*HP(L)*DZC(K)*WIDTHSUP(M) !calculate the force on a cell [m^4/s^2]
        PSUP(L,K)=FSUP*FLOWSPEED(K) !0.5*C_d*Asup*u^3 [m^5/s^3]
        AWEIGHTXW=DYU(L)*HU(L)/(DYU(L)*HU(L)+DYU(LE)*HU(LE));AWEIGHTXE=1.0-AWEIGHTXW !area-weight for west/east face
        AWEIGHTYS=DXV(L)*HV(L)/(DXV(L)*HV(L)+DXV(LN)*HV(LN));AWEIGHTYN=1.0-AWEIGHTYS !area-weight for south/north face
        FXSUP(L ,K)=FXSUP(L ,K)+AWEIGHTXW*SUB(L )*FSUP*UVEC/FLOWSPEED(K) !x-component [m^4/s^2] Note that FLOWSPEED(K) is multiplied in but divided back out again when normalizing by UVEC/FLOWSPEED(K)
        FXSUP(LE,K)=FXSUP(LE,K)+AWEIGHTXE*SUB(LE)*FSUP*UVEC/FLOWSPEED(K) !distribute forces on both U-faces of the cell 
        FYSUP(L ,K)=FYSUP(L ,K)+AWEIGHTYS*SVB(L )*FSUP*VVEC/FLOWSPEED(K) !y-component [m^4/s^2] Note that FLOWSPEED(K) is multiplied in but divided back out again when normalizing by VVEC/FLOWSPEED(K)
        FYSUP(LN,K)=FYSUP(LN,K)+AWEIGHTYN*SVB(LN)*FSUP*VVEC/FLOWSPEED(K) !distribute forces on both V-faces of the cell
        IF(BSC>0.0)THEN !if variable density, take it into account
          PSUP(L,K)=PSUP(L,K)*(B(L,K)+1.0)*1000.0
        ELSE
          PSUP(L,K)=PSUP(l,K)*1024.
       ENDIF
     ENDIF
    ENDDO
    DO K=1,KC
      IF(SUMLAYS==0.0)THEN !No total force can be added to the internal-mode solution the way this is written, the sum across layers is zero. Internal forces are directional, so the sum of FXMHK,FYMHK,FXSUP,FYSUP are used
        FX(L ,K)=FX(L ,K)+PB_COEF*SUM(FXMHK(L ,1:KC))*(LAYFRACM(K)/SUMLAYM-NEGLAYFRACM(K)/SUMNEGLAYM) !pull x-force out of MHK layer for internal mode (no support structure) - push forces in other layers
        FX(LE,K)=FX(LE,K)+PB_COEF*SUM(FXMHK(LE,1:KC))*(LAYFRACM(K)/SUMLAYM-NEGLAYFRACM(K)/SUMNEGLAYM) !pull x-force out of MHK layer for internal mode (east face) - push forces in other layers
        FY(L ,K)=FY(L ,K)+PB_COEF*SUM(FYMHK(L ,1:KC))*(LAYFRACM(K)/SUMLAYM-NEGLAYFRACM(K)/SUMNEGLAYM) !pull y-force out of MHK layer for internal mode (no support structure) - push forces in other layers
        FY(LN,K)=FY(LN,K)+PB_COEF*SUM(FYMHK(LN,1:KC))*(LAYFRACM(K)/SUMLAYM-NEGLAYFRACM(K)/SUMNEGLAYM) !pull y-force out of MHK layer for internal mode (north face) - push forces in other layers
      ELSE
        FX(L ,K)=FX(L ,K)+(PB_COEF*SUM(FXMHK(L ,1:KC))*(LAYFRACM(K)/SUMLAYM-NEGLAYFRACM(K)/SUMNEGLAYM)+SUM(FXSUP(L ,1:KC))*(LAYFRACS(K)/SUMLAYS-NEGLAYFRACS(K)/SUMNEGLAYS)) !pull x-force out of MHK/support layer for internal mode - push forces in other layers
        FX(LE,K)=FX(LE,K)+(PB_COEF*SUM(FXMHK(LE,1:KC))*(LAYFRACM(K)/SUMLAYM-NEGLAYFRACM(K)/SUMNEGLAYM)+SUM(FXSUP(LE,1:KC))*(LAYFRACS(K)/SUMLAYS-NEGLAYFRACS(K)/SUMNEGLAYS)) !pull x-force out of MHK/support layer for internal mode - push forces in other layers
        FY(L ,K)=FY(L ,K)+(PB_COEF*SUM(FYMHK(L ,1:KC))*(LAYFRACM(K)/SUMLAYM-NEGLAYFRACM(K)/SUMNEGLAYM)+SUM(FYSUP(L ,1:KC))*(LAYFRACS(K)/SUMLAYS-NEGLAYFRACS(K)/SUMNEGLAYS)) !pull x-force out of MHK/support layer for internal mode - push forces in other layers
        FY(LN,K)=FY(LN,K)+(PB_COEF*SUM(FYMHK(LN,1:KC))*(LAYFRACM(K)/SUMLAYM-NEGLAYFRACM(K)/SUMNEGLAYM)+SUM(FYSUP(LN,1:KC))*(LAYFRACS(K)/SUMLAYS-NEGLAYFRACS(K)/SUMNEGLAYS)) !pull x-force out of MHK/support layer for internal mode - push forces in other layers
      ENDIF
    ENDDO
    FXMHKE(L)=FXMHKE(L)+SUM(ABS(FXMHK(L,1:KC)));FXMHKE(LE)=FXMHKE(LE)+SUM(ABS(FXMHK(LE,1:KC))) !Sum layer force magnitudes for external mode solution (need absolute value of forces)
    FYMHKE(L)=FYMHKE(L)+SUM(ABS(FYMHK(L,1:KC)));FYMHKE(LN)=FYMHKE(LN)+SUM(ABS(FYMHK(LN,1:KC))) !Sum layer force magnitudes for external mode solution (these forces are later multiplied by a directional velocity so they need to be absolute values here)
    FXSUPE(L)=FXSUPE(L)+SUM(ABS(FXSUP(L,1:KC)));FXSUPE(LE)=FXSUPE(LE)+SUM(ABS(FXSUP(LE,1:KC))) !Sum layer force magnitudes for external mode solution (absolute values because these are later multiplied by the local velocity to apply a direction)
    FYSUPE(L)=FYSUPE(L)+SUM(ABS(FYSUP(L,1:KC)));FYSUPE(LN)=FYSUPE(LN)+SUM(ABS(FYSUP(LN,1:KC))) !Sum layer force magnitudes for external mode solution
    AVGSPD=SUM(FLOWSPEED(1:KC)*DZC(1:KC))
    IF(AVGSPD==0.0)CYCLE
!CALEXP2T is expecting units of [m^3/s] for FXMHKE, which is the sum of absolute values of FXMHK
!CALEXP2T divides by water-column volume before passing this "force" onto FUHDYE (in units of [1/s]), which is used for momentum conservation in CALPUV
!Units of FXMHKE (etc) are same as FX and FXMHK (etc) [m^4/s^2] so they must be divided by the average speed in this water column
!Because these forces should be close to zero when the U or V speeds are zero, we need to include a normalized form as either U/AVGSPD or V/AVGSPD.  Without this “directional normalization,” when AVGSPD is small, external mode forces are too large.
!Multiply by a directional velocity normalized by AVGSPD and then divide by AVGSPD to get the units correct
    FXMHKE(L )=FXMHKE(L )*ABS(SUM(U(L ,1:KC)*DZC(1:KC)))/AVGSPD/AVGSPD !external mode solution units of [m^3/s]
    FXSUPE(L )=FXSUPE(L )*ABS(SUM(U(L ,1:KC)*DZC(1:KC)))/AVGSPD/AVGSPD !external mode solution units of [m^3/s]
    FXMHKE(LE)=FXMHKE(LE)*ABS(SUM(U(LE,1:KC)*DZC(1:KC)))/AVGSPD/AVGSPD !external mode solution units of [m^3/s]
    FXSUPE(LE)=FXSUPE(LE)*ABS(SUM(U(LE,1:KC)*DZC(1:KC)))/AVGSPD/AVGSPD !external mode solution units of [m^3/s]
    FYMHKE(L )=FYMHKE(L )*ABS(SUM(V(L ,1:KC)*DZC(1:KC)))/AVGSPD/AVGSPD !external mode solution units of [m^3/s] 
    FYSUPE(L )=FYSUPE(L )*ABS(SUM(V(L ,1:KC)*DZC(1:KC)))/AVGSPD/AVGSPD !external mode solution units of [m^3/s]
    FYMHKE(LN)=FYMHKE(LN)*ABS(SUM(V(LN,1:KC)*DZC(1:KC)))/AVGSPD/AVGSPD !external mode solution units of [m^3/s]
    FYSUPE(LN)=FYSUPE(LN)*ABS(SUM(V(LN,1:KC)*DZC(1:KC)))/AVGSPD/AVGSPD !external mode solution units of [m^3/s]
    EMHK(MHKCOUNT,L)=EMHK(MHKCOUNT,L)+DT*SUM(PMHK(L,1:KC))*2.7778E-10 !factor converts to MW-hr
    ESUP(MHKCOUNT,L)=ESUP(MHKCOUNT,L)+DT*SUM(PSUP(L,1:KC))*2.7778E-10 !factor converts to MW-hr
  ENDDO
  RETURN
END SUBROUTINE MHKPWRDIS



