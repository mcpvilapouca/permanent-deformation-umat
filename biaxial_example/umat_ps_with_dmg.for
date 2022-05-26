C********************************************************************
C Record of revisions:                                              |
C        Date        Programmer        Description of change        |
C        ====        ==========        =====================        |
C     09/12/2019  Maria Vila Pouca       HGO + GROWTH MODEL         |
C                                       GROWTH IN ONLY ONE          |
C                                         FIBER DIRECTION:4         |
C     16/01/2020  Maria Vila Pouca     corrected so that when       |
C                                     STG4=STGMAX4 STG4=STG40       |
C     17/02/2021  Maria Vila Pouca    Permanent set test, evolution |
C                                       equal to growth            |
C     25/02/2021  Maria Vila Pouca     Non-local fatigue Dmg added  |
C     23/03/2021 Maria Vila Pouca    Change Permament set evolution |
C                                    dependent of Damage            |
C     06/04/2021   Maria Vila Pouca   Corrections on DDSDDE         |
C     07/04/2021  Maria Vila Pouca    Final Corrections             |
C     12/04/2021  Maria Vila Pouca    Final Corrections             |
C     19/04/2021  Maria Vila Pouca    PHI_D added                   |
C     20/04/2021  Maria  Vila Pouca   beta0 added to PS evol        |
C     05/05/2021  Maria Vila Pouca    Corrections on DDSDDE         |
C     29/06/2021  Maria Vila Pouca    correction: zd update ps evol |
C     08/07/2021  Maria Vila Pouca    remove (1-D) contribution     |
C     12/07/2021  Maria Vila Pouca     D=1 maximum ps evolution     |
C     13/08/2021  Maria Vila Pouca    small corrections when D=1    |
C--------------------------------------------------------------------
C     Description:
C     UEXTERNALDB: mount database; fibers directions;.inc, .inp files
C     UMAT: NH matrix + Holzapfel fibers + permanent deformation
C--------------------------------------------------------------------
      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
      INCLUDE 'PARAM_UMAT.F'
C       this subroutine get the fibers directions resorting the
C      .inc files and read the fibers direction
C
C     UEXTERNAL just called once; work in parallel computing
      COMMON /KFIB/FIBORI
      COMMON /KFIB/FIBORI2
      DIMENSION TIME(2)
      REAL*8 FIBORI(NELEM,4),FIBORI2(NELEM,4)
      CHARACTER(256) FILENAME
      CHARACTER(256) JOBDIR
C     LOP=0 --> START OF THE ANALYSIS
      IF(LOP.EQ.0.OR.LOP.EQ.4) THEN
C
        CALL GETOUTDIR(JOBDIR,LENJOBDIR)
C
          FILENAME=JOBDIR(:LENJOBDIR)//'/'//DIR1
          OPEN(15,FILE=FILENAME)
          DO I=1,NELEM
             READ(15,*) (FIBORI(I,J),J=1,4)
          END DO
           CLOSE(15)
C
          FILENAME=JOBDIR(:LENJOBDIR)//'/'//DIR2
          OPEN(15,FILE=FILENAME)
          DO I=1,NELEM
             READ(15,*) (FIBORI2(I,J),J=1,4)
          END DO
           CLOSE(15)
      END IF
C
      RETURN
C
      END
C---------------------------------------------------------------------
C---------------------------------------------------------------------
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.F'
C
C
      CHARACTER*8 CMNAME
      DOUBLE PRECISION STRESS(NTENS),STATEV(NSTATEV),
     1 DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),
     2 TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3),
     3 DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),DDSDDE(NTENS,NTENS),
     4 FIBORI(NELEM,4),FIBORI2(NELEM,4),CCMAT1(NTENS,NTENS)
C
      DOUBLE PRECISION SSE, SPD, SCD, RPL, DRPLDT, DTIME, TEMP,
     1   DTEMP,PNEWDT,CELENT
C
      INTEGER NDI, NSHR, NTENS, NSTATEV, NPROPS, NOEL, NPT,
     1 LAYER, KSPT, KSTEP, KINC
C
C
      COMMON /KFIB/FIBORI
      COMMON /KFIB/FIBORI2
C---------------------------------------------------------------------
C     LOCAL ARRAYS
C---------------------------------------------------------------------
C     STATEV - STATE VARIABLES ARRAY
C     UNIT   - IDENTITY TENSOR
C     II1    - POINTERS VECTOR
C     DFGRD1E - ELASTIC DEFORMATION GRADIENT
C     FG4 - GROWTH TENSOR
C     DISTGR - DEVIATORIC DEFORMATION GRADIENT (DISTORTION TENSOR)
C     C      -  RIGHT CAUCHY-GREEN TENSOR
C     CE      -  ELASTIC RIGHT CAUCHY-GREEN TENSOR
C     CBAR   - DEVIATORIC RIGHT CAUCHY-GREEN TENSOR
C     CBARE   - ELASTIC DEVIATORIC RIGHT CAUCHY-GREEN TENSOR
C     VORIF  - FIBER ORIENTATION IN UNDERFORMED CONFIGURATION (FAMILY 1)
C     VORIF2 - FIBER ORIENTATION IN UNDERFORMED CONFIGURATION (FAMILY 2)
C     M0,N0  - STRUCTURAL FIBERS TENSOR IN UNDEFORMED  CONFIGURATION
C     M,MM,N,NN - STRUCTURAL FIBERS TENSOR IN DEFORMED  CONFIGURATION
C     MM0,NN0 -FIBERS 4ORDER TENSOR IN UNDEFORMED CONFIGURATION
C     VD     - FIBER ORIENTATION IN DERFORMED CONFIGURATION (FAMILY 1)
C     VD2    - FIBER ORIENTATION IN DERFORMED CONFIGURATION (FAMILY 2)
C     SSX    - AUXILIAR TENSORS
C     SX     - PIOLA KIRCHHOFF STRESSES OF EACH CONTRIBUTION: SVOL,SISOM,SISOF4
C     SIGMA  - TOTAL CAUCHY STRESS
C     CCMAT  - MATERIAL ELASTICITY TENSOR
C     CCMATE  - ELASTIC MATERIAL ELASTICITY TENSOR
C     CCSPATIAL - SPATIAL ELASTICITY TENSOR
C     DDSDDE - JACOBIAN MATERIAL CONTRIBUTIONS
C---------------------------------------------------------------------
C     LOCAL VARIABLES
C---------------------------------------------------------------------
C     JAC       - DETERMINANT OF THE TOTAL DEFORMATION GRADIENT
C     JACE      - DETERMINANT OF THE ELASTIC DEFORMATION GRADIENT
C     I1BARE     - 1ST INVARIANT OF CBARE
C     I4BAR     - SQUARE OF THE FIBER1 STRETCH
C     STRETCH4B - DEVIATORIC STRETCH OF FIBER 1
C     STRETCH4  - STRETCH OF FIBER 1
C     I6BAR     - SQUARE OF THE FIBER2 STRETCH
C     STRETCH6B - DEVIATORIC STRETCH OF FIBER 2
C     STRETCH6  - STRETCH OF FIBER 2
C     DUDXX     - 1ST DERIVATIVE OF PHI RESPECTIVE TO XX
C     D2UDXXX - 2ND DERIVATIVES OF PHI RESPECTIVE TO XXX
C---------------------------------------------------------------------
C                COUNTERS AND STATE INDICATORS
       INTEGER  II1(6),II2(6),ISTAT,INOEL,I,J,K,I1,J1,L,
     1         PSFLAG,ISTATC,MAXIT
C
C                KINEMATIC ARRAY VARIABLES
       DOUBLE PRECISION  UNIT(3,3),CE(3,3),CBARE(3,3),
     1        DISTGRE(3,3),DFGRD1TE(3,3),DISTGRTE(3,3),CINVE(3,3),DETCE,
     1        CBARINVE(3,3),DETCBARE,DFGRD1E(3,3),
     1        DFGRD1T(3,3),C(3,3),CINV(3,3),DETC,CBAR(3,3),
     1        CBARINV(3,3),DISTGRT(3,3),DISTGR(3,3),
     1        DETCBAR,DFGRD1M(3,3),SCALE,
C                FIBERS STRUCTURE VARIABLES
     2        VORIF(3),VD(3),M(3,3),M0(3,3),
     3        VORIF2(3),VD2(3),N(3,3),N0(3,3),MM0(3,3,3,3),
     4        NN0(3,3,3,3),S(3,3),
C                STRESS VARIABLES
     5        SBAR(3,3),SBARM(3,3),SBARF4(3,3),SBARF6(3,3),SISO(3,3),
     6        SISOM(3,3),SISOF4(3,3),SISOF6(3,3),SVOL(3,3),
     7        SIGMA(3,3),JAC,JACE,JACM,SE(3,3),SAUX(3,3),
     8        SISOM0(3,3),SISOF40(3,3),SISOF60(3,3),
     9        PK(3,3),SISO0(3,3),
C                MATERIAL TANGENT VARIABLES
     8        PP(3,3,3,3),PPMOD(3,3,3,3),PPT(3,3,3,3),CCISO(3,3,3,3),
     9        CCISO1(3,3,3,3),CCISO2(3,3,3,3),CCISO3(3,3,3,3),
     9        CCVOL(3,3,3,3),CCMAT(3,3,3,3),CCSPATIAL(3,3,3,3),
     9        DDSIGDDE(3,3,3,3),ESC,AAUX(3,3,3,3),
     9        DDSDDEJR(3,3,3,3),CCBAR(3,3,3,3),AA(3,3,3,3),
     9        A1(3,3,3,3),A2(3,3,3,3),A3(3,3,3,3),B1(3,3,3,3),
     9        B2(3,3,3,3),BB(3,3,3,3),AC(3,3,3,3),A3C(3,3,3,3),
     9        A4(3,3,3,3),PPP(3,3,3,3),DSDC2(3,3,3,3),DSDFP4(3,3,3,3),
     9        DSTPDC(3,3),E1(3,3),CCMATE(3,3,3,3),CCMATE2I(3,3,3,3),
     9        E2(3,3,3,3),DFP4DSTP(3,3),
     9        SS0(3,3,3,3)
C
C                MATERIAL PARAMETERS VARIABLES
          DOUBLE PRECISION  D1,C10,K11,K12,K21,K22,
     1          BETA,STPCRIT4,SHP4,BETA0,
C                KINEMATIC SCALAR VARIABLES
     2        VAR0,VAR1,VAR2,VAR3,VAR4,
     3        I1BARE,PP1,PP2,A,STRETCH4B,STRETCH6B,
     3        STRETCH6B2, SUM1, SUM2,STE4,STE42,STE4INV,
     5        DNORM,DNORM2,ST4INV,ST4INV2, ST6INV,ST6INV2,
     6        I4BAR,I6BAR,J23,J43,I1BAR,STRETCH4B2,STE,STE2,
C               DAMAGE SCALAR VARIABLES
     4        ECC,ECC0,DMG,DMG0,CD,MAXACC,MAXD,
     4        GDER0,GDER,SEF,SE0,SEF1,SEF2,SEM,
C               PERMANENT SET VARIABLES
     1         STP,STP0,PHI,FP4(3,3),FP4INV(3,3),FP4INVT(3,3),
     1         STP4INV,TOL,DSTPDC1(3,3),DCBARDC(3,3,3,3),St(3,3),
C                STRAIN-ENERGY DERIVATIVES VARIABLES
     6        DUDI1, DUDJ, DUDST4, DUDST6, D2UDST4, D2UDST6, D2DJ
C
C
C----------------------------------------------------------------------
C     MATERIAL CONSTANTS
C----------------------------------------------------------------------
       D1 = PROPS(1)
       C10 = PROPS(2)
       K11  = PROPS(3)
       K12  = PROPS(4)
       K21  = PROPS(5)
       K22  = PROPS(6)
C     DAMAGE
       CD      = PROPS(7)
       MAXACC  = PROPS(8)
       MAXD   = PROPS(9)
C     PERMANENT SET
       STPCRIT4 = PROPS(10)
       SHP4 = PROPS(11)
       BETA = PROPS(12)
       PSFLAG = PROPS(13)
       TOL = PROPS(14)
       MAXIT = PROPS(15)
       BETA0 = PROPS(16)
C
c
C----------------------------------------------------------------------
C     INITIALIZATIONS
C----------------------------------------------------------------------
C     AUXILIAR VARIABLES
       JAC=ZERO
       JACE=ZERO
       I1BARE=ZERO
       I1BAR=ZERO
       CE=ZERO
       CBARE=ZERO
       CINVE=ZERO
       CBARINVE=ZERO
       C=ZERO
       CBAR=ZERO
       CINV=ZERO
       CBARINV=ZERO
       DUDJ=ZERO
       D2DJ=ZERO
       DUDI1=ZERO
C
       I4BAR=ZERO
       I6BAR=ZERO
       STRETCH4B=ONE
       STRETCH4B2=ONE
       STE4=ONE
       STE42=ONE
       STRETCH6B=ONE
       STRETCH6B2=ONE
       VAR0=ZERO
       VAR1=ZERO
       VAR2=ZERO
       DUDST4=ZERO
       D2UDST4=ZERO
       VAR3=ZERO
       VAR4=ZERO
       DUDST6=ZERO
       D2UDST6=ZERO
       ST4INV=ZERO
       ST6INV=ZERO
       J23=ZERO
       J43=ZERO
       ST4INV2=ZERO
       ST6INV2=ZERO
C
       ISTAT=1
       PP1=ZERO
       PP2=ZERO
C
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C     IDENTITY TENSOR DEFINITION                                       *
C----------------------------------------------------------------------
C
       CAll ONEM(UNIT)
C     INICIALIZAÇAO STATEV
      IF (TIME(1).EQ.0.D0.AND.KSTEP.EQ.1) THEN
        STATEV(3)=ONE
        STATEV(4)=ONE
        STATEV(5)=ONE
        STATEV(6)=ONE
        DO I=7,22
         STATEV(I)=ZERO
        END DO

        END IF
C
C----------------------------------------------------------------------
C     POINTERS DEFINITION
C----------------------------------------------------------------------
C     POINTERS TO STORAGE STRESS AND DDSDDE ARRAYS
       II1(1)=1
       II1(2)=2
       II1(3)=3
       II1(4)=1
       II1(5)=1
       II1(6)=2
C
       II2(1)=1
       II2(2)=2
       II2(3)=3
       II2(4)=2
       II2(5)=3
       II2(6)=3
C
C
C
C----------------------------------------------------------------------
C     UNIT VECTOR IN THE DIRECTION OF THE UNDEFORMED FIBERS
C----------------------------------------------------------------------
            CALL FIB_VEC(NOEL,NELEM,INOEL,FIBORI,FIBORI2,NDI,
     1         VORIF,VORIF2,M0,N0)
C----------------------------------------------------------------------
C     UNIT VECTOR IN THE DIRECTION OF THE DEFORMED FIBERS
C----------------------------------------------------------------------
C
       DO I=1,NDI
          SUM1=ZERO
          SUM2=ZERO
         DO J=1,NDI
          SUM1=SUM1+DFGRD1(I,J)*VORIF(J)
          SUM2=SUM2+DFGRD1(I,J)*VORIF2(J)
         ENDDO
C
C     FIBER DIRECTIONS IN THE DEFORMED CONFIGURATION
C               -FAMILY 1
         VD(I)=SUM1
C
C               -FAMILY 2
         VD2(I)=SUM2
       ENDDO
       DNORM=DSQRT(VD(1)*VD(1)+
     1             VD(2)*VD(2)+
     2             VD(3)*VD(3))
       DNORM2=DSQRT(VD2(1)*VD2(1)+
     1             VD2(2)*VD2(2)+
     2             VD2(3)*VD2(3))
C           COSINE OF THE ANGLE BETWEEN FIBERS
      DO I=1,NDI
       VD(I)=VD(I)/DNORM
       VD2(I)=VD2(I)/DNORM2
       A=A+VD(I)*VD2(I)
      END DO
C
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C      TOTAL DEFORMATION GRADIENT AND TOTAL JACOBIAN
C----------------------------------------------------------------------
C
      JAC=ZERO
      JAC = DFGRD1(1,1) * DFGRD1(2,2) * DFGRD1(3,3)
     1    - DFGRD1(1,2) * DFGRD1(2,1) * DFGRD1(3,3)
C
      IF (NSHR .EQ. 3) THEN
          JAC = JAC + DFGRD1(1,2) * DFGRD1(2,3) * DFGRD1(3,1)
     1              + DFGRD1(1,3) * DFGRD1(3,2) * DFGRD1(2,1)
     2              - DFGRD1(1,3) * DFGRD1(3,1) * DFGRD1(2,2)
     3              - DFGRD1(2,3) * DFGRD1(3,2) * DFGRD1(1,1)
      END IF
C
      DFGRD1T=TRANSPOSE(DFGRD1)
      CALL M3MULT(DFGRD1T,DFGRD1,C)
      CALL MATINV3D(C,CINV,DETC,ISTATC)
C
C     AUXILIAR VARIABLES
      SCALE=JAC**(-ONE/THREE)
C
      DISTGR=SCALE*DFGRD1
c
      CBAR=(JAC**(-TWO/THREE))*C
      CALL TRACE(CBAR,I1BAR)
C
        I4BAR=ZERO
        I6BAR=ZERO
      DO I=1,NDI
        DO J=1, NDI
            I4BAR=I4BAR+CBAR(I,J)*M0(I,J)
            I6BAR=I6BAR+CBAR(I,J)*N0(I,J)
        ENDDO
      ENDDO
C
C     FAMILY 1
      STRETCH4B=DSQRT(I4BAR)
C     FAMILY 2
      STRETCH6B=DSQRT(I6BAR)
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C     GROWTH AND DAMAGE INITIALIZATION
C----------------------------------------------------------------------
C
      SE0=STATEV(10)
      ECC0=STATEV(11)
      DMG0=STATEV(12)
      GDER0=STATEV(13)
      STP0=STATEV(5)
C
C
C----------------------------------------------------------------------
C     GROWTH AND DAMAGE EVOLUTION
C----------------------------------------------------------------------
C
c
       CALL PS_EVOL(PROPS,ECC0,DMG0,STP0,SE0,DTIME,M0,N0,DFGRD1,
     1    NSHR,JAC,STRETCH4B,STRETCH6B,STP,DSTPDC,DMG,ECC,SEF,GDER,
     2    GDER0)
C
C
C----------------------------------------------------------------------
C     UPDATE PERMANENT SET DEFORMATION GRADIENT
C----------------------------------------------------------------------
      STP4INV=STP**(-ONE)
      FP4=UNIT+(STP-ONE)*M0
      FP4INV=UNIT-((STP-ONE)*STP4INV)*M0
      FP4INVT=TRANSPOSE(FP4INV)
C
       CALL ELASTIC_TENSORS(DFGRD1,FP4INV,DFGRD1E,NSHR,
     1      CE,CBARE,I1BARE,JAC,JACE,CINVE)
C
C     SEF SEM DANO E COM PS
      SEM=C10*(I1BARE-THREE)
C     ANISOTROPIC PART
      STE=STRETCH4B/STP
      STE2=STE**TWO
      SEF1=(K11/(TWO*K12))*(DEXP(K12*(STE2-ONE)**TWO)-ONE)
      STRETCH6B2=STRETCH6B**TWO
      SEF2=(K21/(TWO*K22))*(DEXP(K22*(STRETCH6B2-ONE)**TWO)-ONE)
      SEF=SEM+SEF1+SEF2
c
      STRETCH6B2=STRETCH6B**TWO
C     UPDATE STE4
      STE4=STRETCH4B*STP4INV
C     UPDATE VARIABLES
      STATEV(10)=SEF
C     DEFINE ACCUMULATED SEF VARIABLES FROM THE PREVIOUS STATE
      STATEV(11)=ECC
C     DEFINE DAMAGE VARIABLES FROM THE PREVIOUS STATE
      STATEV(12)=DMG
      STATEV(13)=GDER
C----------------------------------------------------------------------
C     STRAIN-ENERGY DERIVATIVES WITH RESPECT TO ELASTIC INVARIANTS
C----------------------------------------------------------------------
C     FIRST AND SECOND DERIVATIVES PHI - VOLUME-PRESERVING PART
      DUDJ=(ONE/D1)*(JAC-ONE/JAC)
      D2DJ=(ONE/D1)*(ONE+ONE/(JAC**TWO))
C----------------------------------------------------------------------
C     FIRST AND SECOND DERIVATIVES PHI
C     MATRIX CONTRIBUTION
      DUDI1=C10
C     FIBERS CONTRIBUTION
C     FAMILY 1
      STE42=STE4**TWO
      VAR0=STE42-ONE
      VAR1=VAR0**TWO
      IF (STE4.GE.ONE) THEN
C
      DUDST4=TWO*K11*VAR0*STE4*EXP(K12*VAR1)
C
      D2UDST4=FOUR*K11*STE42*EXP(K12*VAR1)+
     +       TWO*K11*VAR0*EXP(K12*VAR1)+
     +       (8.D0)*K11*VAR1*STE42*
     *        K12*EXP(K12*VAR1)
      ELSE
      DUDST4=ZERO
      D2UDST4=ZERO
      END IF
c
C     FAMILY 2
      STRETCH6B2=STRETCH6B**TWO
      VAR3=(STRETCH6B2-ONE)
      VAR4=VAR3**TWO
      IF(I6BAR.GE.ONE) THEN
C
      DUDST6=TWO*K21*VAR3*STRETCH6B*EXP(K22*VAR4)
C
      D2UDST6=FOUR*K21*STRETCH6B2*EXP(K22*VAR4)+
     +       TWO*K21*VAR3*EXP(K22*VAR4)+
     +       (8.D0)*K21*VAR4*STRETCH6B2*K22*EXP(K22*VAR4)
      ELSE
      DUDST6=ZERO
      D2UDST6=ZERO
      END IF
c
C----------------------------------------------------------------------
C     SECOND PIOLA-KIRCHHOFF STRESS TENSOR
C----------------------------------------------------------------------
C     AUXILIAR VARIABLES
      STE4INV=(STE4)**(-ONE)
      ST6INV=(STRETCH6B)**(-ONE)
      J23=JAC**(-TWO/THREE)
C
C     ISOTROPIC CONTRIBUTION (SBAR)
      SBARM=TWO*DUDI1*UNIT
      SBARF4=DUDST4*(STE4INV)*M0
      SBARF6=DUDST6*ST6INV*N0
      SBAR=SBARM+SBARF4+SBARF6
C     ISOTROPIC CONTRIBUTION (SISO)
      SISOM=DUDI1*(TWO*J23*UNIT-(TWO/THREE)*I1BARE*CINVE)
      SISOF4=DUDST4*(J23*STE4INV*M0-(ONE/THREE)*STE4*CINVE)
      SISOF6=DUDST6*(J23*ST6INV*N0-(ONE/THREE)*STRETCH6B*CINVE)
C
C     UNDAMAGED SECOND PIOLA-KIRCHHOFF STRESS TENSOR (S)
      SISOM0=SISOM
      SISOF40=SISOF4
      SISOF60=SISOF6
      SISO0=SISOM0+SISOF40+SISOF60
C     SECOND PIOLA-KIRCHHOFF STRESS TENSOR (S)
      SISO=(SISOM+SISOF4+SISOF6)
C     VOLUMETRIC CONTRIBUTION (SVOL)
      SVOL=JAC*DUDJ*CINV
C
C     ELASTIC SECOND PIOLA-KIRCHHOFF STRESS TENSOR
      SE=SISO
C     TOTAL SECOND PIOLA-KIRCHHOFF STRESSES
      CALL M3MULT(FP4INV,SE,SAUX)
      CALL M3MULT(SAUX,FP4INVT,St)
      S=St+SVOL
c
C----------------------------------------------------------------------
C     CAUCHY STRESS TENSOR
C----------------------------------------------------------------------
      CALL PUSH2(S,DFGRD1,JAC,SIGMA,NDI)
C
C      FIRST PIOLA KIRCHHOFF STRESS TENSOR
       CALL M3MULT(DFGRD1,S,PK)
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C     JACOBIAN ELASTIC MATERIAL MATRIX
C----------------------------------------------------------------------
C     AUXILIAR VARIABLES
c
      J43=JAC**(-FOUR/THREE)
      ST4INV2=STE4**(-TWO)
      ST6INV2=STRETCH6B**(-TWO)
C
C     MATERIAL MATRIX(CCBAR)
      CALL TENSORPROD2(M0,M0,MM0,NDI)
      CALL TENSORPROD2(N0,N0,NN0,NDI)
      DO I=1,NDI
       DO J=1,NDI
         DO K=1,NDI
          DO L=1,NDI
          CCBAR(I,J,K,L)=J43*ST4INV2*(D2UDST4-STE4INV*DUDST4)*
     *    MM0(I,J,K,L)+J43*ST6INV2*(D2UDST6-ST6INV*DUDST6)*NN0(I,J,K,L)
          END DO
         END DO
       END DO
      END DO
c
C     MATERIAL MATRIX(CCISO)
      CALL PROJECTION(CE,CINVE,NDI,PP)
      CALL PROJECTIONMOD(CINVE,NDI,PPMOD)
      CALL PROJECTIONT(CE,CINVE,NDI,PPT)
C
      CALL DOUBLEDOT4(PP,CCBAR,AAUX,NDI)
      CALL DOUBLEDOT4(AAUX,PPT,CCISO1,NDI)
C
      CALL DOUBLEDOT2(SBAR,CE,ESC,NDI)
C
      CALL TENSORPROD2(CINVE,SISO0,CCISO2,NDI)
      CALL TENSORPROD2(SISO0,CINVE,CCISO3,NDI)
      CALL TENSORPROD2(SISO0,SISO0,SS0,NDI)
C
      CCISO=CCISO1+(TWO/THREE)*J23*ESC*PPMOD-(TWO/THREE)*(CCISO2+CCISO3)
c
C
c
C     VOLUMETRIC CONTRIBUTION
      DO I=1,NDI
       DO J=1,NDI
         DO K=1,NDI
          DO L=1,NDI
          CCVOL(I,J,K,L)=(TWO/D1)*(JAC**TWO)*(CINV(I,J)*CINV(K,L))-
     -     (ONE/D1)*((JAC**TWO)-ONE)*(CINV(I,K)*CINV(J,L)+CINV(I,L)*
     *     CINV(J,K))
          END DO
         END DO
       END DO
      END DO
c
c
      CCMATE=CCISO
c
C
C----------------------------------------------------------------------
C     JACOBIAN MATERIAL MATRIX
C     CCMAT=2*DSDC+2[DSDFP4:DFP4DSTP](X)DSTPDC   (X) TENSORPROD
C----------------------------------------------------------------------
C     DSDFP4-----------------------------------------------------------
      IF (PSFLAG.EQ.1) THEN
      CCMATE2I=0.5D0*CCMATE
      CALL TENSORPRODBARUP(FP4INV,St,A1,NDI)
      CALL TENSORPRODBARDOWN(St,FP4INV,A2,NDI)
      CALL TENSORPRODBARUP(FP4INV,FP4INV,A3,NDI)
      AA=-(A1+A2)-A3
      CALL DOUBLEDOT4(AA,CCMATE2I,AC,NDI)
      CALL TENSORPRODBARDOWN(FP4INVT,CE,B1,NDI)
c
      CALL TENSORPRODBARUP(CE,FP4INVT,B2,NDI)
      BB=B1+B2
      CALL DOUBLEDOT4(AC,BB,DSDFP4,NDI)
C     DSDC2=2*DSDC------------------------------------------------------
      CALL DOUBLEDOT4(A3,CCMATE,A3C,NDI)
      CALL TENSORPRODBARUP(FP4INVT,FP4INVT,A4,NDI)
      CALL DOUBLEDOT4(A3C,A4,DSDC2,NDI)
C     DFG4DSTG4---------------------------------------------------------
      DFP4DSTP=M0
C     DSTPsDC-----------------------------------------------------------
C     OBTAINED IN PS_EVOL
C     [DSDFP4:DFP4DSTP]------------------------------------------------
      CALL DOUBLEDOT42(DSDFP4,DFP4DSTP,E1,NDI)
      E1=TWO*E1
C     2[DSDFG4:DFP4DSTP](X)DSTG4DC
c
      CALL TENSORPROD2(E1,DSTPDC,E2,NDI)
C
       CCMAT=DSDC2+E2
       ELSE
       CCMAT=CCMATE
       END IF
C
      CCMAT=CCMAT+CCVOL
C----------------------------------------------------------------------
      CALL PUSH4(CCMAT,DFGRD1,JAC,CCSPATIAL,NDI)
C
c
C     SPATIAL TANGENT MODULI BASED ON THE JAUMANN
C                        RATE OF THE KIRCHHOFF STRESS TENSOR
      DO I=1,NDI
         DO J=1,NDI
            DO K=1,NDI
               DO L=1,NDI
C         ------JAUMMAN RATE PART------
                   DDSDDEJR(I,J,K,L)=
     +             (ONE/TWO)*(UNIT(I,K)*SIGMA(J,L)
     +             +SIGMA(I,K)*UNIT(J,L)+UNIT(I,L)*SIGMA(J,K)
     +             +SIGMA(I,L)*UNIT(J,K))
C        -----SPATIAL TANGENT MODULI------
                   DDSIGDDE(I,J,K,L)=DDSDDEJR(I,J,K,L)+
     +                               CCSPATIAL(I,J,K,L)
               END DO
             END DO
          END DO
       END DO
c
C
C
C----------------------------------------------------------------------
C     STRESS AND JACOBIAN MATRIX STORAGE
C----------------------------------------------------------------------
c
      DO I1=1,NTENS
C       STRESS VECTOR
         STRESS(I1)=SIGMA(II1(I1),II2(I1))
         DO J1=1,NTENS
C       DDSDDE - FULLY SIMMETRY IMPOSED
            PP1=DDSIGDDE(II1(I1),II2(I1),II1(J1),II2(J1))
            PP2=DDSIGDDE(II1(I1),II2(I1),II2(J1),II1(J1))
            DDSDDE(I1,J1)=(ONE/TWO)*(PP1+PP2)
         END DO
      END DO
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C     STATE VARIABLES
C----------------------------------------------------------------------
        STATEV(1) = JACE
        STATEV(2) = JAC
        STATEV(3) = STRETCH4B
        STATEV(4) = STRETCH6B
        STATEV(5) = STP
        STATEV(6) = STE4
        STATEV(7) = VD(1)
        STATEV(8) = VD(2)
        STATEV(9) = VD(3)
        STATEV(14)= PK(1,1)
        STATEV(15)= PK(1,2)
        STATEV(16)= PK(1,2)
        STATEV(17)= PK(2,1)
        STATEV(18)= PK(2,2)
        STATEV(19)= PK(2,3)
        STATEV(20)= PK(3,1)
        STATEV(21)= PK(3,2)
        STATEV(22)= PK(3,3)
C
C----------------------------------------------------------------------
      RETURN
      END
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C END OF MAIN UMAT ROUTINE
C****************************************************************************
C     UTILITY SUBROUTINES
C****************************************************************************
C
C
      SUBROUTINE FIB_VEC(NOEL,NELEM,INOEL,FIBORI,FIBORI2,NDI,VORIF,
     1                   VORIF2,M0,N0)
C
       Implicit None
C
       INTEGER I,J,NDI,INOEL,NELEM,NOEL
C
       DOUBLE PRECISION FIBORI(NELEM,4), FIBORI2(NELEM,4),DNORM,DNORM2,
     1    VORIF(3),VORIF2(3), M0(3,3), N0(3,3)
        INOEL=0
        I=0
        DO I=1,NELEM
C               ELEMENT IDENTIFICATION
            IF(NOEL.EQ.INT(FIBORI(I,1))) THEN
                INOEL=I
            ENDIF
        ENDDO
C
C     FIBORI - FIBER ORIENTATION - FAMILY 1
             DNORM=DSQRT(FIBORI(INOEL,2)*FIBORI(INOEL,2)+
     1                   FIBORI(INOEL,3)*FIBORI(INOEL,3)+
     2                   FIBORI(INOEL,4)*FIBORI(INOEL,4))
C
C      FIBORI2 - FIBER ORIENTATION - FAMILY 2
             DNORM2=DSQRT(FIBORI2(INOEL,2)*FIBORI2(INOEL,2)+
     1                   FIBORI2(INOEL,3)*FIBORI2(INOEL,3)+
     2                   FIBORI2(INOEL,4)*FIBORI2(INOEL,4))
C
C       UNDERFORMED FIBER ORIENTATION TENSOR
C
        DO I=1,NDI
        J=I+1
C       FIBER ORIENTATION NORMALIZED VECTOR - FAMILY 1
        VORIF(I)=FIBORI(INOEL,J)/DNORM
C       FIBER ORIENTATION NORMALIZED VECTOR - FAMILY 2
        VORIF2(I)=FIBORI2(INOEL,J)/DNORM2
        END DO
C
      DO I=1,NDI
       DO J=1,NDI
C       STRUCTURAL TENSOR - FAMILY 1
       M0(I,J)=VORIF(I)*VORIF(J)
C       STRUCTURAL TENSOR - FAMILY 2
       N0(I,J)=VORIF2(I)*VORIF2(J)
       END DO
      END DO
c
      RETURN
      END SUBROUTINE FIB_VEC
C-------------------------------------------------------------------------------
       SUBROUTINE DOUBLEDOT2(A,B,C,NDI)
C
       Implicit None
C
       INTEGER I,J,NDI
C
       DOUBLE PRECISION A(3,3),B(3,3),C
C
      C=0.D0
      DO I=1,NDI
       DO J=1,NDI
         C=C+A(I,J)*B(J,I)
       END DO
      END DO
C
      RETURN
C
      end SUBROUTINE DOUBLEDOT2
C-------------------------------------------------------------------------------
       SUBROUTINE DOUBLEDOT42(A,B,C,NDI)
C
       Implicit None
C
       INTEGER I,J,K,L,NDI
C
       DOUBLE PRECISION A(3,3,3,3),B(3,3),C(3,3),AUX
C
      C=0.D0
      DO I=1,NDI
       DO J=1,NDI
         AUX=0.D0
         DO K=1,NDI
           DO L=1,NDI
             AUX=AUX+A(I,J,K,L)*B(K,L)
           END DO
         END DO
         C(I,J)=AUX
       END DO
      END DO
C
      RETURN
C
      end SUBROUTINE DOUBLEDOT42
C-------------------------------------------------------------------------------
       SUBROUTINE DOUBLEDOT4(A,B,C,NDI)
C
       Implicit None
C
       INTEGER I,J,K,L,NDI, M, N
C
       DOUBLE PRECISION A(3,3,3,3),B(3,3,3,3),C(3,3,3,3),SUM
C
      DO I=1,NDI
       DO J=1,NDI
         DO K=1,NDI
          DO L=1,NDI
           SUM=0.D0
           DO M=1,NDI
            DO N=1,NDI
             SUM=SUM+A(I,J,M,N)*B(M,N,K,L)
            END DO
           END DO
           C(I,J,K,L)=SUM
          END DO
         END DO
       END DO
      END DO
C
      RETURN
C
      end SUBROUTINE DOUBLEDOT4
C-------------------------------------------------------------------------------
       SUBROUTINE PS_EVOL(PROPS,ECC0,DMG0,STP0,SE0,DTIME,M0,
     1    N0,DFGRD1,NSHR,JAC,STRETCH4B,STRETCH6B,STP,DSTPDC,DMG,ECC,
     2    SE,GDER,GDER0)
C
       IMPLICIT NONE
       INCLUDE 'PARAM_UMAT.F'
C
       INTEGER MAXIT,IT,NSHR,PSFLAG
C
       DOUBLE PRECISION TOL,PHI_IT,STRETCH4B,STRETCH6B,STP,
     3     ECC0,DMG0,STP0,PHID_IT,DPHID,
     3  SE0,SE,DTIME,DFGRD1(3,3),DSTPDC(3,3),DMG,ECC,GDER,
     2  PROPS(13),M0(3,3),N0(3,3),C10,K11,K12,K21,K22,
     2  CD,MAXACC,MAXD,PSCRIT,SHP,BETA,STPN1,STPINV(3,3),FP(3,3),
     2  FPT(3,3),FPINV(3,3),FPINVT(3,3),
     3  UNIT(3,3),R,KP,DPHI,DKDD,DDDACC,STRETCH4B2,STP2,STP5,
     3  DSEFFDSTP,DFPTDSTP(3,3),A1(3,3),A2(3,3),A3(3,3),A4(3,3),
     3  CBARE(3,3),DCEDSTP(3,3),DSEFMDSTP,DSEFDSTP,DKP,KK,STPNEW,DMGIT,
     3  ECCIT,GDERIT,SEIT,STPIT,DACCDST,DFGRD1E(3,3),ZD,
     3  DIF,GDER0,I1BARE,JAC,JACE,STE2,DIF_D,
     3  CE(3,3),DFPDSTP(3,3),A5,DSTDC(3,3),CINVE(3,3),
     3  SEM, SEF1,SEF2,STE,STRETCH6B2,DSTPDST,BETA0,ZDIT
C
      C10 = PROPS(2)
      K11  = PROPS(3)
      K12  = PROPS(4)
      K21  = PROPS(5)
      K22  = PROPS(6)
C     DAMAGE
      CD       = PROPS(7)
      MAXACC  = PROPS(8)
      MAXD    = PROPS(9)
C     PERMANENT SET
      PSCRIT = PROPS(10)
      SHP = PROPS(11)
      BETA = PROPS(12)
      PSFLAG = PROPS(13)
      TOL=PROPS(14)
      MAXIT=PROPS(15)
      BETA0=PROPS(16)
C
      CAll ONEM(UNIT)
C
         STPINV=STP0**(-ONE)
         FP=UNIT+(STP0-ONE)*M0
         FPT=TRANSPOSE(FP)
         FPINV=UNIT-((STP0-ONE)*STPINV)*M0
         FPINVT=TRANSPOSE(FPINV)
c
c
            CALL ELASTIC_TENSORS(DFGRD1,FPINV,DFGRD1E,NSHR,
     1    CE,CBARE,I1BARE,JAC,JACE,CINVE)
c
c
C     INITIAL STRAIN-ENERGY
      SEM=C10*(I1BARE-THREE)
C     ANISOTROPIC PART
      STE=STRETCH4B/STP0
      STE2=STE**TWO
      SEF1=(K11/(TWO*K12))*(DEXP(K12*(STE2-ONE)**TWO)-ONE)
      STRETCH6B2=STRETCH6B**TWO
      SEF2=(K21/(TWO*K22))*(DEXP(K22*(STRETCH6B2-ONE)**TWO)-ONE)
      SE=SEM+SEF1+SEF2
C
C
      IF (PSFLAG.EQ.1) THEN
      IT=0
      ZD=SE-SE0
C
      DO WHILE (IT.LT.MAXIT)
C
      IT=IT+1
c
c
      IF (IT.EQ.1) THEN
          STPN1=STP0
          CALL DMG_EVOL(SE0,ECC0,DMG0,GDER0,STPN1,CD,STRETCH4B,
     1   STRETCH6B,I1BARE,C10,K11,K12,K21,K22,MAXACC,MAXD,GDERIT,
     1   DMG,ECC,SE)
         ZD=SE-SE0
         ELSE
          STPN1=STPIT
          DMG=DMGIT
          ZD=ZDIT

         END IF
         PHI_IT=(STRETCH4B/STPN1)-PSCRIT
C
            PHID_IT=DMG
c
       IF ((PHI_IT.GE.ZERO).AND.(ZD.GT.ZERO).AND.(PHID_IT.GT.ZERO)) THEN
C
         STPINV=STPN1**(-ONE)
         FP=UNIT+(STPN1-ONE)*M0
         FPT=TRANSPOSE(FP)
         FPINV=UNIT-((STPN1-ONE)*STPINV)*M0
         FPINVT=TRANSPOSE(FPINV)
C
      CALL ELASTIC_TENSORS(DFGRD1,FPINV,DFGRD1E,NSHR,
     1      CE,CBARE,I1BARE,JAC,JACE,CINVE)
c
c
       PHI_IT=(STRETCH4B/STPN1)-PSCRIT
c
c
       KP=(BETA0+(BETA*DMG))**SHP
c
C
       R=STPN1-STP0-KP*PHI_IT*DTIME
C
       DPHI=(-STRETCH4B)*(STPN1**(-TWO))
C
       DKDD=BETA*SHP*((DMG*BETA+BETA0)**(SHP-ONE))
       IF (DMG.GE.ONE) THEN
            ECC=0.9999*MAXACC
       END IF
       DDDACC=(MAXD/(MAXACC*CD))*((ONE-(ECC/MAXACC))**((ONE/CD)-ONE))
C
       STRETCH4B2=STRETCH4B*STRETCH4B
       STP2=STPN1*STPN1
       STP5=STPN1**(-5.0D0)
       DSEFFDSTP=TWO*K11*STRETCH4B2*STP5*(STP2-STRETCH4B2)*
     *    DEXP(K12*(((STRETCH4B2/STP2)-ONE)**TWO))

       DFPTDSTP=TRANSPOSE(M0)
       DFPDSTP=M0
       CALL M3MULT(FPINVT,DFPTDSTP,A1)
       CALL M3MULT(A1,CBARE,A2)
       CALL M3MULT(CBARE,DFPDSTP,A3)
       CALL M3MULT(A3,FPINV,A4)
       DCEDSTP=-A2-A4
       CALL DOUBLEDOT2(UNIT,DCEDSTP,A5,3)
       DSEFMDSTP=C10*A5

       DSEFDSTP=DSEFFDSTP+DSEFMDSTP
c
       DKP=DKDD*DDDACC*DSEFDSTP
c
c
C
       KK=ONE-(DKP*PHI_IT+KP*DPHI)*DTIME
c
C
       STPNEW=STPN1-R*(KK**(-ONE))
c
c
         STPINV=STPNEW**(-ONE)
         FP=UNIT+(STPNEW-ONE)*M0
         FPT=TRANSPOSE(FP)
         FPINV=UNIT-((STPNEW-ONE)*STPINV)*M0
         FPINVT=TRANSPOSE(FPINV)
C
      CALL ELASTIC_TENSORS(DFGRD1,FPINV,DFGRD1E,NSHR,
     1      CE,CBARE,I1BARE,JAC,JACE,CINVE)

      CALL DMG_EVOL(SE0,ECC0,DMG0,GDER0,STPNEW,CD,STRETCH4B,
     1       STRETCH6B,I1BARE,C10,K11,K12,K21,K22,MAXACC,MAXD,GDERIT,
     1        DMGIT,ECCIT,SEIT)
c
      ZDIT=SEIT-SE0
c
      DIF=ABS(R)
      DIF_D=ABS(DMGIT-DMG)
C
      IF ((DIF.LE.TOL).AND.(DIF_D.LE.TOL)) THEN
       STP=STPNEW
       DMG=DMGIT
       ECC=ECCIT
       GDER=GDERIT
       SE=SEIT
C
      STPINV=STP**(-ONE)
      FP=UNIT+(STP-ONE)*M0
      FPT=TRANSPOSE(FP)
      FPINV=UNIT-((STP-ONE)*STPINV)*M0
      FPINVT=TRANSPOSE(FPINV)

      CALL ELASTIC_TENSORS(DFGRD1,FPINV,DFGRD1E,NSHR,
     1     CE,CBARE,I1BARE,JAC,JACE,CINVE)


      PHI_IT=(STRETCH4B/STP)-PSCRIT

      KP=(BETA0+(BETA*DMG))**SHP

      DPHI=(-STRETCH4B)*(STP**(-TWO))
      IF (DMG.GE.ONE) THEN
            ECC=0.9999*MAXACC
      END IF
      DKDD=BETA*SHP*((DMG*BETA+BETA0)**(SHP-ONE))
      DDDACC=(MAXD/(MAXACC*CD))*((ONE-(ECC/MAXACC))**((ONE/CD)-ONE))

      STRETCH4B2=STRETCH4B*STRETCH4B
      STP2=STP*STP
      STP5=STP**(-5.0D0)
      DSEFFDSTP=TWO*K11*STRETCH4B2*STP5*(STP2-STRETCH4B2)*
     *    DEXP(K12*(((STRETCH4B2/STP2)-ONE)**TWO))
      DFPTDSTP=TRANSPOSE(M0)
      DFPDSTP=M0
      CALL M3MULT(FPINVT,DFPTDSTP,A1)
      CALL M3MULT(A1,CBARE,A2)
      CALL M3MULT(CBARE,DFPDSTP,A3)
      CALL M3MULT(A3,FPINV,A4)
      DCEDSTP=-A2-A4
      CALL DOUBLEDOT2(UNIT,DCEDSTP,A5,3)
      DSEFMDSTP=C10*A5

      DSEFDSTP=DSEFFDSTP+DSEFMDSTP

      DKP=DKDD*DDDACC*DSEFDSTP

      KK=ONE-(DKP*PHI_IT+KP*DPHI)*DTIME
c
      DSTPDST=(KP*DTIME)/(KK*STP)

      DSTDC=((TWO*STRETCH4B)**(-ONE))*M0
c
      DSTPDC=DSTPDST*DSTDC
      EXIT
      ELSE
      STPIT=STPNEW
      END IF
c
C
      ELSE
C
      STP=STP0
C
       PHI_IT=(STRETCH4B/STP)-PSCRIT
       STPINV=STP**(-ONE)
       FP=UNIT+(STP-ONE)*M0
       FPT=TRANSPOSE(FP)
       FPINV=UNIT-((STP-ONE)*STPINV)*M0
       FPINVT=TRANSPOSE(FPINV)
c
C
      CALL ELASTIC_TENSORS(DFGRD1,FPINV,DFGRD1E,NSHR,
     1      CE,CBARE,I1BARE,JAC,JACE,CINVE)
c
      CALL DMG_EVOL(SE0,ECC0,DMG0,GDER0,STP,CD,STRETCH4B,
     1       STRETCH6B,I1BARE,C10,K11,K12,K21,K22,MAXACC,MAXD,GDER,
     1        DMG,ECC,SE)
c
        IF (DMG.GE.ONE) THEN
            ECC=0.9999*MAXACC
        END IF
c
       KP=(BETA0+(BETA*DMG))**SHP
c
C
       DPHI=(-STRETCH4B)*(STP**(-TWO))
C
       DKDD=BETA*SHP*((DMG*BETA+BETA0)**(SHP-ONE))
       DDDACC=(MAXD/(MAXACC*CD))*((ONE-(ECC/MAXACC))**((ONE/CD)-ONE))
C
       STRETCH4B2=STRETCH4B*STRETCH4B
       STP2=STP*STP
       STP5=STP**(-5.0D0)
       DSEFFDSTP=TWO*K11*STRETCH4B2*STP5*(STP2-STRETCH4B2)*
     *    DEXP(K12*(((STRETCH4B2/STP2)-ONE)**TWO))

       DFPTDSTP=TRANSPOSE(M0)
       DFPDSTP=M0
       CALL M3MULT(FPINVT,DFPTDSTP,A1)
       CALL M3MULT(A1,CBARE,A2)
       CALL M3MULT(CBARE,DFPDSTP,A3)
       CALL M3MULT(A3,FPINV,A4)
       DCEDSTP=-A2-A4
       CALL DOUBLEDOT2(UNIT,DCEDSTP,A5,3)
       DSEFMDSTP=C10*A5

       DSEFDSTP=DSEFFDSTP+DSEFMDSTP
c
       DKP=DKDD*DDDACC*DSEFDSTP
C
       KK=ONE-(DKP*PHI_IT+KP*DPHI)*DTIME
c
      DSTPDST=(KP*DTIME)/(KK*STP)

      DSTDC=((TWO*STRETCH4B)**(-ONE))*M0
c
      DSTPDC=DSTPDST*DSTDC
c
c
      EXIT
C
      END IF

      END DO
c
C
      IF (IT.GT.MAXIT) THEN
      WRITE(*,*) 'WARNING GROWTH EVOL: IT>MAXIT'
      STOP
      END IF
      ELSE
      STP=STP0
C
       STPINV=STP**(-ONE)
       FP=UNIT+(STP-ONE)*M0
       FPT=TRANSPOSE(FP)
       FPINV=UNIT-((STP-ONE)*STPINV)*M0
       FPINVT=TRANSPOSE(FPINV)
C
      CALL ELASTIC_TENSORS(DFGRD1,FPINV,DFGRD1E,NSHR,
     1      CE,CBARE,I1BARE,JAC,JACE,CINVE)
c
      CALL DMG_EVOL(SE0,ECC0,DMG0,GDER0,STP,CD,STRETCH4B,
     1       STRETCH6B,I1BARE,C10,K11,K12,K21,K22,MAXACC,MAXD,GDER,
     1        DMG,ECC,SE)
C
       DSTPDC=ZERO

      END IF
C
c
      RETURN
C
      END SUBROUTINE PS_EVOL
C-------------------------------------------------------------------------------
      SUBROUTINE M3MULT(A,B,C)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION A(3,3),B(3,3), C(3,3), SUM
      INTEGER I,J,K
C
      DO I=1,3
       DO J=1,3
        SUM=0.D0
          DO K=1,3
           SUM=SUM+A(I,K)*B(K,J)
          END DO
        C(I,J)=SUM
       END DO
      END DO
C
      RETURN
C
      END SUBROUTINE M3MULT
C-------------------------------------------------------------------------------
       SUBROUTINE PROJECTION(C,CINV,NDI,PP)
C
       Implicit None
C
       INTEGER I,J,K,L,NDI
C
       DOUBLE PRECISION C(3,3),CINV(3,3),PP(3,3,3,3)
C
      DO I=1,NDI
       DO J=1,NDI
         DO K=1,NDI
          DO L=1,NDI
          IF ((I.EQ.K).AND.(J.EQ.L)) THEN
          PP(I,J,K,L)=1.D0-((1.D0)/(3.D0))*CINV(I,J)*C(K,L)
          ELSE
          PP(I,J,K,L)=-((1.D0)/(3.D0))*CINV(I,J)*C(K,L)
          END IF
          END DO
         END DO
       END DO
      END DO
C
      RETURN
C
      end SUBROUTINE PROJECTION
C-------------------------------------------------------------------------------
       SUBROUTINE PROJECTIONMOD(CINV,NDI,PPMOD)
C
       Implicit None
C
       INTEGER I,J,K,L,NDI
C
       DOUBLE PRECISION CINV(3,3),PPMOD(3,3,3,3)
C
      DO I=1,NDI
       DO J=1,NDI
         DO K=1,NDI
          DO L=1,NDI
          PPMOD(I,J,K,L)=((1.D0)/(2.D0))*(CINV(I,K)*CINV(J,L)+
     +      CINV(I,L)*CINV(J,K))-((1.D0)/(3.D0))*(CINV(I,J)*CINV(K,L))
          END DO
         END DO
       END DO
      END DO
C
      RETURN
C
      end SUBROUTINE PROJECTIONMOD
C-------------------------------------------------------------------------------
       SUBROUTINE PROJECTIONT(C,CINV,NDI,PPT)
C
       Implicit None
C
       INTEGER I,J,K,L,NDI
C
       DOUBLE PRECISION C(3,3),CINV(3,3),PPT(3,3,3,3)
C
      DO I=1,NDI
       DO J=1,NDI
         DO K=1,NDI
          DO L=1,NDI
          IF ((I.EQ.K).AND.(J.EQ.L)) THEN
          PPT(I,J,K,L)=1.D0-((1.D0)/(3.D0))*C(I,J)*CINV(K,L)
          ELSE
          PPT(I,J,K,L)=-((1.D0)/(3.D0))*C(I,J)*CINV(K,L)
          END IF
          END DO
         END DO
       END DO
      END DO
C
      RETURN
C
      end SUBROUTINE PROJECTIONT
C-------------------------------------------------------------------------------
       SUBROUTINE PULL2(A,DFGRD1,DET,PULL2A,NDI)
C
       Implicit None
C
       INTEGER NDI,I,J,K,L,ISTAT
C
       DOUBLE PRECISION A(NDI,NDI),DFGRD1(NDI,NDI),PULL2A(NDI,NDI),DET,
     1                  DFGRD1_INV(NDI,NDI),DET_DFGRD1,CONT
C
       CAll MATINV3D(DFGRD1,DFGRD1_INV,DET_DFGRD1,ISTAT)
C
      DO I=1,NDI
       DO J=1,NDI
       CONT=0.D0
         DO K=1,NDI
          DO L=1,NDI
          CONT=CONT+DET*(DFGRD1_INV(I,K)*A(K,L)*
     *      DFGRD1_INV(J,L))
          END DO
         END DO
         PULL2A(I,J)=CONT
       END DO
      END DO
C
      RETURN
C
      end SUBROUTINE
C-------------------------------------------------------------------------------
      SUBROUTINE PULL4(CC,DFGRD1,DET,PULL4CC,NDI)
C
       Implicit None
C
       INTEGER NDI,I,J,K,L,A,B,C,D,ISTAT
C
       DOUBLE PRECISION CC(NDI,NDI,NDI,NDI),DFGRD1(NDI,NDI),
     1                  PULL4CC(NDI,NDI,NDI,NDI),DET,
     1                  DFGRD1_INV(NDI,NDI),DET_DFGRD1,CONT
C
       CAll MATINV3D(DFGRD1,DFGRD1_INV,DET_DFGRD1,ISTAT)
C
      DO I=1,NDI
       DO J=1,NDI
         DO K=1,NDI
          DO L=1,NDI
          CONT=0.D0
          DO A=1,NDI
           DO B=1,NDI
            DO C=1,NDI
             DO D=1,NDI
             CONT=CONT+DET*(DFGRD1_INV(I,A)*DFGRD1_INV(J,B)*
     *       DFGRD1_INV(K,C)*DFGRD1_INV(L,D)*CC(A,B,C,D))
             END DO
            END DO
           END DO
          END DO
         PULL4CC(I,J,K,L)=CONT
         END DO
        END DO
       END DO
      END DO
C
      RETURN
C
      end SUBROUTINE
C-------------------------------------------------------------------------------
       SUBROUTINE PUSH2(A,DFGRD1,DET,PUSH2A,NDI)
C
       Implicit None
C
       INTEGER NDI,I,J,K,L
C
       DOUBLE PRECISION A(NDI,NDI),DFGRD1(NDI,NDI),PUSH2A(NDI,NDI),DET,
     1                  CONT
C
      DO I=1,NDI
       DO J=1,NDI
       CONT=0.D0
         DO K=1,NDI
          DO L=1,NDI
          CONT=CONT+(1.d0/DET)*(DFGRD1(I,K)*A(K,L)*
     *      DFGRD1(J,L))
          END DO
         END DO
         PUSH2A(I,J)=CONT
       END DO
      END DO
C
      RETURN
C
      end SUBROUTINE
C-------------------------------------------------------------------------------
      SUBROUTINE PUSH4(CC,DFGRD1,DET,PUSH4CC,NDI)
C
       Implicit None
C
       INTEGER NDI,I,J,K,L,A,B,C,D
C
       DOUBLE PRECISION CC(NDI,NDI,NDI,NDI),DFGRD1(NDI,NDI),
     1                  PUSH4CC(NDI,NDI,NDI,NDI),DET,CONT
C
      DO I=1,NDI
       DO J=1,NDI
         DO K=1,NDI
          DO L=1,NDI
          CONT=0.D0
          DO A=1,NDI
           DO B=1,NDI
            DO C=1,NDI
             DO D=1,NDI
             CONT=CONT+(1.D0/DET)*(DFGRD1(I,A)*DFGRD1(J,B)*
     *       DFGRD1(K,C)*DFGRD1(L,D)*CC(A,B,C,D))
             END DO
            END DO
           END DO
          END DO
         PUSH4CC(I,J,K,L)=CONT
         END DO
        END DO
       END DO
      END DO
C
      RETURN
C
      end SUBROUTINE
C----------------------------------------------------------------------
       SUBROUTINE TENSORPRODBARDOWN(A,B,C,NDI)
C
       Implicit None
C
       INTEGER I,J,K,L,NDI
C
       DOUBLE PRECISION A(3,3),B(3,3),C(3,3,3,3)
C
      DO I=1,NDI
       DO J=1,NDI
         DO K=1,NDI
          DO L=1,NDI
          C(I,J,K,L)=A(I,L)*B(J,K)
          END DO
         END DO
       END DO
      END DO
C
      RETURN
C
      end SUBROUTINE TENSORPRODBARDOWN
C----------------------------------------------------------------------
       SUBROUTINE TENSORPRODBARUP(A,B,C,NDI)
C
       Implicit None
C
       INTEGER I,J,K,L,NDI
C
       DOUBLE PRECISION A(3,3),B(3,3),C(3,3,3,3)
C
      DO I=1,NDI
       DO J=1,NDI
         DO K=1,NDI
          DO L=1,NDI
          C(I,J,K,L)=A(I,K)*B(J,L)
          END DO
         END DO
       END DO
      END DO
C
      RETURN
C
      end SUBROUTINE TENSORPRODBARUP
C-------------------------------------------------------------------------------
       SUBROUTINE TENSORPROD2(A,B,C,NDI)
C
       Implicit None
C
       INTEGER I,J,K,L,NDI
C
       DOUBLE PRECISION A(3,3),B(3,3),C(3,3,3,3)
C
      DO I=1,NDI
       DO J=1,NDI
         DO K=1,NDI
          DO L=1,NDI
          C(I,J,K,L)=A(I,J)*B(K,L)
          END DO
         END DO
       END DO
      END DO
C
      RETURN
C
      end SUBROUTINE TENSORPROD2
C-------------------------------------------------------------------------------
      SUBROUTINE TRACE(A,SUM)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION A(3,3),SUM
      INTEGER I,J
C
      SUM=0.D0
      DO I=1,3
       DO J=1,3
       IF (I.EQ.J) THEN
        SUM=SUM+A(I,J)
       END IF
       END DO
      END DO
C
      RETURN
C
      END SUBROUTINE TRACE
C-------------------------------------------------------------------------------
      SUBROUTINE MATINV3D(A,A_INV,DET_A,ISTAT)
C
C      RETURNS A_INV, THE INVERSE AND DET_A, THE DETERMINANT
C     DET OF THE ORIGINAL MATRIX, NOT OF THE INVERSE
C      RECURSIVE EXPRESSION OBTAINED FROM MATHEMATICA
      IMPLICIT NONE
C
      INTEGER ISTAT
C
      REAL*8 A(3,3),A_INV(3,3),DET_A,DET_A_INV
C
C
      ISTAT = 1
C
      DET_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     +        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     +        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))
C
      IF (DET_A .LE. 0.D0) THEN
        WRITE(*,*) 'WARNING: SUBROUTINE MATINV3D:'
        WRITE(*,*) 'WARNING: DET OF MAT=',DET_A
        ISTAT = 0
        RETURN
      END IF
C
      DET_A_INV = 1.D0/DET_A
C
      A_INV(1,1) = DET_A_INV*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      A_INV(1,2) = DET_A_INV*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
      A_INV(1,3) = DET_A_INV*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
      A_INV(2,1) = DET_A_INV*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
      A_INV(2,2) = DET_A_INV*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
      A_INV(2,3) = DET_A_INV*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
      A_INV(3,1) = DET_A_INV*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
      A_INV(3,2) = DET_A_INV*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
      A_INV(3,3) = DET_A_INV*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
C
C
      RETURN
      END SUBROUTINE MATINV3D
!****************************************************************************
C
      SUBROUTINE MATINV2D(A,A_INV,DET_A,ISTAT)
C
C     RETURNS A_INV, THE INVERSE, AND DET_A, THE DETERMINANT
C     NOTE THAT THE DET IS OF THE ORIGINAL MATRIX, NOT THE
C     INVERSE
C
      IMPLICIT NONE
C
      INTEGER ISTAT
C
      REAL*8 A(2,2),A_INV(2,2),DET_A,DET_A_INV
C
C
      ISTAT = 1.D0
C
      DET_A = A(1,1)*A(2,2) - A(1,2)*A(2,1)
C
      IF (DET_A .LE. 0.D0) THEN
        WRITE(*,*) 'WARNING: SUBROUTINE MATINV2D:'
        WRITE(*,*) 'WARNING: DET OF MAT=',DET_A
        ISTAT = 0
        RETURN
      END IF
C
      DET_A_INV = 1.D0/DET_A
C
      A_INV(1,1) =  DET_A_INV*A(2,2)
      A_INV(1,2) = -DET_A_INV*A(1,2)
      A_INV(2,1) = -DET_A_INV*A(2,1)
      A_INV(2,2) =  DET_A_INV*A(1,1)
C
C
      RETURN
      END SUBROUTINE MATINV2D
C
!****************************************************************************
C
      SUBROUTINE MDET(A,DET)
C
C      THIS SUBROUTINE CALCULATES THE DETERMINANT
C      OF A 3 BY 3 MATRIX [A]
C
      IMPLICIT NONE
C
      REAL*8  A(3,3),DET
C
C
      DET = A(1,1)*A(2,2)*A(3,3)
     +	  + A(1,2)*A(2,3)*A(3,1)
     +	  + A(1,3)*A(2,1)*A(3,2)
     +	  - A(3,1)*A(2,2)*A(1,3)
     +	  - A(3,2)*A(2,3)*A(1,1)
     +	  - A(3,3)*A(2,1)*A(1,2)
C
C
      RETURN
      END SUBROUTINE MDET
C
!****************************************************************************
C
      SUBROUTINE ONEM(A)
C
C      THIS SUBROUTINE STORES THE IDENTITY MATRIX IN THE
C      3 BY 3 MATRIX [A]
C
      IMPLICIT NONE
C
      INTEGER I,J
C
      DOUBLE PRECISION A(3,3)
C
C
      DO I=1,3
         DO J=1,3
	           IF (I .EQ. J) THEN
              A(I,J) = 1.D0
            ELSE
              A(I,J) = 0.D0
            END IF
         END DO
      END DO
C
C
      RETURN
      END SUBROUTINE ONEM
C
C-------------------------------------------------------------------------------
            SUBROUTINE DMG_EVOL(SE0,ECC0,DMG0,GDER0,STP,C,STRETCH4B,
     1       STRETCH6B,I1BARE,C10,K11,K12,K21,K22,ECCMAX,MAXD,GDER,
     1        DMG,ECC,SEISO)
c
            Implicit None
c
            DOUBLE PRECISION SEISO, ECC,SE0,ZD,ECC0,C,ECCMAX,GDER,
     1        DMG,AUX,AUX1,MAXD,ONE,DMG0,GDER0,AUX2,AUX3,AUX4,AUX5,
     1        C10,K11,K12,K21,K22,STRETCH4B,STRETCH6B,STE,STP,STE2,
     1        I1BARE,ZERO,SEM,SEF1,SEF2,STRETCH6B2,TWO,THREE
C
            ONE=1.0D0
            ZERO=0.0D0
            TWO=2.0D0
            THREE=3.0D0
C
C
C     STRAIN-ENERGY
C     ISOTROPIC PART
      SEM=C10*(I1BARE-THREE)
C     ANISOTROPIC PART
      STE=STRETCH4B/STP
      STE2=STE**TWO
      SEF1=(K11/(TWO*K12))*(DEXP(K12*(STE2-ONE)**TWO)-ONE)
      STRETCH6B2=STRETCH6B**TWO
      SEF2=(K21/(TWO*K22))*(DEXP(K22*(STRETCH6B2-ONE)**TWO)-ONE)
      SEISO=SEM+SEF1+SEF2
c
C      CRITERIA FOR ACCUMULATE ECC
            ZD=SEISO-SE0
c
C
           IF ((ZD.GT.ZERO).AND.(DMG0.NE.ONE)) THEN
                        IF (ECC0.GE.ECCMAX) THEN
                              DMG=MAXD*ONE
                              ECC=ABS(SEISO-SE0)+ECC0
                        ELSE
                         ECC=ABS(SEISO-SE0)+ECC0
                              IF (ECC.GE.ECCMAX) THEN
                                DMG=MAXD*ONE
                              ELSE
                                AUX=ECC/ECCMAX
                                AUX1=ONE/C
                                AUX2=ABS(ONE-AUX)
                                DMG=MAXD*(ONE-((AUX2)**AUX1))
                              END IF
                        END IF
                        AUX=ECC/ECCMAX
                        IF ((C.GT.ONE).AND.(AUX.EQ.ONE)) THEN
                              AUX=0.999
                              AUX1=ONE/C
                              AUX3=ABS(AUX-ONE)
                              AUX4=AUX-ONE
                              AUX5=AUX1-ONE
                    GDER=-(MAXD*AUX3**AUX5)*(AUX4/AUX3)*(ONE/(C*ECCMAX))
                        ELSE
                              AUX=ECC/ECCMAX
                              AUX1=ONE/C
                              AUX3=ABS(AUX-ONE)
                              AUX4=AUX-ONE
                              AUX5=AUX1-ONE
                    GDER=-(MAXD*AUX3**AUX5)*(AUX4/AUX3)*(ONE/(C*ECCMAX))
                        END IF
           ELSE
            ECC=ECC0
            GDER=GDER0
            DMG=DMG0
           END IF

           RETURN
C
           end SUBROUTINE DMG_EVOL
C
C-----------------------------------------------------------------------------
       SUBROUTINE ELASTIC_TENSORS(DFGRD1,FP4INV,DFGRD1E,NSHR,
     1      CE,CBARE,I1BARE,JAC,JACE,CINVE)
C
       IMPLICIT NONE
       INCLUDE 'PARAM_UMAT.F'
C
            INTEGER ISTATC,NSHR,ISTATCE,ISTATCBARE
c
            DOUBLE PRECISION JAC,JACE,DFGRD1(3,3),
     1  SCALE,DISTGR(3,3),CBAR(3,3),
     1  DFGRD1E(3,3),DFGRD1TE(3,3),DISTGRTE(3,3),
     2  CE(3,3),CBARE(3,3),CINVE(3,3),CBARINVE(3,3),
     3  DETCBARE,FP4INV(3,3),I1BARE,DETCE,DISTGRE(3,3)
C
C     ELASTIC DEFORMATION GRADIENT
      CALL M3MULT(DFGRD1,FP4INV,DFGRD1E)

C     ELASTIC JACOBIAN AND ISOCHORIC ELASTIC DEFORMATION TENSOR
C----------------------------------------------------------------------
C     ELASTIC JACOBIAN
      JACE=ZERO
      JACE = DFGRD1E(1,1) * DFGRD1E(2,2) * DFGRD1E(3,3)
     1    - DFGRD1E(1,2) * DFGRD1E(2,1) * DFGRD1E(3,3)
C
      IF (NSHR .EQ. 3) THEN
          JACE = JACE + DFGRD1E(1,2) * DFGRD1E(2,3) * DFGRD1E(3,1)
     1              + DFGRD1E(1,3) * DFGRD1E(3,2) * DFGRD1E(2,1)
     2              - DFGRD1E(1,3) * DFGRD1E(3,1) * DFGRD1E(2,2)
     3              - DFGRD1E(2,3) * DFGRD1E(3,2) * DFGRD1E(1,1)
      END IF
C
      SCALE=JAC**(-ONE/THREE)
      DISTGRE=SCALE*DFGRD1E
C
C     ELASTIC RIGHT CAUCHY-GREEN TENSOR AND ELASTIC INVARIANTS
C
C     TRANSPOSE DFGRD1 AND DISTGR
      DFGRD1TE=TRANSPOSE(DFGRD1E)
      DISTGRTE=TRANSPOSE(DISTGRE)
C
C     RIGHT CAUCHY-GREEN TENSORS
      CALL M3MULT(DFGRD1TE,DFGRD1E,CE)
      CALL M3MULT(DISTGRTE,DISTGRE,CBARE)
C
C     INVERSE OF ELASTIC LEFT CAUCHY-GREEN TENSOR
      CALL MATINV3D(CE,CINVE,DETCE,ISTATCE)
      CALL MATINV3D(CBARE,CBARINVE,DETCBARE,ISTATCBARE)
C
C     INVARIANTS OF C AND CBAR
      CALL TRACE(CBARE,I1BARE)
C
       RETURN
       END SUBROUTINE ELASTIC_TENSORS
C----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE MULT_VEC_MATRIX(V,M,R)
C
        INCLUDE 'ABA_PARAM.INC'
C
        DOUBLE PRECISION AUX,M(3,3),V(3),R(3)
        INTEGER I,J
C
        DO  I=1,3
        AUX=0.0D0
            DO J=1,3
                AUX=AUX+M(I,J)*V(J)
            END DO
            R(I)=AUX
        END DO
        RETURN
      END SUBROUTINE MULT_VEC_MATRIX
C-----------------------------------------------------------------------
      SUBROUTINE SUB_ARRAY(A1,A2,R)
C
        INCLUDE 'ABA_PARAM.INC'
C
        DOUBLE PRECISION A1(3),A2(3),R(3)
C
        INTEGER I
C
        DO  I=1,3
            R(I)=A1(I)-A2(I)
        END DO
        RETURN
      END SUBROUTINE SUB_ARRAY
C-----------------------------------------------------------------------