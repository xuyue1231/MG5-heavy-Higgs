      SUBROUTINE SBORN(P1,ANS)
C     
C     Generated by MadGraph5_aMC@NLO v. 2.6.6, 2018-06-28
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C     AND HELICITIES
C     FOR THE POINT IN PHASE SPACE P1(0:3,NEXTERNAL-1)
C     
C     Process: u d~ > w+ x0 > w+ e+ e- z WEIGHTED<=4 [ all = QCD ] / a
C     Process: c s~ > w+ x0 > w+ e+ e- z WEIGHTED<=4 [ all = QCD ] / a
C     Process: u d~ > w+ x0 > w+ mu+ mu- z WEIGHTED<=4 [ all = QCD ] /
C      a
C     Process: c s~ > w+ x0 > w+ mu+ mu- z WEIGHTED<=4 [ all = QCD ] /
C      a
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INCLUDE 'nexternal.inc'
      INCLUDE 'born_nhel.inc'
      INTEGER     NCOMB
      PARAMETER ( NCOMB=  144 )
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*4)
      INTEGER NGRAPHS
      PARAMETER (NGRAPHS=   1)
C     
C     ARGUMENTS 
C     
      REAL*8 P1(0:3,NEXTERNAL-1)
      COMPLEX*16 ANS(2)
C     
C     LOCAL VARIABLES 
C     
      INTEGER IHEL,IDEN,I,J,JJ,GLU_IJ
      REAL*8 BORN,BORNS(2)
      COMPLEX*16 BORNTILDE
      INTEGER NTRY(4)
      DATA NTRY /4*0/
      INTEGER NHEL(NEXTERNAL-1,NCOMB)
      DATA (NHEL(I,   1),I=1,6) / 1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,   2),I=1,6) / 1,-1,-1, 1,-1, 0/
      DATA (NHEL(I,   3),I=1,6) / 1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,   4),I=1,6) / 1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,   5),I=1,6) / 1,-1,-1, 1, 1, 0/
      DATA (NHEL(I,   6),I=1,6) / 1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,   7),I=1,6) / 1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,   8),I=1,6) / 1,-1,-1,-1,-1, 0/
      DATA (NHEL(I,   9),I=1,6) / 1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  10),I=1,6) / 1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  11),I=1,6) / 1,-1,-1,-1, 1, 0/
      DATA (NHEL(I,  12),I=1,6) / 1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  13),I=1,6) / 1,-1, 0, 1,-1,-1/
      DATA (NHEL(I,  14),I=1,6) / 1,-1, 0, 1,-1, 0/
      DATA (NHEL(I,  15),I=1,6) / 1,-1, 0, 1,-1, 1/
      DATA (NHEL(I,  16),I=1,6) / 1,-1, 0, 1, 1,-1/
      DATA (NHEL(I,  17),I=1,6) / 1,-1, 0, 1, 1, 0/
      DATA (NHEL(I,  18),I=1,6) / 1,-1, 0, 1, 1, 1/
      DATA (NHEL(I,  19),I=1,6) / 1,-1, 0,-1,-1,-1/
      DATA (NHEL(I,  20),I=1,6) / 1,-1, 0,-1,-1, 0/
      DATA (NHEL(I,  21),I=1,6) / 1,-1, 0,-1,-1, 1/
      DATA (NHEL(I,  22),I=1,6) / 1,-1, 0,-1, 1,-1/
      DATA (NHEL(I,  23),I=1,6) / 1,-1, 0,-1, 1, 0/
      DATA (NHEL(I,  24),I=1,6) / 1,-1, 0,-1, 1, 1/
      DATA (NHEL(I,  25),I=1,6) / 1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  26),I=1,6) / 1,-1, 1, 1,-1, 0/
      DATA (NHEL(I,  27),I=1,6) / 1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  28),I=1,6) / 1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  29),I=1,6) / 1,-1, 1, 1, 1, 0/
      DATA (NHEL(I,  30),I=1,6) / 1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  31),I=1,6) / 1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  32),I=1,6) / 1,-1, 1,-1,-1, 0/
      DATA (NHEL(I,  33),I=1,6) / 1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  34),I=1,6) / 1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  35),I=1,6) / 1,-1, 1,-1, 1, 0/
      DATA (NHEL(I,  36),I=1,6) / 1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  37),I=1,6) / 1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  38),I=1,6) / 1, 1,-1, 1,-1, 0/
      DATA (NHEL(I,  39),I=1,6) / 1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  40),I=1,6) / 1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  41),I=1,6) / 1, 1,-1, 1, 1, 0/
      DATA (NHEL(I,  42),I=1,6) / 1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  43),I=1,6) / 1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  44),I=1,6) / 1, 1,-1,-1,-1, 0/
      DATA (NHEL(I,  45),I=1,6) / 1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  46),I=1,6) / 1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  47),I=1,6) / 1, 1,-1,-1, 1, 0/
      DATA (NHEL(I,  48),I=1,6) / 1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  49),I=1,6) / 1, 1, 0, 1,-1,-1/
      DATA (NHEL(I,  50),I=1,6) / 1, 1, 0, 1,-1, 0/
      DATA (NHEL(I,  51),I=1,6) / 1, 1, 0, 1,-1, 1/
      DATA (NHEL(I,  52),I=1,6) / 1, 1, 0, 1, 1,-1/
      DATA (NHEL(I,  53),I=1,6) / 1, 1, 0, 1, 1, 0/
      DATA (NHEL(I,  54),I=1,6) / 1, 1, 0, 1, 1, 1/
      DATA (NHEL(I,  55),I=1,6) / 1, 1, 0,-1,-1,-1/
      DATA (NHEL(I,  56),I=1,6) / 1, 1, 0,-1,-1, 0/
      DATA (NHEL(I,  57),I=1,6) / 1, 1, 0,-1,-1, 1/
      DATA (NHEL(I,  58),I=1,6) / 1, 1, 0,-1, 1,-1/
      DATA (NHEL(I,  59),I=1,6) / 1, 1, 0,-1, 1, 0/
      DATA (NHEL(I,  60),I=1,6) / 1, 1, 0,-1, 1, 1/
      DATA (NHEL(I,  61),I=1,6) / 1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  62),I=1,6) / 1, 1, 1, 1,-1, 0/
      DATA (NHEL(I,  63),I=1,6) / 1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  64),I=1,6) / 1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  65),I=1,6) / 1, 1, 1, 1, 1, 0/
      DATA (NHEL(I,  66),I=1,6) / 1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  67),I=1,6) / 1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  68),I=1,6) / 1, 1, 1,-1,-1, 0/
      DATA (NHEL(I,  69),I=1,6) / 1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  70),I=1,6) / 1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  71),I=1,6) / 1, 1, 1,-1, 1, 0/
      DATA (NHEL(I,  72),I=1,6) / 1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  73),I=1,6) /-1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  74),I=1,6) /-1,-1,-1, 1,-1, 0/
      DATA (NHEL(I,  75),I=1,6) /-1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  76),I=1,6) /-1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  77),I=1,6) /-1,-1,-1, 1, 1, 0/
      DATA (NHEL(I,  78),I=1,6) /-1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  79),I=1,6) /-1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  80),I=1,6) /-1,-1,-1,-1,-1, 0/
      DATA (NHEL(I,  81),I=1,6) /-1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  82),I=1,6) /-1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  83),I=1,6) /-1,-1,-1,-1, 1, 0/
      DATA (NHEL(I,  84),I=1,6) /-1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  85),I=1,6) /-1,-1, 0, 1,-1,-1/
      DATA (NHEL(I,  86),I=1,6) /-1,-1, 0, 1,-1, 0/
      DATA (NHEL(I,  87),I=1,6) /-1,-1, 0, 1,-1, 1/
      DATA (NHEL(I,  88),I=1,6) /-1,-1, 0, 1, 1,-1/
      DATA (NHEL(I,  89),I=1,6) /-1,-1, 0, 1, 1, 0/
      DATA (NHEL(I,  90),I=1,6) /-1,-1, 0, 1, 1, 1/
      DATA (NHEL(I,  91),I=1,6) /-1,-1, 0,-1,-1,-1/
      DATA (NHEL(I,  92),I=1,6) /-1,-1, 0,-1,-1, 0/
      DATA (NHEL(I,  93),I=1,6) /-1,-1, 0,-1,-1, 1/
      DATA (NHEL(I,  94),I=1,6) /-1,-1, 0,-1, 1,-1/
      DATA (NHEL(I,  95),I=1,6) /-1,-1, 0,-1, 1, 0/
      DATA (NHEL(I,  96),I=1,6) /-1,-1, 0,-1, 1, 1/
      DATA (NHEL(I,  97),I=1,6) /-1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  98),I=1,6) /-1,-1, 1, 1,-1, 0/
      DATA (NHEL(I,  99),I=1,6) /-1,-1, 1, 1,-1, 1/
      DATA (NHEL(I, 100),I=1,6) /-1,-1, 1, 1, 1,-1/
      DATA (NHEL(I, 101),I=1,6) /-1,-1, 1, 1, 1, 0/
      DATA (NHEL(I, 102),I=1,6) /-1,-1, 1, 1, 1, 1/
      DATA (NHEL(I, 103),I=1,6) /-1,-1, 1,-1,-1,-1/
      DATA (NHEL(I, 104),I=1,6) /-1,-1, 1,-1,-1, 0/
      DATA (NHEL(I, 105),I=1,6) /-1,-1, 1,-1,-1, 1/
      DATA (NHEL(I, 106),I=1,6) /-1,-1, 1,-1, 1,-1/
      DATA (NHEL(I, 107),I=1,6) /-1,-1, 1,-1, 1, 0/
      DATA (NHEL(I, 108),I=1,6) /-1,-1, 1,-1, 1, 1/
      DATA (NHEL(I, 109),I=1,6) /-1, 1,-1, 1,-1,-1/
      DATA (NHEL(I, 110),I=1,6) /-1, 1,-1, 1,-1, 0/
      DATA (NHEL(I, 111),I=1,6) /-1, 1,-1, 1,-1, 1/
      DATA (NHEL(I, 112),I=1,6) /-1, 1,-1, 1, 1,-1/
      DATA (NHEL(I, 113),I=1,6) /-1, 1,-1, 1, 1, 0/
      DATA (NHEL(I, 114),I=1,6) /-1, 1,-1, 1, 1, 1/
      DATA (NHEL(I, 115),I=1,6) /-1, 1,-1,-1,-1,-1/
      DATA (NHEL(I, 116),I=1,6) /-1, 1,-1,-1,-1, 0/
      DATA (NHEL(I, 117),I=1,6) /-1, 1,-1,-1,-1, 1/
      DATA (NHEL(I, 118),I=1,6) /-1, 1,-1,-1, 1,-1/
      DATA (NHEL(I, 119),I=1,6) /-1, 1,-1,-1, 1, 0/
      DATA (NHEL(I, 120),I=1,6) /-1, 1,-1,-1, 1, 1/
      DATA (NHEL(I, 121),I=1,6) /-1, 1, 0, 1,-1,-1/
      DATA (NHEL(I, 122),I=1,6) /-1, 1, 0, 1,-1, 0/
      DATA (NHEL(I, 123),I=1,6) /-1, 1, 0, 1,-1, 1/
      DATA (NHEL(I, 124),I=1,6) /-1, 1, 0, 1, 1,-1/
      DATA (NHEL(I, 125),I=1,6) /-1, 1, 0, 1, 1, 0/
      DATA (NHEL(I, 126),I=1,6) /-1, 1, 0, 1, 1, 1/
      DATA (NHEL(I, 127),I=1,6) /-1, 1, 0,-1,-1,-1/
      DATA (NHEL(I, 128),I=1,6) /-1, 1, 0,-1,-1, 0/
      DATA (NHEL(I, 129),I=1,6) /-1, 1, 0,-1,-1, 1/
      DATA (NHEL(I, 130),I=1,6) /-1, 1, 0,-1, 1,-1/
      DATA (NHEL(I, 131),I=1,6) /-1, 1, 0,-1, 1, 0/
      DATA (NHEL(I, 132),I=1,6) /-1, 1, 0,-1, 1, 1/
      DATA (NHEL(I, 133),I=1,6) /-1, 1, 1, 1,-1,-1/
      DATA (NHEL(I, 134),I=1,6) /-1, 1, 1, 1,-1, 0/
      DATA (NHEL(I, 135),I=1,6) /-1, 1, 1, 1,-1, 1/
      DATA (NHEL(I, 136),I=1,6) /-1, 1, 1, 1, 1,-1/
      DATA (NHEL(I, 137),I=1,6) /-1, 1, 1, 1, 1, 0/
      DATA (NHEL(I, 138),I=1,6) /-1, 1, 1, 1, 1, 1/
      DATA (NHEL(I, 139),I=1,6) /-1, 1, 1,-1,-1,-1/
      DATA (NHEL(I, 140),I=1,6) /-1, 1, 1,-1,-1, 0/
      DATA (NHEL(I, 141),I=1,6) /-1, 1, 1,-1,-1, 1/
      DATA (NHEL(I, 142),I=1,6) /-1, 1, 1,-1, 1,-1/
      DATA (NHEL(I, 143),I=1,6) /-1, 1, 1,-1, 1, 0/
      DATA (NHEL(I, 144),I=1,6) /-1, 1, 1,-1, 1, 1/
      INTEGER IDEN_VALUES(4)
      DATA IDEN_VALUES /36, 36, 36, 36/
      INTEGER IJ_VALUES(4)
      DATA IJ_VALUES /1, 2, 1, 2/
C     
C     GLOBAL VARIABLES
C     
      DOUBLE PRECISION AMP2(1), JAMP2(0:1)
      COMMON/TO_AMPS/  AMP2,       JAMP2
      DATA JAMP2(0) /   1/
      LOGICAL GOODHEL(NCOMB,4)
      COMMON /C_GOODHEL/GOODHEL
      DOUBLE COMPLEX SAVEAMP(NGRAPHS,MAX_BHEL)
      COMMON/TO_SAVEAMP/SAVEAMP
      DOUBLE PRECISION SAVEMOM(NEXTERNAL-1,2)
      COMMON/TO_SAVEMOM/SAVEMOM
      DOUBLE PRECISION HEL_FAC
      INTEGER GET_HEL,SKIP(4)
      COMMON/CBORN/HEL_FAC,GET_HEL,SKIP
      LOGICAL CALCULATEDBORN
      COMMON/CCALCULATEDBORN/CALCULATEDBORN
      INTEGER NFKSPROCESS
      COMMON/C_NFKSPROCESS/NFKSPROCESS
      DOUBLE PRECISION       WGT_ME_BORN,WGT_ME_REAL
      COMMON /C_WGT_ME_TREE/ WGT_ME_BORN,WGT_ME_REAL
      LOGICAL COND_IJ
C     ----------
C     BEGIN CODE
C     ----------
      IDEN=IDEN_VALUES(NFKSPROCESS)
      GLU_IJ = IJ_VALUES(NFKSPROCESS)
      NTRY(NFKSPROCESS)=NTRY(NFKSPROCESS)+1
      IF (NTRY(NFKSPROCESS).LT.2) THEN
        IF (GLU_IJ.EQ.0) THEN
          SKIP(NFKSPROCESS)=0
        ELSE
          SKIP(NFKSPROCESS)=1
          DO WHILE(NHEL(GLU_IJ ,SKIP(NFKSPROCESS)).NE.-NHEL(GLU_IJ ,1))
            SKIP(NFKSPROCESS)=SKIP(NFKSPROCESS)+1
          ENDDO
          SKIP(NFKSPROCESS)=SKIP(NFKSPROCESS)-1
        ENDIF
      ENDIF
      DO JJ=1,NGRAPHS
        AMP2(JJ)=0D0
      ENDDO
      DO JJ=1,INT(JAMP2(0))
        JAMP2(JJ)=0D0
      ENDDO
      IF (CALCULATEDBORN) THEN
        DO J=1,NEXTERNAL-1
          IF (SAVEMOM(J,1).NE.P1(0,J) .OR. SAVEMOM(J,2).NE.P1(3,J))
     $      THEN
            CALCULATEDBORN=.FALSE.
            WRITE (*,*) 'momenta not the same in Born'
            STOP
          ENDIF
        ENDDO
      ENDIF
      IF (.NOT.CALCULATEDBORN) THEN
        DO J=1,NEXTERNAL-1
          SAVEMOM(J,1)=P1(0,J)
          SAVEMOM(J,2)=P1(3,J)
        ENDDO
        DO J=1,MAX_BHEL
          DO JJ=1,NGRAPHS
            SAVEAMP(JJ,J)=(0D0,0D0)
          ENDDO
        ENDDO
      ENDIF
      ANS(1) = 0D0
      ANS(2) = 0D0
      HEL_FAC=1D0
      DO IHEL=1,NCOMB
          ! the following lines are to avoid segfaults when glu_ij=0
        COND_IJ=SKIP(NFKSPROCESS).EQ.0
        IF (.NOT.COND_IJ) COND_IJ=COND_IJ.OR.NHEL(GLU_IJ,IHEL)
     $   .EQ.NHEL(GLU_IJ,1)
          !if (nhel(glu_ij,ihel).EQ.NHEL(GLU_IJ,1).or.skip(nfksprocess).eq.0) then
        IF (COND_IJ) THEN
          IF ((GOODHEL(IHEL,NFKSPROCESS) .OR. GOODHEL(IHEL
     $     +SKIP(NFKSPROCESS),NFKSPROCESS) .OR. NTRY(NFKSPROCESS) .LT.
     $      2) ) THEN
            ANS(1)=ANS(1)+BORN(P1,NHEL(1,IHEL),IHEL,BORNTILDE,BORNS)
            ANS(2)=ANS(2)+BORNTILDE
            IF ( BORNS(1).NE.0D0 .AND. .NOT. GOODHEL(IHEL,NFKSPROCESS)
     $        ) THEN
              GOODHEL(IHEL,NFKSPROCESS)=.TRUE.
            ENDIF
            IF ( BORNS(2).NE.0D0 .AND. .NOT. GOODHEL(IHEL
     $       +SKIP(NFKSPROCESS),NFKSPROCESS) ) THEN
              GOODHEL(IHEL+SKIP(NFKSPROCESS),NFKSPROCESS)=.TRUE.
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      ANS(1)=ANS(1)/DBLE(IDEN)
      ANS(2)=ANS(2)/DBLE(IDEN)
      WGT_ME_BORN=DBLE(ANS(1))
      CALCULATEDBORN=.TRUE.
      END


      REAL*8 FUNCTION BORN(P,NHEL,HELL,BORNTILDE,BORNS)
C     
C     Generated by MadGraph5_aMC@NLO v. 2.6.6, 2018-06-28
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C     FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL-1)

C     Process: u d~ > w+ x0 > w+ e+ e- z WEIGHTED<=4 [ all = QCD ] / a
C     Process: c s~ > w+ x0 > w+ e+ e- z WEIGHTED<=4 [ all = QCD ] / a
C     Process: u d~ > w+ x0 > w+ mu+ mu- z WEIGHTED<=4 [ all = QCD ] /
C      a
C     Process: c s~ > w+ x0 > w+ mu+ mu- z WEIGHTED<=4 [ all = QCD ] /
C      a
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS,    NEIGEN
      PARAMETER (NGRAPHS=   1,NEIGEN=  1)
      INTEGER    NWAVEFUNCS, NCOLOR
      PARAMETER (NWAVEFUNCS=9, NCOLOR=1)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1
      PARAMETER (IMAG1 = (0D0,1D0))
      INCLUDE 'nexternal.inc'
      INCLUDE 'born_nhel.inc'
      INCLUDE 'coupl.inc'
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL-1),BORNS(2)
      INTEGER NHEL(NEXTERNAL-1), HELL
      COMPLEX*16 BORNTILDE
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J,IHEL,BACK_HEL,GLU_IJ
      INTEGER IC(NEXTERNAL-1),NMO
      PARAMETER (NMO=NEXTERNAL-1)
      DATA IC /NMO*1/
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 ZTEMP, AMP(NGRAPHS), JAMP(NCOLOR), W(8,NWAVEFUNCS),
     $  JAMPH(2, NCOLOR)
C     
C     GLOBAL VARIABLES
C     
      DOUBLE PRECISION AMP2(NGRAPHS), JAMP2(0:NCOLOR)
      COMMON/TO_AMPS/  AMP2,       JAMP2
      DOUBLE COMPLEX SAVEAMP(NGRAPHS,MAX_BHEL)
      COMMON/TO_SAVEAMP/SAVEAMP
      DOUBLE PRECISION HEL_FAC
      INTEGER GET_HEL,SKIP(4)
      COMMON/CBORN/HEL_FAC,GET_HEL,SKIP
      LOGICAL CALCULATEDBORN
      COMMON/CCALCULATEDBORN/CALCULATEDBORN
      INTEGER NFKSPROCESS
      COMMON/C_NFKSPROCESS/NFKSPROCESS
      INTEGER STEP_HEL
      LOGICAL COND_IJ
      INTEGER IJ_VALUES(4)
      DATA IJ_VALUES /1, 2, 1, 2/
C     
C     COLOR DATA
C     
      DATA DENOM(1)/1/
      DATA (CF(I,  1),I=  1,  1) /    3/
C     1 T(2,1)
C     ----------
C     BEGIN CODE
C     ----------
      GLU_IJ = IJ_VALUES(NFKSPROCESS)
      BORN = 0D0
      BORNTILDE = (0D0,0D0)
      BACK_HEL = NHEL(GLU_IJ)
      BORNS(1) = 0D0
      BORNS(2) = 0D0
      IF (GLU_IJ.NE.0) THEN
        BACK_HEL = NHEL(GLU_IJ)
        IF (BACK_HEL.NE.0) THEN
          STEP_HEL=-2*BACK_HEL
        ELSE
          STEP_HEL=1
        ENDIF
      ELSE
        BACK_HEL=0
        STEP_HEL=1
      ENDIF
      DO IHEL=BACK_HEL,-BACK_HEL,STEP_HEL
        IF (GLU_IJ.NE.0) THEN
          COND_IJ=IHEL.EQ.BACK_HEL.OR.NHEL(GLU_IJ).NE.0
        ELSE
          COND_IJ=IHEL.EQ.BACK_HEL
        ENDIF
        IF (COND_IJ) THEN
          IF (GLU_IJ.NE.0) THEN
            IF (NHEL(GLU_IJ).NE.0) NHEL(GLU_IJ) = IHEL
          ENDIF
          IF (.NOT. CALCULATEDBORN) THEN
            CALL IXXXXX(P(0,1),ZERO,NHEL(1),+1*IC(1),W(1,1))
            CALL OXXXXX(P(0,2),ZERO,NHEL(2),-1*IC(2),W(1,2))
            CALL VXXXXX(P(0,3),MDL_MW,NHEL(3),+1*IC(3),W(1,3))
            CALL IXXXXX(P(0,4),ZERO,NHEL(4),-1*IC(4),W(1,4))
            CALL OXXXXX(P(0,5),ZERO,NHEL(5),+1*IC(5),W(1,5))
            CALL VXXXXX(P(0,6),MDL_MZ,NHEL(6),+1*IC(6),W(1,6))
            CALL FFV2_3(W(1,1),W(1,2),GC_47,MDL_MW,MDL_WW,W(1,7))
            CALL FFV2_4_3(W(1,4),W(1,5),-GC_22,GC_24,MDL_MZ,MDL_WZ,W(1
     $       ,2))
            CALL VVS10_1_11_2_12_3(W(1,7),W(1,3),GC_3005A,GC_3005H1
     $       ,GC_3005H2,GC_3005H3,GC_3005H4,MDL_MX0,MDL_WX0,W(1,5))
C           Amplitude(s) for diagram number 1
            CALL VVS10_1_2_13_0(W(1,2),W(1,6),W(1,5),GC_3007A
     $       ,GC_3007H1,GC_3007H2,GC_3007H3,AMP(1))
            DO I=1,NGRAPHS
              IF(IHEL.EQ.BACK_HEL)THEN
                SAVEAMP(I,HELL)=AMP(I)
              ELSEIF(IHEL.EQ.-BACK_HEL)THEN
                SAVEAMP(I,HELL+SKIP(NFKSPROCESS))=AMP(I)
              ELSE
                WRITE(*,*) 'ERROR #1 in born.f'
                STOP
              ENDIF
            ENDDO
          ELSEIF (CALCULATEDBORN) THEN
            DO I=1,NGRAPHS
              IF(IHEL.EQ.BACK_HEL)THEN
                AMP(I)=SAVEAMP(I,HELL)
              ELSEIF(IHEL.EQ.-BACK_HEL)THEN
                AMP(I)=SAVEAMP(I,HELL+SKIP(NFKSPROCESS))
              ELSE
                WRITE(*,*) 'ERROR #1 in born.f'
                STOP
              ENDIF
            ENDDO
          ENDIF
          JAMP(1)=+AMP(1)
          DO I = 1, NCOLOR
            ZTEMP = (0.D0,0.D0)
            DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
            ENDDO
            BORNS(2-(1+BACK_HEL*IHEL)/2)=BORNS(2-(1+BACK_HEL*IHEL)/2)
     $       +ZTEMP*DCONJG(JAMP(I))/DENOM(I)
          ENDDO
          DO I = 1, NGRAPHS
            AMP2(I)=AMP2(I)+AMP(I)*DCONJG(AMP(I))
          ENDDO
          DO I = 1, NCOLOR
            JAMP2(I)=JAMP2(I)+JAMP(I)*DCONJG(JAMP(I))
            JAMPH(2-(1+BACK_HEL*IHEL)/2,I)=JAMP(I)
          ENDDO
        ENDIF
      ENDDO
      BORN=BORNS(1)+BORNS(2)
      DO I = 1, NCOLOR
        ZTEMP = (0.D0,0.D0)
        DO J = 1, NCOLOR
          ZTEMP = ZTEMP + CF(J,I)*JAMPH(2,J)
        ENDDO
        BORNTILDE = BORNTILDE + ZTEMP*DCONJG(JAMPH(1,I))/DENOM(I)
      ENDDO
      IF (GLU_IJ.NE.0) NHEL(GLU_IJ) = BACK_HEL
      END


      BLOCK DATA GOODHELS
      INTEGER     NCOMB
      PARAMETER ( NCOMB=  144 )
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*4)
      LOGICAL GOODHEL(NCOMB,4)
      COMMON /C_GOODHEL/GOODHEL
      DATA GOODHEL/THEL*.FALSE./
      END
