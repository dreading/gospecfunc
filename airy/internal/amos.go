// Copyright 2019 Infin IT Pty Ltd. All rights reserved.
//
// Use of this source code is governed by a BSe-style
// license that can be found in the LICENSE file.
//
// The original Fortran code are from
//  ALGORITHM 644, TRANSACTIONS ON MATHEMATICAL SOFTWARE,
//	VOL. 21, NO. 4, December, 1995, P.  388--393.
//  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES

package amos

import (
	"math"
	"math/cmplx"
)

var I1MACH = []int{-0, 5, 6, 0, 0, 32, 4, 2, 31, 2147483647, 2, 24, -125, 127, 53, -1021, 1023}
var D1MACH = []float64{math.NaN(), 2.23e-308, 1.79e-308, 1.11e-16, 2.22e-16, 0.30103000998497009}

// DGAMLN COMPUTES THE NATURAL LOG OF THE GAMMA FUNCTION
func DGAMLN(z float64) float64 {
	r, _ := math.Lgamma(z)
	return r
}

// ZABS COMPUTES THE ABSOLUTE VALUE OR MAGNITUDE OF A COMPLEX NUMBER REAL AND IMAG PARTS
func ZABS(ZR float64, ZI float64) float64 {
	return cmplx.Abs(complex(ZR, ZI))
}

func ZFLOATS(Z complex128) (float64, float64) {
	return real(Z), imag(Z)
}

// COMPLEX LOGARITHM
func ZLOG(AR float64, AI float64) (float64, float64) {
	return ZFLOATS(cmplx.Log(complex(AR, AI)))
}

// COMPLEX MULTIPLY
func ZMLT(AR float64, AI float64, BR float64, BI float64) (float64, float64) {
	return ZFLOATS(complex(AR, AI) * complex(BR, BI))
}

// COMPLEX SQUARE ROOT
func ZSQRT(AR float64, AI float64) (float64, float64) {
	return ZFLOATS(cmplx.Sqrt(complex(AR, AI)))
}

// COMPLEX DIVIDE
func ZDIV(AR float64, AI float64, BR float64, BI float64) (float64, float64) {
	return ZFLOATS(complex(AR, AI) / complex(BR, BI))
}

// COMPLEX EXPONENTIAL FUNCTION
func ZEXP(AR float64, AI float64) (float64, float64) {
	return ZFLOATS(cmplx.Exp(complex(AR, AI)))
}

// ABS INT
func ABS(a int) int {
	if a >= 0 {
		return a
	}
	return -a
}

// MIN INT
func MIN(a int, b int) int {
	if a < b {
		return a
	}
	return b
}

// MAX INT
func MAX(a int, b int) int {
	if a > b {
		return a
	}
	return b
}

//ZSHCH COMPUTES THE COMPLEX HYPERBOLIC FUNCTIONS SINH(X+I*Y) AND COSH(X+I*Y)
//RETURNS REAL AND IMAG PARTS OF SINH AND COSH
func ZSHCH(ZR float64, ZI float64) (float64, float64, float64, float64) {
	SH := math.Sinh(ZR)
	CH := math.Cosh(ZR)
	SN, CN := math.Sincos(ZI)
	return SH * CN, CH * SN, CH * CN, SH * SN
}

// ZS1S2 TESTS FOR A POSSIBLE UNDERFLOW RESULTING FROM THE ADDITION OF THE I AND K FUNCTIONS IN THE ANALYTIC CON-
// TINUATION FORMULA WHERE S1=K FUNCTION AND S2=I FUNCTION.
// ON KODE=1 THE I AND K FUNCTIONS ARE DIFFERENT ORDERS OF MAGNITUDE, BUT FOR KODE=2 THEY CAN BE OF THE SAME ORDER
// OF MAGNITUDE AND THE MAXIMUM MUST BE AT LEAST ONE PRECISION ABOVE THE UNDERFLOW LIMIT.
func ZS1S2(ZRR float64, ZRI float64, S1R float64, S1I float64, S2R float64, S2I float64, NZ int, ASCLE float64,
	ALIM float64, IUF int) (float64, float64, float64, float64, float64, float64, int, float64, float64, int) {

	const (
		ZEROR = 0.0e0
		ZEROI = 0.0e0
	)

	var AA, ALN, AS1, AS2, C1I, C1R, S1DI, S1DR float64

	NZ = 0
	AS1 = ZABS(S1R, S1I)
	AS2 = ZABS(S2R, S2I)
	if S1R == 0.0e0 && S1I == 0.0e0 {
		goto L10
	}
	if AS1 == 0.0e0 {
		goto L10
	}
	ALN = -ZRR - ZRR + math.Log(AS1)
	S1DR = S1R
	S1DI = S1I
	S1R = ZEROR
	S1I = ZEROI
	AS1 = ZEROR
	if ALN < -ALIM {
		goto L10
	}
	C1R, C1I = ZLOG(S1DR, S1DI)
	C1R = C1R - ZRR - ZRR
	C1I = C1I - ZRI - ZRI
	S1R, S1I = ZEXP(C1R, C1I)
	AS1 = ZABS(S1R, S1I)
	IUF = IUF + 1
L10:
	AA = math.Max(AS1, AS2)
	if AA > ASCLE {
		return ZRR, ZRI, S1R, S1I, S2R, S2I, NZ, ASCLE, ALIM, IUF
	}
	S1R = ZEROR
	S1I = ZEROI
	S2R = ZEROR
	S2I = ZEROI
	NZ = 1
	IUF = 0
	return ZRR, ZRI, S1R, S1I, S2R, S2I, NZ, ASCLE, ALIM, IUF
}

// ZUCHK TAKES A SCALED QUANTITY Y WHOSE MAGNITUDE IS GREATER THAN EXP(-ALIM)=ASCLE=1.0E+3*D1MACH(1)/TOL.
// THE TEST IS MADE TO SEE if THE MAGNITUDE OF THE REAL OR IMAGINARY PART WOULD UNDERFLOW
// WHEN Y IS SCALED (BY TOL) TO ITS PROPER VALUE. Y IS ACCEPTED if THE UNDERFLOW IS AT LEAST ONE PRECISION
// BELOW THE MAGNITUDE OF THE LARGEST COMPONENT; OTHERWISE THE PHASE ANGLE DOES NOT HAVE
// ABSOLUTE ACCURACY AND AN UNDERFLOW IS ASSUMED.
func ZUCHK(YR float64, YI float64, NZ int, ASCLE float64, TOL float64) (float64, float64, int, float64, float64) {

	var SS, ST, WR, WI float64

	NZ = 0
	WR = math.Abs(YR)
	WI = math.Abs(YI)
	ST = math.Min(WR, WI)
	if ST > ASCLE {
		return YR, YI, NZ, ASCLE, TOL
	}
	SS = math.Max(WR, WI)
	ST = ST / TOL
	if SS < ST {
		NZ = 1
	}
	return YR, YI, NZ, ASCLE, TOL
}

// ZACAI APPLIES THE ANALYTIC CONTINUATION FORMULA
//
//    K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
//            MP=PI*MR*CMPLX(0.0,1.0)
//
// TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
// HALF Z PLANE FOR USE WITH ZAIRY WHERE FNU=1/3 OR 2/3 AND N=1.
// ZACAI IS THE SAME AS ZACON WITH THE PARTS FOR LARGER ORDERS AND
// RECURRENCE REMOVED. A RECURSIVE CALL TO ZACON CAN RESULT if ZACON
// IS CALLED FROM ZAIRY.
func ZACAI(ZR float64, ZI float64, FNU float64, KODE int, MR int, N int,
	YR []float64, YI []float64, NZ int, RL float64, TOL float64, ELIM float64, ALIM float64) (float64, float64, float64, int, int, int,
	[]float64, []float64, int, float64, float64, float64, float64) {

	// COMPLEX CSGN,CSPN,C1,C2,Y,Z,ZN,CY
	var ARG, ASCLE, AZ, CSGNR, CSGNI, CSPNR, CSPNI, C1R, C1I, C2R, C2I, DFNU, FMR, SGN, YY, ZNR, ZNI float64
	var INU, IUF, NN, NW int

	CYR := []float64{math.NaN(), 0, 0}
	CYI := []float64{math.NaN(), 0, 0}

	NZ = 0
	ZNR = -ZR
	ZNI = -ZI
	AZ = ZABS(ZR, ZI)
	NN = N
	DFNU = FNU + float64(float32(N-1))
	if AZ <= 2.0e0 {
		goto L10
	}
	if AZ*AZ*0.25e0 > DFNU+1.0e0 {
		goto L20
	}
L10:
	// POWER SERIES FOR THE I FUNCTION
	ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, TOL, ELIM, ALIM = ZSERI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, TOL, ELIM, ALIM)
	goto L40
L20:
	if AZ < RL {
		goto L30
	}
	// ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION
	ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, RL, TOL, ELIM, ALIM = ZASYI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, RL, TOL, ELIM, ALIM)
	if NW < 0 {
		goto L80
	}
	goto L40
L30:
	// MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION
	ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, TOL = ZMLRI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, TOL)
	if NW < 0 {
		goto L80
	}
L40:
	// ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
	ZNR, ZNI, FNU, KODE, _, CYR, CYI, NW, TOL, ELIM, ALIM = ZBKNU(ZNR, ZNI, FNU, KODE, 1, CYR, CYI, NW, TOL, ELIM, ALIM)
	if NW != 0 {
		goto L80
	}
	FMR = float64(float32(MR))
	SGN = -math.Copysign(math.Pi, FMR)
	CSGNR = 0.0e0
	CSGNI = SGN
	if KODE == 1 {
		goto L50
	}
	YY = -ZNI
	CSGNR = -CSGNI * math.Sin(YY)
	CSGNI = CSGNI * math.Cos(YY)
L50:
	// CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE WHEN FNU IS LARGE
	INU = int(float32(FNU))
	ARG = (FNU - float64(float32(INU))) * SGN
	CSPNR = math.Cos(ARG)
	CSPNI = math.Sin(ARG)
	if INU%2 == 0 {
		goto L60
	}
	CSPNR = -CSPNR
	CSPNI = -CSPNI
L60:
	C1R = CYR[1]
	C1I = CYI[1]
	C2R = YR[1]
	C2I = YI[1]
	if KODE == 1 {
		goto L70
	}
	IUF = 0
	ASCLE = 1000.0e0 * D1MACH[1] / TOL
	ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF = ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF)
	NZ = NZ + NW
L70:
	YR[1] = CSPNR*C1R - CSPNI*C1I + CSGNR*C2R - CSGNI*C2I
	YI[1] = CSPNR*C1I + CSPNI*C1R + CSGNR*C2I + CSGNI*C2R
	return ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, RL, TOL, ELIM, ALIM
L80:
	NZ = -1
	if NW == -2 {
		NZ = -2
	}
	return ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, RL, TOL, ELIM, ALIM
}

// ZAIRY COMPUTES THE AIRY FUNCTIONS AI(Z) AND DAI(Z) FOR COMPLEX Z
func ZAIRY(ZR float64, ZI float64, ID int, KODE int) (float64, float64, int, int) {

	//    ON KODE=1, ZAIRY COMPUTES THE COMPLEX AIRY FUNCTION AI(Z) OR
	//    ITS DERIVATIVE DAI(Z)/DZ ON ID=0 OR ID=1 RESPECTIVELY. ON
	//    KODE=2, A SCALING OPTION CEXP(ZTA)*AI(Z) OR CEXP(ZTA)*
	//    DAI(Z)/DZ IS PROVIDED TO REMOVE THE EXPONENTIAL DECAY IN
	//    -PI/3 < ARG(Z) < PI/3 AND THE EXPONENTIAL GROWTH IN
	//    PI/3 < ABS(ARG(Z)) < PI WHERE ZTA=(2/3)*Z*CSQRT(Z).
	//
	//    WHILE THE AIRY FUNCTIONS AI(Z) AND DAI(Z)/DZ ARE ANALYTIC IN
	//    THE WHOLE Z PLANE, THE CORRESPONDING SCALED FUNCTIONS DEFINED
	//    FOR KODE=2 HAVE A CUT ALONG THE NEGATIVE REAL AXIS.
	//    DEFINTIONS AND NOTATION ARE FOUND IN THE NBS HANDBOOK OF
	//    MATHEMATICAL FUNCTIONS (REF. 1).
	//
	//    INPUT      ZR,ZI ARE DOUBLE PRECISION
	//      ZR,ZI  - Z=CMPLX(ZR,ZI)
	//      ID     - ORDER OF DERIVATIVE, ID=0 OR ID=1
	//      KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
	//               KODE= 1  RETURNS
	//                        AI=AI(Z)                ON ID=0 OR
	//                        AI=DAI(Z)/DZ            ON ID=1
	//                   = 2  RETURNS
	//                        AI=CEXP(ZTA)*AI(Z)       ON ID=0 OR
	//                        AI=CEXP(ZTA)*DAI(Z)/DZ   ON ID=1 WHERE
	//                        ZTA=(2/3)*Z*CSQRT(Z)
	//
	//    OUTPUT     AIR,AII ARE DOUBLE PRECISION
	//      AIR,AII- COMPLEX ANSWER DEPENDING ON THE CHOICES FOR ID AND
	//               KODE
	//      NZ     - UNDERFLOW INDICATOR
	//               NZ= 0   , NORMAL RETURN
	//               NZ= 1   , AI=CMPLX(0.0e0,0.0e0) DUE TO UNDERFLOW IN
	//                         -PI/3 < ARG(Z) < PI/3 ON KODE=1
	//      IERR   - ERROR FLAG
	//               IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
	//               IERR=1, INPUT ERROR   - NO COMPUTATION
	//               IERR=2, OVERFLOW      - NO COMPUTATION, REAL(ZTA)
	//                       TOO LARGE ON KODE=1
	//               IERR=3, CABS(Z) LARGE      - COMPUTATION COMPLETED
	//                       LOSSES OF SIGNIFCANCE BY ARGUMENT REDUCTION
	//                       PRODUCE LESS THAN HALF OF MACHINE ACCURACY
	//               IERR=4, CABS(Z) TOO LARGE  - NO COMPUTATION
	//                       COMPLETE LOSS OF ACCURACY BY ARGUMENT
	//                       REDUCTION
	//               IERR=5, ERROR              - NO COMPUTATION,
	//                       ALGORITHM TERMINATION CONDITION NOT MET
	//
	//    AI AND DAI ARE COMPUTED FOR CABS(Z) > 1.0 FROM THE K BESSEL
	//    FUNCTIONS BY
	//
	//       AI(Z)=C*SQRT(Z)*K(1/3,ZTA) , DAI(Z)=-C*Z*K(2/3,ZTA)
	//                      C=1.0/(PI*SQRT(3.0))
	//                       ZTA=(2/3)*Z**(3/2)
	//
	//    WITH THE POWER SERIES FOR CABS(Z) <= 1.0.
	//
	//    IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
	//    MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z IS LARGE, LOSSES
	//    OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR. CONSEQUENTLY, IF
	//    THE MAGNITUDE OF ZETA=(2/3)*Z**1.5 EXCEEDS U1=SQRT(0.5/UR),
	//    THEN LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR
	//    FLAG IERR=3 IS TRIGGERED WHERE UR=math.Max(D1MACH(4),1.0D-18) IS
	//    DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
	//    ALSO, if THE MAGNITUDE OF ZETA IS LARGER THAN U2=0.5/UR, THEN
	//    ALL SIGNIFICANCE IS LOST AND IERR=4. IN ORDER TO USE THE INT
	//    FUNCTION, ZETA MUST BE FURTHER RESTRICTED NOT TO EXCEED THE
	//    LARGEST INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF ZETA
	//    MUST BE RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2,
	//    AND U3 ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE
	//    PRECISION ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE
	//    PRECISION ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMIT-
	//    ING IN THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT THE MAG-
	//    NITUDE OF Z CANNOT EXCEED 3.1E+4 IN SINGLE AND 2.1E+6 IN
	//    DOUBLE PRECISION ARITHMETIC. THIS ALSO MEANS THAT ONE CAN
	//    EXPECT TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES,
	//    NO DIGITS IN SINGLE PRECISION AND ONLY 7 DIGITS IN DOUBLE
	//    PRECISION ARITHMETIC. SIMILAR CONSIDERATIONS HOLD FOR OTHER
	//    MACHINES.
	//
	//    THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
	//    BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
	//    ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
	//    SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
	//    ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
	//    ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
	//    CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
	//    HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
	//    ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
	//    SEVERAL ORDERS OF MAGNITUDE. if ONE COMPONENT IS 10**K LARGER
	//    THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
	//    0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
	//    THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
	//    COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
	//    BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
	//    COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
	//    MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
	//    THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
	//    OR -PI/2+P.
	//
	// REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
	//            AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
	//            COMMERCE, 1955.
	//
	//          COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
	//            AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
	//
	//          A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
	//            ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
	//            1018, MAY, 1985
	//
	//          A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
	//            ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, TRANS.
	//            MATH. SOFTWARE, 1986

	const (
		TTH   = 6.66666666666666667e-01
		C1    = 3.55028053887817240e-01
		C2    = 2.58819403792806799e-01
		COEF  = 1.83776298473930683e-01
		ZEROR = 0.0e0
		ZEROI = 0.0e0
		CONER = 1.0e0
		CONEI = 0.0e0
	)

	var AIR, AII, AA, AD, AK, ALIM, ATRM, AZ, AZ3, BK, CC, CK, CSQI, CSQR, DIG, DK, D1, D2, ELIM,
		FID, FNU, PTR, RL, R1M5, SFAC, STI, STR, S1I, S1R, S2I, S2R, TOL, TRM1I, TRM1R, TRM2I,
		TRM2R, ZTAI, ZTAR, Z3I, Z3R, ALAZ, BB float64
	var IFLAG, K, K1, K2, MR, NN, NZ, IERR int

	var CYR = []float64{math.NaN(), 0}
	var CYI = []float64{math.NaN(), 0}

	IERR = 0
	NZ = 0
	if ID < 0 || ID > 1 {
		IERR = 1
	}
	if KODE < 1 || KODE > 2 {
		IERR = 1
	}
	if IERR != 0 {
		return AIR, AII, NZ, IERR
	}
	AZ = ZABS(ZR, ZI)
	TOL = math.Max(D1MACH[4], 1.0e-18)
	FID = float64(float32(ID))
	if AZ > 1.0e0 {
		goto L70
	}

	// POWER SERIES FOR CABS(Z) <= 1.
	S1R = CONER
	S1I = CONEI
	S2R = CONER
	S2I = CONEI
	if AZ < TOL {
		goto L170
	}
	AA = AZ * AZ
	if AA < TOL/AZ {
		goto L40
	}
	TRM1R = CONER
	TRM1I = CONEI
	TRM2R = CONER
	TRM2I = CONEI
	ATRM = 1.0e0
	STR = ZR*ZR - ZI*ZI
	STI = ZR*ZI + ZI*ZR
	Z3R = STR*ZR - STI*ZI
	Z3I = STR*ZI + STI*ZR
	AZ3 = AZ * AA
	AK = 2.0e0 + FID
	BK = 3.0e0 - FID - FID
	CK = 4.0e0 - FID
	DK = 3.0e0 + FID + FID
	D1 = AK * DK
	D2 = BK * CK
	AD = math.Min(D1, D2)
	AK = 24.0e0 + 9.0e0*FID
	BK = 30.0e0 - 9.0e0*FID
	for K = 1; K <= 25; K++ {
		STR = (TRM1R*Z3R - TRM1I*Z3I) / D1
		TRM1I = (TRM1R*Z3I + TRM1I*Z3R) / D1
		TRM1R = STR
		S1R = S1R + TRM1R
		S1I = S1I + TRM1I
		STR = (TRM2R*Z3R - TRM2I*Z3I) / D2
		TRM2I = (TRM2R*Z3I + TRM2I*Z3R) / D2
		TRM2R = STR
		S2R = S2R + TRM2R
		S2I = S2I + TRM2I
		ATRM = ATRM * AZ3 / AD
		D1 = D1 + AK
		D2 = D2 + BK
		AD = math.Min(D1, D2)
		if ATRM < TOL*AD {
			goto L40
		}
		AK = AK + 18.0e0
		BK = BK + 18.0e0
	}
L40:
	if ID == 1 {
		goto L50
	}
	AIR = S1R*C1 - C2*(ZR*S2R-ZI*S2I)
	AII = S1I*C1 - C2*(ZR*S2I+ZI*S2R)
	if KODE == 1 {
		return AIR, AII, NZ, IERR
	}
	STR, STI = ZSQRT(ZR, ZI)
	ZTAR = TTH * (ZR*STR - ZI*STI)
	ZTAI = TTH * (ZR*STI + ZI*STR)
	STR, STI = ZEXP(ZTAR, ZTAI)
	PTR = AIR*STR - AII*STI
	AII = AIR*STI + AII*STR
	AIR = PTR
	return AIR, AII, NZ, IERR
L50:
	AIR = -S2R * C2
	AII = -S2I * C2
	if AZ <= TOL {
		goto L60
	}
	STR = ZR*S1R - ZI*S1I
	STI = ZR*S1I + ZI*S1R
	CC = C1 / (1.0e0 + FID)
	AIR = AIR + CC*(STR*ZR-STI*ZI)
	AII = AII + CC*(STR*ZI+STI*ZR)
L60:
	if KODE == 1 {
		return AIR, AII, NZ, IERR
	}
	STR, STI = ZSQRT(ZR, ZI)
	ZTAR = TTH * (ZR*STR - ZI*STI)
	ZTAI = TTH * (ZR*STI + ZI*STR)
	STR, STI = ZEXP(ZTAR, ZTAI)
	PTR = STR*AIR - STI*AII
	AII = STR*AII + STI*AIR
	AIR = PTR
	return AIR, AII, NZ, IERR

	// CASE FOR CABS(Z) > 1.0

L70:
	FNU = (1.0e0 + FID) / 3.0e0

	// SET PARAMETERS RELATED TO MACHINE CONSTANTS.
	// TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0D-18.
	// ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
	// EXP(-ELIM) < EXP(-ALIM)=EXP(-ELIM)/TOL    AND
	// EXP(ELIM) > EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
	// UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
	// RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
	// DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).

	K1 = I1MACH[15]
	K2 = I1MACH[16]
	R1M5 = D1MACH[5]
	K = MIN(ABS(K1), ABS(K2))
	ELIM = 2.303e0 * (float64(float32(K))*R1M5 - 3.0e0)
	K1 = I1MACH[14] - 1
	AA = R1M5 * float64(float32(K1))
	DIG = math.Min(AA, 18.0e0)
	AA = AA * 2.303e0
	ALIM = ELIM + math.Max(-AA, -41.45e0)
	RL = 1.2e0*DIG + 3.0e0
	ALAZ = math.Log(AZ)

	// TEST FOR PROPER RANGE
	AA = 0.5e0 / TOL
	BB = float64(float32(I1MACH[9])) * 0.5e0
	AA = math.Min(AA, BB)
	AA = math.Pow(AA, TTH)
	if AZ > AA {
		goto L260
	}
	AA = math.Sqrt(AA)
	if AZ > AA {
		IERR = 3
	}
	CSQR, CSQI = ZSQRT(ZR, ZI)
	ZTAR = TTH * (ZR*CSQR - ZI*CSQI)
	ZTAI = TTH * (ZR*CSQI + ZI*CSQR)

	// RE(ZTA) <= 0 WHEN RE(Z) < 0, ESPECIALLY WHEN IM(Z) IS SMALL
	IFLAG = 0
	SFAC = 1.0e0
	AK = ZTAI
	if ZR >= 0.0e0 {
		goto L80
	}
	BK = ZTAR
	CK = -math.Abs(BK)
	ZTAR = CK
	ZTAI = AK
L80:
	if ZI != 0.0e0 {
		goto L90
	}
	if ZR > 0.0e0 {
		goto L90
	}
	ZTAR = 0.0e0
	ZTAI = AK
L90:
	AA = ZTAR
	if AA >= 0.0e0 && ZR > 0.0e0 {
		goto L110
	}
	if KODE == 2 {
		goto L100
	}

	// OVERFLOW TEST
	if AA > -ALIM {
		goto L100
	}
	AA = -AA + 0.25e0*ALAZ
	IFLAG = 1
	SFAC = TOL
	if AA > ELIM {
		goto L270
	}
L100:
	// CBKNU AND CACON RETURN EXP(ZTA)*K(FNU,ZTA) ON KODE=2
	MR = 1
	if ZI < 0.0e0 {
		MR = -1
	}
	ZTAR, ZTAI, FNU, KODE, MR, _, CYR, CYI, NN, RL, TOL, ELIM, ALIM = ZACAI(ZTAR, ZTAI, FNU, KODE, MR, 1, CYR, CYI, NN, RL, TOL, ELIM, ALIM)
	if NN < 0 {
		goto L280
	}
	NZ = NZ + NN
	goto L130
L110:
	if KODE == 2 {
		goto L120
	}

	// UNDERFLOW TEST
	if AA < ALIM {
		goto L120
	}
	AA = -AA - 0.25e0*ALAZ
	IFLAG = 2
	SFAC = 1.0e0 / TOL
	if AA < -ELIM {
		goto L210
	}
L120:
	ZTAR, ZTAI, FNU, KODE, _, CYR, CYI, NZ, TOL, ELIM, ALIM = ZBKNU(ZTAR, ZTAI, FNU, KODE, 1, CYR, CYI, NZ, TOL, ELIM, ALIM)
L130:
	S1R = CYR[1] * COEF
	S1I = CYI[1] * COEF
	if IFLAG != 0 {
		goto L150
	}
	if ID == 1 {
		goto L140
	}
	AIR = CSQR*S1R - CSQI*S1I
	AII = CSQR*S1I + CSQI*S1R
	return AIR, AII, NZ, IERR
L140:
	AIR = -(ZR*S1R - ZI*S1I)
	AII = -(ZR*S1I + ZI*S1R)
	return AIR, AII, NZ, IERR
L150:
	S1R = S1R * SFAC
	S1I = S1I * SFAC
	if ID == 1 {
		goto L160
	}
	STR = S1R*CSQR - S1I*CSQI
	S1I = S1R*CSQI + S1I*CSQR
	S1R = STR
	AIR = S1R / SFAC
	AII = S1I / SFAC
	return AIR, AII, NZ, IERR
L160:
	STR = -(S1R*ZR - S1I*ZI)
	S1I = -(S1R*ZI + S1I*ZR)
	S1R = STR
	AIR = S1R / SFAC
	AII = S1I / SFAC
	return AIR, AII, NZ, IERR
L170:
	AA = 1000 * D1MACH[1]
	S1R = ZEROR
	S1I = ZEROI
	if ID == 1 {
		goto L190
	}
	if AZ <= AA {
		goto L180
	}
	S1R = C2 * ZR
	S1I = C2 * ZI
L180:
	AIR = C1 - S1R
	AII = -S1I
	return AIR, AII, NZ, IERR
L190:
	AIR = -C2
	AII = 0.0e0
	AA = math.Sqrt(AA)
	if AZ <= AA {
		goto L200
	}
	S1R = 0.5e0 * (ZR*ZR - ZI*ZI)
	S1I = ZR * ZI
L200:
	AIR = AIR + C1*S1R
	AII = AII + C1*S1I
	return AIR, AII, NZ, IERR
L210:
	NZ = 1
	AIR = ZEROR
	AII = ZEROI
	return AIR, AII, NZ, IERR
L270:
	NZ = 0
	IERR = 2
	return AIR, AII, NZ, IERR
L280:
	if NN == -1 {
		goto L270
	}
	NZ = 0
	IERR = 5
	return AIR, AII, NZ, IERR
L260:
	IERR = 4
	NZ = 0
	return AIR, AII, NZ, IERR
}

// ZASYI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z) >= 0.0 BY MEANS OF THE ASYMPTOTIC EXPANSION
// FOR LARGE CABS(Z) IN THE REGION CABS(Z) > MAX(RL,FNU*FNU/2). NZ=0 IS A NORMAL RETURN.
// NZ < 0 INDICATES AN OVERFLOW ON KODE=1.
func ZASYI(ZR float64, ZI float64, FNU float64, KODE int, N int, YR []float64, YI []float64, NZ int,
	RL float64, TOL float64, ELIM float64, ALIM float64) (float64, float64, float64, int, int, []float64, []float64, int,
	float64, float64, float64, float64) {

	const (
		ZEROR = 0.0e0
		ZEROI = 0.0e0
		CONER = 1.0e0
		CONEI = 0.0e0
		RTPI  = 0.159154943091895336e0
	)

	var AA, AEZ, AK, AK1I, AK1R, ARG, ARM, ATOL, AZ, BB, BK, CKI, CKR, CS1I, CS1R, CS2I, CS2R, CZI, CZR, DFNU, DKI,
		DKR, DNU2, EZI, EZR, FDN, P1I, P1R, RAZ, RTR1, RZI, RZR, S, SGN, SQK, STI, STR, S2I, S2R, TZI, TZR float64
	var I, IB, IL, INU, J, JL, K, KODED, M, NN int

	NZ = 0
	AZ = ZABS(ZR, ZI)
	ARM = 1000 * D1MACH[1]
	RTR1 = math.Sqrt(ARM)
	IL = MIN(2, N)
	DFNU = FNU + float64(float32(N-IL))

	// OVERFLOW TEST
	RAZ = 1.0e0 / AZ
	STR = ZR * RAZ
	STI = -ZI * RAZ
	AK1R = RTPI * STR * RAZ
	AK1I = RTPI * STI * RAZ
	AK1R, AK1I = ZSQRT(AK1R, AK1I)
	CZR = ZR
	CZI = ZI
	if KODE != 2 {
		goto L10
	}
	CZR = ZEROR
	CZI = ZI
L10:
	if math.Abs(CZR) > ELIM {
		goto L100
	}
	DNU2 = DFNU + DFNU
	KODED = 1
	if math.Abs(CZR) > ALIM && N > 2 {
		goto L20
	}
	KODED = 0
	STR, STI = ZEXP(CZR, CZI)
	AK1R, AK1I = ZMLT(AK1R, AK1I, STR, STI)
L20:
	FDN = 0.0e0
	if DNU2 > RTR1 {
		FDN = DNU2 * DNU2
	}
	EZR = ZR * 8.0e0
	EZI = ZI * 8.0e0

	// WHEN Z IS IMAGINARY, THE ERROR TEST MUST BE MADE RELATIVE TO THE FIRST RECIPROCAL POWER
	// SINCE THIS IS THE LEADING TERM OF THE EXPANSION FOR THE IMAGINARY PART.
	AEZ = 8.0e0 * AZ
	S = TOL / AEZ
	JL = int(float32(RL+RL)) + 2
	P1R = ZEROR
	P1I = ZEROI
	if ZI == 0.0e0 {
		goto L30
	}

	// CALCULATE EXP(PI*(0.5+FNU+N-IL)*I) TO MINIMIZE LOSSES OF SIGNIFICANCE WHEN FNU OR N IS LARGE
	INU = int(float32(FNU))
	ARG = (FNU - float64(float32(INU))) * math.Pi
	INU = INU + N - IL
	AK = -math.Sin(ARG)
	BK = math.Cos(ARG)
	if ZI < 0.0e0 {
		BK = -BK
	}
	P1R = AK
	P1I = BK
	if INU%2 == 0 {
		goto L30
	}
	P1R = -P1R
	P1I = -P1I
L30:
	for K = 1; K <= IL; K++ {
		SQK = FDN - 1.0e0
		ATOL = S * math.Abs(SQK)
		SGN = 1.0e0
		CS1R = CONER
		CS1I = CONEI
		CS2R = CONER
		CS2I = CONEI
		CKR = CONER
		CKI = CONEI
		AK = 0.0e0
		AA = 1.0e0
		BB = AEZ
		DKR = EZR
		DKI = EZI
		for J = 1; J <= JL; J++ {
			STR, STI = ZDIV(CKR, CKI, DKR, DKI)
			CKR = STR * SQK
			CKI = STI * SQK
			CS2R = CS2R + CKR
			CS2I = CS2I + CKI
			SGN = -SGN
			CS1R = CS1R + CKR*SGN
			CS1I = CS1I + CKI*SGN
			DKR = DKR + EZR
			DKI = DKI + EZI
			AA = AA * math.Abs(SQK) / BB
			BB = BB + AEZ
			AK = AK + 8.0e0
			SQK = SQK - AK
			if AA <= ATOL {
				goto L50
			}
		}
		goto L110
	L50:
		S2R = CS1R
		S2I = CS1I
		if ZR+ZR >= ELIM {
			goto L60
		}
		TZR = ZR + ZR
		TZI = ZI + ZI
		STR, STI = ZEXP(-TZR, -TZI)
		STR, STI = ZMLT(STR, STI, P1R, P1I)
		STR, STI = ZMLT(STR, STI, CS2R, CS2I)
		S2R = S2R + STR
		S2I = S2I + STI
	L60:
		FDN = FDN + 8.0e0*DFNU + 4.0e0
		P1R = -P1R
		P1I = -P1I
		M = N - IL + K
		YR[M] = S2R*AK1R - S2I*AK1I
		YI[M] = S2R*AK1I + S2I*AK1R
	}

	if N <= 2 {
		return ZR, ZI, FNU, KODE, N, YR, YI, NZ, RL, TOL, ELIM, ALIM
	}
	NN = N
	K = NN - 2
	AK = float64(float32(K))
	STR = ZR * RAZ
	STI = -ZI * RAZ
	RZR = (STR + STR) * RAZ
	RZI = (STI + STI) * RAZ
	IB = 3
	for I = IB; I <= NN; I++ {
		YR[K] = (AK+FNU)*(RZR*YR[K+1]-RZI*YI[K+1]) + YR[K+2]
		YI[K] = (AK+FNU)*(RZR*YI[K+1]+RZI*YR[K+1]) + YI[K+2]
		AK = AK - 1.0e0
		K = K - 1
	}
	if KODED == 0 {
		return ZR, ZI, FNU, KODE, N, YR, YI, NZ, RL, TOL, ELIM, ALIM
	}
	CKR, CKI = ZEXP(CZR, CZI)

	for I = 1; I <= NN; I++ {
		STR = YR[I]*CKR - YI[I]*CKI
		YI[I] = YR[I]*CKI + YI[I]*CKR
		YR[I] = STR
	}
	return ZR, ZI, FNU, KODE, N, YR, YI, NZ, RL, TOL, ELIM, ALIM
L100:
	NZ = -1
	return ZR, ZI, FNU, KODE, N, YR, YI, NZ, RL, TOL, ELIM, ALIM
L110:
	NZ = -2
	return ZR, ZI, FNU, KODE, N, YR, YI, NZ, RL, TOL, ELIM, ALIM
}

// ZBKNU COMPUTES THE K BESSEL FUNCTION IN THE RIGHT HALF Z PLANE.
func ZBKNU(ZR float64, ZI float64, FNU float64, KODE int, N int,
	YR []float64, YI []float64, NZ int, TOL float64, ELIM float64, ALIM float64) (float64, float64, float64, int, int,
	[]float64, []float64, int, float64, float64, float64) {

	const (
		DPI    = 3.14159265358979324e0
		RTHPI  = 1.25331413731550025e0
		SPI    = 1.90985931710274403e0
		HPI    = 1.57079632679489662e0
		FPI    = 1.89769999331517738e0
		TTH    = 6.66666666666666666e-01
		KMAX   = 30
		CZEROR = 0.0e0
		CZEROI = 0.0e0
		CONER  = 1.0e0
		CONEI  = 0.0e0
		CTWOR  = 2.0e0
		R1     = 2.0e0
	)

	var AA, AK, ASCLE, A1, A2, BB, BK, CAZ, CBI, CBR, CCHI, CCHR, CKI, CKR, COEFI, COEFR, CRSCR, CSCLR,
		CSHI, CSHR, CSI, CSR, CZI, CZR, DNU, DNU2, ETEST, FC, FHS, FI, FK, FKS, FMUI, FMUR, FR,
		G1, G2, PI, PR, PTI, PTR, P1I, P1R, P2I, P2M, P2R, QI, QR, RAK, RCAZ, RZI, RZR, S, SMUI, SMUR,
		STI, STR, S1I, S1R, S2I, S2R, TM, T1, T2, ELM, CELMR, ZDR, ZDI, AS, ALAS, HELIM float64

	var I, IFLAG, INU, K, KFLAG, KK, KODED, INUB, J, IC, NW int

	var CSSR, CSRR, BRY [4]float64
	var CYR, CYI [3]float64

	var CC = [9]float64{math.NaN(),
		5.77215664901532861e-01,
		-4.20026350340952355e-02,
		-4.21977345555443367e-02,
		7.21894324666309954e-03,
		-2.15241674114950973e-04,
		-2.01348547807882387e-05,
		1.13302723198169588e-06,
		6.11609510448141582e-09}

	CAZ = ZABS(ZR, ZI)
	CSCLR = 1.0e0 / TOL
	CRSCR = TOL
	CSSR[1] = CSCLR
	CSSR[2] = 1.0e0
	CSSR[3] = CRSCR
	CSRR[1] = CRSCR
	CSRR[2] = 1.0e0
	CSRR[3] = CSCLR
	BRY[1] = 1000.0 * D1MACH[1] / TOL
	BRY[2] = 1.0e0 / BRY[1]
	BRY[3] = D1MACH[2]
	NZ = 0
	IFLAG = 0
	KODED = KODE
	RCAZ = 1.0e0 / CAZ
	STR = ZR * RCAZ
	STI = -ZI * RCAZ
	RZR = (STR + STR) * RCAZ
	RZI = (STI + STI) * RCAZ
	INU = int(float32(FNU + 0.5e0))
	DNU = FNU - float64(float32(INU))
	if math.Abs(DNU) == 0.5e0 {
		goto L110
	}
	DNU2 = 0.0e0
	if math.Abs(DNU) > TOL {
		DNU2 = DNU * DNU
	}
	if CAZ > R1 {
		goto L110
	}

	// SERIES FOR CABS(Z) <= R1
	FC = 1.0e0
	SMUR, SMUI = ZLOG(RZR, RZI)
	FMUR = SMUR * DNU
	FMUI = SMUI * DNU
	CSHR, CSHI, CCHR, CCHI = ZSHCH(FMUR, FMUI)
	if DNU == 0.0e0 {
		goto L10
	}
	FC = DNU * DPI
	FC = FC / math.Sin(FC)
	SMUR = CSHR / DNU
	SMUI = CSHI / DNU
L10:
	A2 = 1.0e0 + DNU

	// GAM(1-Z)*GAM(1+Z)=PI*Z/SIN(PI*Z), T1=1/GAM(1-DNU), T2=1/GAM(1+DNU)
	T2 = math.Exp(-DGAMLN(A2))
	T1 = 1.0e0 / (T2 * FC)
	if math.Abs(DNU) > 0.1e0 {
		goto L40
	}

	// SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
	AK = 1.0e0
	S = CC[1]
	for K = 2; K <= 8; K++ {
		AK = AK * DNU2
		TM = CC[K] * AK
		S = S + TM
		if math.Abs(TM) < TOL {
			goto L30
		}
	}
L30:
	G1 = -S
	goto L50
L40:
	G1 = (T1 - T2) / (DNU + DNU)
L50:
	G2 = (T1 + T2) * 0.5e0
	FR = FC * (CCHR*G1 + SMUR*G2)
	FI = FC * (CCHI*G1 + SMUI*G2)
	STR, STI = ZEXP(FMUR, FMUI)
	PR = 0.5e0 * STR / T2
	PI = 0.5e0 * STI / T2
	PTR, PTI = ZDIV(0.5e0, 0.0e0, STR, STI)
	QR = PTR / T1
	QI = PTI / T1
	S1R = FR
	S1I = FI
	S2R = PR
	S2I = PI
	AK = 1.0e0
	A1 = 1.0e0
	CKR = CONER
	CKI = CONEI
	BK = 1.0e0 - DNU2
	if INU > 0 || N > 1 {
		goto L80
	}

	// GENERATE K(FNU,Z), 0.0e0  <=  FNU  <  0.5e0 AND N=1
	if CAZ < TOL {
		goto L70
	}
	CZR, CZI = ZMLT(ZR, ZI, ZR, ZI)
	CZR = 0.25e0 * CZR
	CZI = 0.25e0 * CZI
	T1 = 0.25e0 * CAZ * CAZ
L60:
	FR = (FR*AK + PR + QR) / BK
	FI = (FI*AK + PI + QI) / BK
	STR = 1.0e0 / (AK - DNU)
	PR = PR * STR
	PI = PI * STR
	STR = 1.0e0 / (AK + DNU)
	QR = QR * STR
	QI = QI * STR
	STR = CKR*CZR - CKI*CZI
	RAK = 1.0e0 / AK
	CKI = (CKR*CZI + CKI*CZR) * RAK
	CKR = STR * RAK
	S1R = CKR*FR - CKI*FI + S1R
	S1I = CKR*FI + CKI*FR + S1I
	A1 = A1 * T1 * RAK
	BK = BK + AK + AK + 1.0e0
	AK = AK + 1.0e0
	if A1 > TOL {
		goto L60
	}
L70:
	YR[1] = S1R
	YI[1] = S1I
	if KODED == 1 {
		return ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM, ALIM
	}
	STR, STI = ZEXP(ZR, ZI)
	YR[1], YI[1] = ZMLT(S1R, S1I, STR, STI)
	return ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM, ALIM

	// GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE
L80:
	if CAZ < TOL {
		goto L100
	}
	CZR, CZI = ZMLT(ZR, ZI, ZR, ZI)
	CZR = 0.25e0 * CZR
	CZI = 0.25e0 * CZI
	T1 = 0.25e0 * CAZ * CAZ
L90:
	FR = (FR*AK + PR + QR) / BK
	FI = (FI*AK + PI + QI) / BK
	STR = 1.0e0 / (AK - DNU)
	PR = PR * STR
	PI = PI * STR
	STR = 1.0e0 / (AK + DNU)
	QR = QR * STR
	QI = QI * STR
	STR = CKR*CZR - CKI*CZI
	RAK = 1.0e0 / AK
	CKI = (CKR*CZI + CKI*CZR) * RAK
	CKR = STR * RAK
	S1R = CKR*FR - CKI*FI + S1R
	S1I = CKR*FI + CKI*FR + S1I
	STR = PR - FR*AK
	STI = PI - FI*AK
	S2R = CKR*STR - CKI*STI + S2R
	S2I = CKR*STI + CKI*STR + S2I
	A1 = A1 * T1 * RAK
	BK = BK + AK + AK + 1.0e0
	AK = AK + 1.0e0
	if A1 > TOL {
		goto L90
	}
L100:
	KFLAG = 2
	A1 = FNU + 1.0e0
	AK = A1 * math.Abs(SMUR)
	if AK > ALIM {
		KFLAG = 3
	}
	STR = CSSR[KFLAG]
	P2R = S2R * STR
	P2I = S2I * STR
	S2R, S2I = ZMLT(P2R, P2I, RZR, RZI)
	S1R = S1R * STR
	S1I = S1I * STR
	if KODED == 1 {
		goto L210
	}
	FR, FI = ZEXP(ZR, ZI)
	S1R, S1I = ZMLT(S1R, S1I, FR, FI)
	S2R, S2I = ZMLT(S2R, S2I, FR, FI)
	goto L210

	// IFLAG=0 MEANS NO UNDERFLOW OCCURRED
	// IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH
	// KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD
	// RECURSION
L110:
	STR, STI = ZSQRT(ZR, ZI)
	COEFR, COEFI = ZDIV(RTHPI, CZEROI, STR, STI)
	KFLAG = 2
	if KODED == 2 {
		goto L120
	}
	if ZR > ALIM {
		goto L290
	}

	STR = math.Exp(-ZR) * CSSR[KFLAG]
	STI = -STR * math.Sin(ZI)
	STR = STR * math.Cos(ZI)
	COEFR, COEFI = ZMLT(COEFR, COEFI, STR, STI)
L120:
	if math.Abs(DNU) == 0.5e0 {
		goto L300
	}

	// MILLER ALGORITHM FOR CABS(Z) > R1
	AK = math.Cos(DPI * DNU)
	AK = math.Abs(AK)
	if AK == CZEROR {
		goto L300
	}
	FHS = math.Abs(0.25e0 - DNU2)
	if FHS == CZEROR {
		goto L300
	}

	// COMPUTE R2=F(E). if CABS(Z) >= R2, USE FORWARD RECURRENCE TO
	// DETERMINE THE BACKWARD INDEX K. R2=F(E) IS A STRAIGHT LINE ON
	// 12 <= E <= 60. E IS COMPUTED FROM 2**(-E)=B**(1-I1MACH(14))=
	// TOL WHERE B IS THE BASE OF THE ARITHMETIC.
	T1 = float64(I1MACH[14] - 1)
	T1 = T1 * D1MACH[5] * 3.321928094e0
	T1 = math.Max(T1, 12.0e0)
	T1 = math.Min(T1, 60.0e0)
	T2 = TTH*T1 - 6.0e0
	if ZR != 0.0e0 {
		goto L130
	}
	T1 = HPI
	goto L140
L130:
	T1 = math.Atan(ZI / ZR)
	T1 = math.Abs(T1)
L140:
	if T2 > CAZ {
		goto L170
	}

	// FORWARD RECURRENCE LOOP WHEN CABS(Z) >= R2
	ETEST = AK / (DPI * CAZ * TOL)
	FK = CONER
	if ETEST < CONER {
		goto L180
	}
	FKS = CTWOR
	CKR = CAZ + CAZ + CTWOR
	P1R = CZEROR
	P2R = CONER
	for I = 1; I <= KMAX; I++ {
		AK = FHS / FKS
		CBR = CKR / (FK + CONER)
		PTR = P2R
		P2R = CBR*P2R - P1R*AK
		P1R = PTR
		CKR = CKR + CTWOR
		FKS = FKS + FK + FK + CTWOR
		FHS = FHS + FK + FK
		FK = FK + CONER
		STR = math.Abs(P2R) * FK
		if ETEST < STR {
			goto L160
		}
	}
	goto L310
L160:
	FK = FK + SPI*T1*math.Sqrt(T2/CAZ)
	FHS = math.Abs(0.25e0 - DNU2)
	goto L180
L170:
	// COMPUTE BACKWARD INDEX K FOR CABS(Z) < R2
	A2 = math.Sqrt(CAZ)
	AK = FPI * AK / (TOL * math.Sqrt(A2))
	AA = 3.0e0 * T1 / (1.0e0 + CAZ)
	BB = 14.7e0 * T1 / (28.0e0 + CAZ)
	AK = (math.Log(AK) + CAZ*math.Cos(AA)/(1.0e0+0.008e0*CAZ)) / math.Cos(BB)
	FK = 0.12125e0*AK*AK/CAZ + 1.5e0
L180:
	// BACKWARD RECURRENCE LOOP FOR MILLER ALGORITHM
	K = int(float32(FK))
	FK = float64(float32(K))
	FKS = FK * FK
	P1R = CZEROR
	P1I = CZEROI
	P2R = TOL
	P2I = CZEROI
	CSR = P2R
	CSI = P2I
	for I = 1; I <= K; I++ {
		A1 = FKS - FK
		AK = (FKS + FK) / (A1 + FHS)
		RAK = 2.0e0 / (FK + CONER)
		CBR = (FK + ZR) * RAK
		CBI = ZI * RAK
		PTR = P2R
		PTI = P2I
		P2R = (PTR*CBR - PTI*CBI - P1R) * AK
		P2I = (PTI*CBR + PTR*CBI - P1I) * AK
		P1R = PTR
		P1I = PTI
		CSR = CSR + P2R
		CSI = CSI + P2I
		FKS = A1 - FK + CONER
		FK = FK - CONER
	}

	// COMPUTE (P2/CS)=(P2/CABS(CS))*(CONJG(CS)/CABS(CS)) FOR BETTER
	// SCALING
	TM = ZABS(CSR, CSI)
	PTR = 1.0e0 / TM
	S1R = P2R * PTR
	S1I = P2I * PTR
	CSR = CSR * PTR
	CSI = -CSI * PTR
	STR, STI = ZMLT(COEFR, COEFI, S1R, S1I)
	S1R, S1I = ZMLT(STR, STI, CSR, CSI)
	if INU > 0 || N > 1 {
		goto L200
	}
	ZDR = ZR
	ZDI = ZI
	if IFLAG == 1 {
		goto L270
	}
	goto L240
L200:
	// COMPUTE P1/P2=(P1/CABS(P2)*CONJG(P2)/CABS(P2) FOR SCALING
	TM = ZABS(P2R, P2I)
	PTR = 1.0e0 / TM
	P1R = P1R * PTR
	P1I = P1I * PTR
	P2R = P2R * PTR
	P2I = -P2I * PTR
	PTR, PTI = ZMLT(P1R, P1I, P2R, P2I)
	STR = DNU + 0.5e0 - PTR
	STI = -PTI
	STR, STI = ZDIV(STR, STI, ZR, ZI)
	STR = STR + 1.0e0
	S2R, S2I = ZMLT(STR, STI, S1R, S1I)

	// FORWARD RECURSION ON THE THREE TERM RECURSION WITH RELATION WITH
	// SCALING NEAR EXPONENT EXTREMES ON KFLAG=1 OR KFLAG=3
L210:
	STR = DNU + 1.0e0
	CKR = STR * RZR
	CKI = STR * RZI
	if N == 1 {
		INU = INU - 1
	}
	if INU > 0 {
		goto L220
	}
	if N > 1 {
		goto L215
	}
	S1R = S2R
	S1I = S2I
L215:
	ZDR = ZR
	ZDI = ZI
	if IFLAG == 1 {
		goto L270
	}
	goto L240
L220:
	INUB = 1
	if IFLAG == 1 {
		goto L261
	}
L225:
	P1R = CSRR[KFLAG]
	ASCLE = BRY[KFLAG]
	for I = INUB; I <= INU; I++ {
		STR = S2R
		STI = S2I
		S2R = CKR*STR - CKI*STI + S1R
		S2I = CKR*STI + CKI*STR + S1I
		S1R = STR
		S1I = STI
		CKR = CKR + RZR
		CKI = CKI + RZI
		if KFLAG >= 3 {
			continue
		}
		P2R = S2R * P1R
		P2I = S2I * P1R
		STR = math.Abs(P2R)
		STI = math.Abs(P2I)
		P2M = math.Max(STR, STI)
		if P2M <= ASCLE {
			continue
		}
		KFLAG = KFLAG + 1
		ASCLE = BRY[KFLAG]
		S1R = S1R * P1R
		S1I = S1I * P1R
		S2R = P2R
		S2I = P2I
		STR = CSSR[KFLAG]
		S1R = S1R * STR
		S1I = S1I * STR
		S2R = S2R * STR
		S2I = S2I * STR
		P1R = CSRR[KFLAG]
	}

	if N != 1 {
		goto L240
	}
	S1R = S2R
	S1I = S2I
L240:
	STR = CSRR[KFLAG]
	YR[1] = S1R * STR
	YI[1] = S1I * STR
	if N == 1 {
		return ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM, ALIM
	}
	YR[2] = S2R * STR
	YI[2] = S2I * STR
	if N == 2 {
		return ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM, ALIM
	}
	KK = 2
L250:
	KK = KK + 1
	if KK > N {
		return ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM, ALIM
	}
	P1R = CSRR[KFLAG]
	ASCLE = BRY[KFLAG]
	for I = KK; I <= N; I++ {
		P2R = S2R
		P2I = S2I
		S2R = CKR*P2R - CKI*P2I + S1R
		S2I = CKI*P2R + CKR*P2I + S1I
		S1R = P2R
		S1I = P2I
		CKR = CKR + RZR
		CKI = CKI + RZI
		P2R = S2R * P1R
		P2I = S2I * P1R
		YR[I] = P2R
		YI[I] = P2I
		if KFLAG >= 3 {
			continue
		}
		STR = math.Abs(P2R)
		STI = math.Abs(P2I)
		P2M = math.Max(STR, STI)
		if P2M <= ASCLE {
			continue
		}
		KFLAG = KFLAG + 1
		ASCLE = BRY[KFLAG]
		S1R = S1R * P1R
		S1I = S1I * P1R
		S2R = P2R
		S2I = P2I
		STR = CSSR[KFLAG]
		S1R = S1R * STR
		S1I = S1I * STR
		S2R = S2R * STR
		S2I = S2I * STR
		P1R = CSRR[KFLAG]
	}
	return ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM, ALIM

	// IFLAG=1 CASES, FORWARD RECURRENCE ON SCALED VALUES ON UNDERFLOW
L261:
	HELIM = 0.5e0 * ELIM
	ELM = math.Exp(-ELIM)
	CELMR = ELM
	ASCLE = BRY[1]
	ZDR = ZR
	ZDI = ZI
	IC = -1
	J = 2
	for I = 1; I <= INU; I++ {
		STR = S2R
		STI = S2I
		S2R = STR*CKR - STI*CKI + S1R
		S2I = STI*CKR + STR*CKI + S1I
		S1R = STR
		S1I = STI
		CKR = CKR + RZR
		CKI = CKI + RZI
		AS = ZABS(S2R, S2I)
		ALAS = math.Log(AS)
		P2R = -ZDR + ALAS
		if P2R < -ELIM {
			goto L263
		}
		STR, STI = ZLOG(S2R, S2I)
		P2R = -ZDR + STR
		P2I = -ZDI + STI
		P2M = math.Exp(P2R) / TOL
		P1R = P2M * math.Cos(P2I)
		P1I = P2M * math.Sin(P2I)
		P1R, P1I, NW, ASCLE, TOL = ZUCHK(P1R, P1I, NW, ASCLE, TOL)
		if NW != 0 {
			goto L263
		}
		J = 3 - J
		CYR[J] = P1R
		CYI[J] = P1I
		if IC == (I - 1) {
			goto L264
		}
		IC = I
		continue
	L263:
		if ALAS < HELIM {
			continue
		}
		ZDR = ZDR - ELIM
		S1R = S1R * CELMR
		S1I = S1I * CELMR
		S2R = S2R * CELMR
		S2I = S2I * CELMR
	}
	if N != 1 {
		goto L270
	}
	S1R = S2R
	S1I = S2I
	goto L270
L264:
	KFLAG = 1
	INUB = I + 1
	S2R = CYR[J]
	S2I = CYI[J]
	J = 3 - J
	S1R = CYR[J]
	S1I = CYI[J]
	if INUB <= INU {
		goto L225
	}
	if N != 1 {
		goto L240
	}
	S1R = S2R
	S1I = S2I
	goto L240
L270:
	YR[1] = S1R
	YI[1] = S1I
	if N == 1 {
		goto L280
	}
	YR[2] = S2R
	YI[2] = S2I
L280:
	ASCLE = BRY[1]
	ZDR, ZDI, FNU, N, YR, YI, NZ, RZR, RZI, ASCLE, TOL, ELIM = ZKSCL(ZDR, ZDI, FNU, N, YR, YI, NZ, RZR, RZI, ASCLE, TOL, ELIM)
	INU = N - NZ
	if INU <= 0 {
		return ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM, ALIM
	}
	KK = NZ + 1
	S1R = YR[KK]
	S1I = YI[KK]
	YR[KK] = S1R * CSRR[1]
	YI[KK] = S1I * CSRR[1]
	if INU == 1 {
		return ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM, ALIM
	}
	KK = NZ + 2
	S2R = YR[KK]
	S2I = YI[KK]
	YR[KK] = S2R * CSRR[1]
	YI[KK] = S2I * CSRR[1]
	if INU == 2 {
		return ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM, ALIM
	}
	T2 = FNU + float64(float32(KK-1))
	CKR = T2 * RZR
	CKI = T2 * RZI
	KFLAG = 1
	goto L250
L290:
	// SCALE BY DEXP(Z), IFLAG = 1 CASES
	KODED = 2
	IFLAG = 1
	KFLAG = 2
	goto L120

	// FNU=HALF ODD INTEGER CASE, DNU=-0.5
L300:
	S1R = COEFR
	S1I = COEFI
	S2R = COEFR
	S2I = COEFI
	goto L210

L310:
	NZ = -2
	return ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM, ALIM

}

// SET K FUNCTIONS TO ZERO ON UNDERFLOW, CONTINUE RECURRENCE ON SCALED FUNCTIONS
// UNTIL TWO MEMBERS COME ON SCALE, THEN RETURN WITH MIN(NZ+2,N) VALUES SCALED BY 1/TOL.
func ZKSCL(ZRR float64, ZRI float64, FNU float64, N int, YR []float64, YI []float64, NZ int,
	RZR float64, RZI float64, ASCLE float64, TOL float64, ELIM float64) (float64, float64, float64, int, []float64, []float64, int,
	float64, float64, float64, float64, float64) {

	const (
		ZEROR = 0.0e0
		ZEROI = 0.0e0
	)

	var ACS, AS, CSI, CSR, FN, STR, S1I, S1R, S2I, S2R, CKR, CKI, ZDR, ZDI, CELMR, ELM, HELIM, ALAS float64
	var I, IC, KK, NN, NW int
	var CYR, CYI [3]float64

	NZ = 0
	IC = 0
	NN = MIN(2, N)
	for I = 1; I <= NN; I++ {
		S1R = YR[I]
		S1I = YI[I]
		CYR[I] = S1R
		CYI[I] = S1I
		AS = ZABS(S1R, S1I)
		ACS = -ZRR + math.Log(AS)
		NZ = NZ + 1
		YR[I] = ZEROR
		YI[I] = ZEROI
		if ACS < -ELIM {
			continue
		}
		CSR, CSI = ZLOG(S1R, S1I)
		CSR = CSR - ZRR
		CSI = CSI - ZRI
		STR = math.Exp(CSR) / TOL
		CSR = STR * math.Cos(CSI)
		CSI = STR * math.Sin(CSI)
		CSR, CSI, NW, ASCLE, TOL = ZUCHK(CSR, CSI, NW, ASCLE, TOL)
		if NW != 0 {
			continue
		}
		YR[I] = CSR
		YI[I] = CSI
		IC = I
		NZ = NZ - 1
	}
	if N == 1 {
		return ZRR, ZRI, FNU, N, YR, YI, NZ, RZR, RZI, ASCLE, TOL, ELIM
	}
	if IC > 1 {
		goto L20
	}
	YR[1] = ZEROR
	YI[1] = ZEROI
	NZ = 2
L20:
	if N == 2 {
		return ZRR, ZRI, FNU, N, YR, YI, NZ, RZR, RZI, ASCLE, TOL, ELIM
	}
	if NZ == 0 {
		return ZRR, ZRI, FNU, N, YR, YI, NZ, RZR, RZI, ASCLE, TOL, ELIM
	}
	FN = FNU + 1.0e0
	CKR = FN * RZR
	CKI = FN * RZI
	S1R = CYR[1]
	S1I = CYI[1]
	S2R = CYR[2]
	S2I = CYI[2]
	HELIM = 0.5e0 * ELIM
	ELM = math.Exp(-ELIM)
	CELMR = ELM
	ZDR = ZRR
	ZDI = ZRI

	// FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE RECURRENCE IF
	// S2 GETS LARGER THAN EXP(ELIM/2)
	for I = 3; I <= NN; I++ {

		KK = I
		CSR = S2R
		CSI = S2I
		S2R = CKR*CSR - CKI*CSI + S1R
		S2I = CKI*CSR + CKR*CSI + S1I
		S1R = CSR
		S1I = CSI
		CKR = CKR + RZR
		CKI = CKI + RZI
		AS = ZABS(S2R, S2I)
		ALAS = math.Log(AS)
		ACS = -ZDR + ALAS
		NZ = NZ + 1
		YR[I] = ZEROR
		YI[I] = ZEROI
		if ACS < -ELIM {
			goto L25
		}
		CSR, CSI = ZLOG(S2R, S2I)
		CSR = CSR - ZDR
		CSI = CSI - ZDI
		STR = math.Exp(CSR) / TOL
		CSR = STR * math.Cos(CSI)
		CSI = STR * math.Sin(CSI)
		CSR, CSI, NW, ASCLE, TOL = ZUCHK(CSR, CSI, NW, ASCLE, TOL)
		if NW != 0 {
			goto L25
		}
		YR[I] = CSR
		YI[I] = CSI
		NZ = NZ - 1
		if IC == KK-1 {
			goto L40
		}
		IC = KK
		continue
	L25:
		if ALAS < HELIM {
			continue
		}
		ZDR = ZDR - ELIM
		S1R = S1R * CELMR
		S1I = S1I * CELMR
		S2R = S2R * CELMR
		S2I = S2I * CELMR
	}
	NZ = N
	if IC == N {
		NZ = N - 1
	}
	goto L45
L40:
	NZ = KK - 2
L45:
	for I = 1; I <= NZ; I++ {
		YR[I] = ZEROR
		YI[I] = ZEROI
	}
	return ZRR, ZRI, FNU, N, YR, YI, NZ, RZR, RZI, ASCLE, TOL, ELIM
}

// ZMLRI COMPUTES THE I BESSEL FUNCTION FOR RE(Z) >= 0.0 BY THE MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES.
func ZMLRI(ZR float64, ZI float64, FNU float64, KODE int, N int, YR []float64, YI []float64,
	NZ int, TOL float64) (float64, float64, float64, int, int, []float64, []float64, int, float64) {

	const (
		ZEROR = 0.0e0
		ZEROI = 0.0e0
		CONER = 1.0e0
		CONEI = 0.0e0
	)
	var ACK, AK, AP, AT, AZ, BK, CKI, CKR, CNORMI, CNORMR, FKAP, FKK, FLAM, FNF, PTI, PTR, P1I,
		P1R, P2I, P2R, RAZ, RHO, RHO2, RZI, RZR, SCLE, STI, STR, SUMI, SUMR, TFNF, TST float64
	var I, IAZ, IFNU, INU, ITIME, K, KK, KM, M int

	SCLE = D1MACH[1] / TOL
	NZ = 0
	AZ = ZABS(ZR, ZI)
	IAZ = int(float32(AZ))
	IFNU = int(float32(FNU))
	INU = IFNU + N - 1
	AT = float64(float32(IAZ)) + 1.0e0
	RAZ = 1.0e0 / AZ
	STR = ZR * RAZ
	STI = -ZI * RAZ
	CKR = STR * AT * RAZ
	CKI = STI * AT * RAZ
	RZR = (STR + STR) * RAZ
	RZI = (STI + STI) * RAZ
	P1R = ZEROR
	P1I = ZEROI
	P2R = CONER
	P2I = CONEI
	ACK = (AT + 1.0e0) * RAZ
	RHO = ACK + math.Sqrt(ACK*ACK-1.0e0)
	RHO2 = RHO * RHO
	TST = (RHO2 + RHO2) / ((RHO2 - 1.0e0) * (RHO - 1.0e0))
	TST = TST / TOL

	// COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES
	AK = AT
	for I = 1; I <= 80; I++ {

		PTR = P2R
		PTI = P2I
		P2R = P1R - (CKR*PTR - CKI*PTI)
		P2I = P1I - (CKI*PTR + CKR*PTI)
		P1R = PTR
		P1I = PTI
		CKR = CKR + RZR
		CKI = CKI + RZI
		AP = ZABS(P2R, P2I)
		if AP > TST*AK*AK {
			goto L20
		}
		AK = AK + 1.0e0
	}
	goto L110
L20:
	I = I + 1
	K = 0
	if INU < IAZ {
		goto L40
	}

	// COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS
	P1R = ZEROR
	P1I = ZEROI
	P2R = CONER
	P2I = CONEI
	AT = float64(float32(INU)) + 1.0e0
	STR = ZR * RAZ
	STI = -ZI * RAZ
	CKR = STR * AT * RAZ
	CKI = STI * AT * RAZ
	ACK = AT * RAZ
	TST = math.Sqrt(ACK / TOL)
	ITIME = 1
	for K = 1; K <= 80; K++ {
		PTR = P2R
		PTI = P2I
		P2R = P1R - (CKR*PTR - CKI*PTI)
		P2I = P1I - (CKR*PTI + CKI*PTR)
		P1R = PTR
		P1I = PTI
		CKR = CKR + RZR
		CKI = CKI + RZI
		AP = ZABS(P2R, P2I)
		if AP < TST {
			continue
		}
		if ITIME == 2 {
			goto L40
		}
		ACK = ZABS(CKR, CKI)
		FLAM = ACK + math.Sqrt(ACK*ACK-1.0e0)
		FKAP = AP / ZABS(P1R, P1I)
		RHO = math.Min(FLAM, FKAP)
		TST = TST * math.Sqrt(RHO/(RHO*RHO-1.0e0))
		ITIME = 2
	}
	goto L110
L40:
	// BACKWARD RECURRENCE AND SUM NORMALIZING RELATION
	K = K + 1
	KK = MAX(I+IAZ, K+INU)
	FKK = float64(float32(KK))
	P1R = ZEROR
	P1I = ZEROI

	// SCALE P2 AND SUM BY SCLE
	P2R = SCLE
	P2I = ZEROI
	FNF = FNU - float64(float32(IFNU))
	TFNF = FNF + FNF
	BK = DGAMLN(FKK+TFNF+1.0e0) - DGAMLN(FKK+1.0e0) - DGAMLN(TFNF+1.0e0)
	BK = math.Exp(BK)
	SUMR = ZEROR
	SUMI = ZEROI
	KM = KK - INU
	for I = 1; I <= KM; I++ {
		PTR = P2R
		PTI = P2I
		P2R = P1R + (FKK+FNF)*(RZR*PTR-RZI*PTI)
		P2I = P1I + (FKK+FNF)*(RZI*PTR+RZR*PTI)
		P1R = PTR
		P1I = PTI
		AK = 1.0e0 - TFNF/(FKK+TFNF)
		ACK = BK * AK
		SUMR = SUMR + (ACK+BK)*P1R
		SUMI = SUMI + (ACK+BK)*P1I
		BK = ACK
		FKK = FKK - 1.0e0
	}
	YR[N] = P2R
	YI[N] = P2I
	if N == 1 {
		goto L70
	}
	for I = 2; I <= N; I++ {
		PTR = P2R
		PTI = P2I
		P2R = P1R + (FKK+FNF)*(RZR*PTR-RZI*PTI)
		P2I = P1I + (FKK+FNF)*(RZI*PTR+RZR*PTI)
		P1R = PTR
		P1I = PTI
		AK = 1.0e0 - TFNF/(FKK+TFNF)
		ACK = BK * AK
		SUMR = SUMR + (ACK+BK)*P1R
		SUMI = SUMI + (ACK+BK)*P1I
		BK = ACK
		FKK = FKK - 1.0e0
		M = N - I + 1
		YR[M] = P2R
		YI[M] = P2I
	}
L70:
	if IFNU <= 0 {
		goto L90
	}
	for I = 1; I <= IFNU; I++ {
		PTR = P2R
		PTI = P2I
		P2R = P1R + (FKK+FNF)*(RZR*PTR-RZI*PTI)
		P2I = P1I + (FKK+FNF)*(RZR*PTI+RZI*PTR)
		P1R = PTR
		P1I = PTI
		AK = 1.0e0 - TFNF/(FKK+TFNF)
		ACK = BK * AK
		SUMR = SUMR + (ACK+BK)*P1R
		SUMI = SUMI + (ACK+BK)*P1I
		BK = ACK
		FKK = FKK - 1.0e0
	}
L90:
	PTR = ZR
	PTI = ZI
	if KODE == 2 {
		PTR = ZEROR
	}
	STR, STI = ZLOG(RZR, RZI)
	P1R = -FNF*STR + PTR
	P1I = -FNF*STI + PTI
	AP = DGAMLN(1.0e0 + FNF)
	PTR = P1R - AP
	PTI = P1I

	// THE DIVISION CEXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW
	// IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES
	P2R = P2R + SUMR
	P2I = P2I + SUMI
	AP = ZABS(P2R, P2I)
	P1R = 1.0e0 / AP
	STR, STI = ZEXP(PTR, PTI)
	CKR = STR * P1R
	CKI = STI * P1R
	PTR = P2R * P1R
	PTI = -P2I * P1R
	CNORMR, CNORMI = ZMLT(CKR, CKI, PTR, PTI)
	for I = 1; I <= N; I++ {
		STR = YR[I]*CNORMR - YI[I]*CNORMI
		YI[I] = YR[I]*CNORMI + YI[I]*CNORMR
		YR[I] = STR
	}
	return ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL
L110:
	NZ = -2
	return ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL
}

// ZSERI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z) >= 0.0 BY
// MEANS OF THE POWER SERIES FOR LARGE CABS(Z) IN THE
// REGION CABS(Z) <= 2*SQRT(FNU+1). NZ=0 IS A NORMAL RETURN.
// NZ > 0 MEANS THAT THE LAST NZ COMPONENTS WERE SET TO ZERO
// DUE TO UNDERFLOW. NZ < 0 MEANS UNDERFLOW OCCURRED, BUT THE
// CONDITION CABS(Z) <= 2*SQRT(FNU+1) WAS VIOLATED AND THE
// COMPUTATION MUST BE COMPLETED IN ANOTHER ROUTINE WITH N=N-ABS(NZ).
func ZSERI(ZR float64, ZI float64, FNU float64, KODE int, N int, YR []float64, YI []float64, NZ int,
	TOL float64, ELIM float64, ALIM float64) (float64, float64, float64, int, int, []float64, []float64, int,
	float64, float64, float64) {

	const (
		ZEROR = 0.0e0
		ZEROI = 0.0e0
		CONER = 1.0e0
		CONEI = 0.0e0
	)

	// COMPLEX AK1,CK,COEF,CONE,CRSC,CSCL,CZ,CZERO,HZ,RZ,S1,S2,Y,Z
	var AA, ACZ, AK, AK1I, AK1R, ARM, ASCLE, ATOL, AZ, CKI, CKR, COEFI, COEFR, CRSCR, CZI, CZR, DFNU,
		FNUP, HZI, HZR, RAZ, RS, RTR1, RZI, RZR, S, SS, STI, STR, S1I, S1R, S2I, S2R float64
	var I, IB, IFLAG, IL, K, L, M, NN, NW int
	WR := []float64{math.NaN(), 0, 0}
	WI := []float64{math.NaN(), 0, 0}

	NZ = 0
	AZ = ZABS(ZR, ZI)
	if AZ == 0.0e0 {
		goto L160
	}
	ARM = 1000.0 * D1MACH[1]
	RTR1 = math.Sqrt(ARM)
	CRSCR = 1.0e0
	IFLAG = 0
	if AZ < ARM {
		goto L150
	}
	HZR = 0.5e0 * ZR
	HZI = 0.5e0 * ZI
	CZR = ZEROR
	CZI = ZEROI
	if AZ <= RTR1 {
		goto L10
	}
	CZR, CZI = ZMLT(HZR, HZI, HZR, HZI)
L10:
	ACZ = ZABS(CZR, CZI)
	NN = N
	CKR, CKI = ZLOG(HZR, HZI)
L20:
	DFNU = FNU + float64(float32(NN-1))
	FNUP = DFNU + 1.0e0

	// UNDERFLOW TEST
	AK1R = CKR * DFNU
	AK1I = CKI * DFNU
	AK = DGAMLN(FNUP)
	AK1R = AK1R - AK
	if KODE == 2 {
		AK1R = AK1R - ZR
	}
	if AK1R > -ELIM {
		goto L40
	}
L30:
	NZ = NZ + 1
	YR[NN] = ZEROR
	YI[NN] = ZEROI
	if ACZ > DFNU {
		goto L190
	}
	NN = NN - 1
	if NN == 0 {
		return ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM, ALIM
	}
	goto L20
L40:
	if AK1R > -ALIM {
		goto L50
	}
	IFLAG = 1
	SS = 1.0e0 / TOL
	CRSCR = TOL
	ASCLE = ARM * SS
L50:
	AA = math.Exp(AK1R)
	if IFLAG == 1 {
		AA = AA * SS
	}
	COEFR = AA * math.Cos(AK1I)
	COEFI = AA * math.Sin(AK1I)
	ATOL = TOL * ACZ / FNUP
	IL = MIN(2, NN)
	for I = 1; I <= IL; I++ {
		DFNU = FNU + float64(float32(NN-I))
		FNUP = DFNU + 1.0e0
		S1R = CONER
		S1I = CONEI
		if ACZ < TOL*FNUP {
			goto L70
		}
		AK1R = CONER
		AK1I = CONEI
		AK = FNUP + 2.0e0
		S = FNUP
		AA = 2.0e0
	L60:
		RS = 1.0e0 / S
		STR = AK1R*CZR - AK1I*CZI
		STI = AK1R*CZI + AK1I*CZR
		AK1R = STR * RS
		AK1I = STI * RS
		S1R = S1R + AK1R
		S1I = S1I + AK1I
		S = S + AK
		AK = AK + 2.0e0
		AA = AA * ACZ * RS
		if AA > ATOL {
			goto L60
		}
	L70:
		S2R = S1R*COEFR - S1I*COEFI
		S2I = S1R*COEFI + S1I*COEFR
		WR[I] = S2R
		WI[I] = S2I
		if IFLAG == 0 {
			goto L80
		}
		S2R, S2I, NW, ASCLE, TOL = ZUCHK(S2R, S2I, NW, ASCLE, TOL)
		if NW != 0 {
			goto L30
		}
	L80:
		M = NN - I + 1
		YR[M] = S2R * CRSCR
		YI[M] = S2I * CRSCR
		if I == IL {
			continue
		}
		STR, STI = ZDIV(COEFR, COEFI, HZR, HZI)
		COEFR = STR * DFNU
		COEFI = STI * DFNU
	}
	if NN <= 2 {
		return ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM, ALIM
	}
	K = NN - 2
	AK = float64(float32(K))
	RAZ = 1.0e0 / AZ
	STR = ZR * RAZ
	STI = -ZI * RAZ
	RZR = (STR + STR) * RAZ
	RZI = (STI + STI) * RAZ
	if IFLAG == 1 {
		goto L120
	}
	IB = 3
L100:
	for I = IB; I <= NN; I++ {
		YR[K] = (AK+FNU)*(RZR*YR[K+1]-RZI*YI[K+1]) + YR[K+2]
		YI[K] = (AK+FNU)*(RZR*YI[K+1]+RZI*YR[K+1]) + YI[K+2]
		AK = AK - 1.0e0
		K = K - 1
	}
	return ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM, ALIM

	// RECUR BACKWARD WITH SCALED VALUES
L120:
	// EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION ABOVE THE
	// UNDERFLOW LIMIT = ASCLE = D1MACH(1)*SS*1.0D+3
	S1R = WR[1]
	S1I = WI[1]
	S2R = WR[2]
	S2I = WI[2]
	for L = 3; L <= NN; L++ {
		CKR = S2R
		CKI = S2I
		S2R = S1R + (AK+FNU)*(RZR*CKR-RZI*CKI)
		S2I = S1I + (AK+FNU)*(RZR*CKI+RZI*CKR)
		S1R = CKR
		S1I = CKI
		CKR = S2R * CRSCR
		CKI = S2I * CRSCR
		YR[K] = CKR
		YI[K] = CKI
		AK = AK - 1.0e0
		K = K - 1
		if ZABS(CKR, CKI) > ASCLE {
			goto L140
		}
	}
	return ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM, ALIM
L140:
	IB = L + 1
	if IB > NN {
		return ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM, ALIM
	}
	goto L100
L150:
	NZ = N
	if FNU == 0.0e0 {
		NZ = NZ - 1
	}
L160:
	YR[1] = ZEROR
	YI[1] = ZEROI
	if FNU != 0.0e0 {
		goto L170
	}
	YR[1] = CONER
	YI[1] = CONEI
L170:
	if N == 1 {
		return ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM, ALIM
	}
	for I = 2; I <= N; I++ {
		YR[I] = ZEROR
		YI[I] = ZEROI
	}
	return ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM, ALIM

	// RETURN WITH NZ < 0 if CABS(Z*Z/4) > FNU+N-NZ-1 COMPLETE
	// THE CALCULATION IN CBINU WITH N=N-IABS(NZ)
L190:
	NZ = -NZ
	return ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM, ALIM

}

// ZUNIK COMPUTES PARAMETERS FOR THE UNIFORM ASYMPTOTIC
// EXPANSIONS OF THE I AND K FUNCTIONS ON IKFLG= 1 OR 2
// RESPECTIVELY BY
//
// W(FNU,ZR) = PHI*EXP(ZETA)*SUM
//
// WHERE       ZETA=-ZETA1 + ZETA2       OR
//                   ZETA1 - ZETA2
//
// THE FIRST CALL MUST HAVE INIT=0. SUBSEQUENT CALLS WITH THE
// SAME ZR AND FNU WILL return  THE I OR K FUNCTION ON IKFLG=
// 1 OR 2 WITH NO CHANGE IN INIT. CWRK IS A COMPLEX WORK
// ARRAY. IPMTR=0 COMPUTES ALL PARAMETERS. IPMTR=1 COMPUTES PHI,
// ZETA1,ZETA2.
func ZUNIK(ZRR float64, ZRI float64, FNU float64, IKFLG int, IPMTR int, TOL float64, INIT int, PHIR float64,
	PHII float64, ZETA1R float64, ZETA1I float64, ZETA2R float64, ZETA2I float64,
	SUMR float64, SUMI float64, CWRKR []float64, CWRKI []float64) (float64, float64, float64, int, int, float64, int, float64,
	float64, float64, float64, float64, float64, float64, float64, []float64, []float64) {

	const (
		ZEROR = 0.0e0
		ZEROI = 0.0e0
		CONER = 1.0e0
		CONEI = 0.0e0
	)

	var AC, CRFNI, CRFNR, RFN, SI, SR, SRI, SRR, STI, STR,
		TEST, TI, TR, T2I, T2R, ZNI, ZNR float64
	var I, J, K, L int
	var CON = [3]float64{math.NaN(), 3.98942280401432678e-01, 1.25331413731550025e+00}
	//var CWRKR, CWRKI [17]float64

	var C = [121]float64{math.NaN(),
		1.00000000000000000e+00, -2.08333333333333333e-01,
		1.25000000000000000e-01, 3.34201388888888889e-01,
		-4.01041666666666667e-01, 7.03125000000000000e-02,
		-1.02581259645061728e+00, 1.84646267361111111e+00,
		-8.91210937500000000e-01, 7.32421875000000000e-02,
		4.66958442342624743e+00, -1.12070026162229938e+01,
		8.78912353515625000e+00, -2.36408691406250000e+00,
		1.12152099609375000e-01, -2.82120725582002449e+01,
		8.46362176746007346e+01, -9.18182415432400174e+01,
		4.25349987453884549e+01, -7.36879435947963170e+00,
		2.27108001708984375e-01, 2.12570130039217123e+02,
		-7.65252468141181642e+02, 1.05999045252799988e+03,
		-6.99579627376132541e+02, 2.18190511744211590e+02,
		-2.64914304869515555e+01, 5.72501420974731445e-01,
		-1.91945766231840700e+03, 8.06172218173730938e+03,
		-1.35865500064341374e+04, 1.16553933368645332e+04,
		-5.30564697861340311e+03, 1.20090291321635246e+03,
		-1.08090919788394656e+02, 1.72772750258445740e+00,
		2.02042913309661486e+04, -9.69805983886375135e+04,
		1.92547001232531532e+05, -2.03400177280415534e+05,
		1.22200464983017460e+05, -4.11926549688975513e+04,
		7.10951430248936372e+03, -4.93915304773088012e+02,
		6.07404200127348304e+00, -2.42919187900551333e+05,
		1.31176361466297720e+06, -2.99801591853810675e+06,
		3.76327129765640400e+06, -2.81356322658653411e+06,
		1.26836527332162478e+06, -3.31645172484563578e+05,
		4.52187689813627263e+04, -2.49983048181120962e+03,
		2.43805296995560639e+01, 3.28446985307203782e+06,
		-1.97068191184322269e+07, 5.09526024926646422e+07,
		-7.41051482115326577e+07, 6.63445122747290267e+07,
		-3.75671766607633513e+07, 1.32887671664218183e+07,
		-2.78561812808645469e+06, 3.08186404612662398e+05,
		-1.38860897537170405e+04, 1.10017140269246738e+02,
		-4.93292536645099620e+07, 3.25573074185765749e+08,
		-9.39462359681578403e+08, 1.55359689957058006e+09,
		-1.62108055210833708e+09, 1.10684281682301447e+09,
		-4.95889784275030309e+08, 1.42062907797533095e+08,
		-2.44740627257387285e+07, 2.24376817792244943e+06,
		-8.40054336030240853e+04, 5.51335896122020586e+02,
		8.14789096118312115e+08, -5.86648149205184723e+09,
		1.86882075092958249e+10, -3.46320433881587779e+10,
		4.12801855797539740e+10, -3.30265997498007231e+10,
		1.79542137311556001e+10, -6.56329379261928433e+09,
		1.55927986487925751e+09, -2.25105661889415278e+08,
		1.73951075539781645e+07, -5.49842327572288687e+05,
		3.03809051092238427e+03, -1.46792612476956167e+10,
		1.14498237732025810e+11, -3.99096175224466498e+11,
		8.19218669548577329e+11, -1.09837515608122331e+12,
		1.00815810686538209e+12, -6.45364869245376503e+11,
		2.87900649906150589e+11, -8.78670721780232657e+10,
		1.76347306068349694e+10, -2.16716498322379509e+09,
		1.43157876718888981e+08, -3.87183344257261262e+06,
		1.82577554742931747e+04, 2.86464035717679043e+11,
		-2.40629790002850396e+12, 9.10934118523989896e+12,
		-2.05168994109344374e+13, 3.05651255199353206e+13,
		-3.16670885847851584e+13, 2.33483640445818409e+13,
		-1.23204913055982872e+13, 4.61272578084913197e+12,
		-1.19655288019618160e+12, 2.05914503232410016e+11,
		-2.18229277575292237e+10, 1.24700929351271032e+09,
		-2.91883881222208134e+07, 1.18838426256783253e+05}

	if INIT != 0 {
		goto L40
	}

	RFN = 1.0e0 / FNU

	// OVERFLOW TEST (ZR/FNU TOO SMALL)
	TEST = D1MACH[1] * 1.0e+3
	AC = FNU * TEST
	if math.Abs(ZRR) > AC || math.Abs(ZRI) > AC {
		goto L15
	}
	ZETA1R = 2.0e0*math.Abs(math.Log(TEST)) + FNU
	ZETA1I = 0.0e0
	ZETA2R = FNU
	ZETA2I = 0.0e0
	PHIR = 1.0e0
	PHII = 0.0e0
	return ZRR, ZRI, FNU, IKFLG, IPMTR, TOL, INIT, PHIR, PHII, ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI
L15:
	TR = ZRR * RFN
	TI = ZRI * RFN
	SR = CONER + (TR*TR - TI*TI)
	SI = CONEI + (TR*TI + TI*TR)
	SRR, SRI = ZSQRT(SR, SI)
	STR = CONER + SRR
	STI = CONEI + SRI
	ZNR, ZNI = ZDIV(STR, STI, TR, TI)
	STR, STI = ZLOG(ZNR, ZNI)
	ZETA1R = FNU * STR
	ZETA1I = FNU * STI
	ZETA2R = FNU * SRR
	ZETA2I = FNU * SRI
	TR, TI = ZDIV(CONER, CONEI, SRR, SRI)
	SRR = TR * RFN
	SRI = TI * RFN
	CWRKR[16], CWRKI[16] = ZSQRT(SRR, SRI)
	PHIR = CWRKR[16] * CON[IKFLG]
	PHII = CWRKI[16] * CON[IKFLG]
	if IPMTR != 0 {
		return ZRR, ZRI, FNU, IKFLG, IPMTR, TOL, INIT, PHIR, PHII, ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI
	}
	T2R, T2I = ZDIV(CONER, CONEI, SR, SI)
	CWRKR[1] = CONER
	CWRKI[1] = CONEI
	CRFNR = CONER
	CRFNI = CONEI
	AC = 1.0e0
	L = 1
	for K = 2; K <= 15; K++ {
		SR = ZEROR
		SI = ZEROI
		for J = 1; J <= K; J++ {
			L = L + 1
			STR = SR*T2R - SI*T2I + C[L]
			SI = SR*T2I + SI*T2R
			SR = STR
		}
		STR = CRFNR*SRR - CRFNI*SRI
		CRFNI = CRFNR*SRI + CRFNI*SRR
		CRFNR = STR
		CWRKR[K] = CRFNR*SR - CRFNI*SI
		CWRKI[K] = CRFNR*SI + CRFNI*SR
		AC = AC * RFN
		TEST = math.Abs(CWRKR[K]) + math.Abs(CWRKI[K])
		if AC < TOL && TEST < TOL {
			goto L30
		}
	}
	K = 15
L30:
	INIT = K
L40:
	if IKFLG == 2 {
		goto L60
	}
	// COMPUTE SUM FOR THE I FUNCTION
	SR = ZEROR
	SI = ZEROI
	for I = 1; I < INIT; I++ {
		SR = SR + CWRKR[I]
		SI = SI + CWRKI[I]
	}
	SUMR = SR
	SUMI = SI
	PHIR = CWRKR[16] * CON[1]
	PHII = CWRKI[16] * CON[1]
	return ZRR, ZRI, FNU, IKFLG, IPMTR, TOL, INIT, PHIR, PHII, ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI
L60:
	// COMPUTE SUM FOR THE K FUNCTION
	SR = ZEROR
	SI = ZEROI
	TR = CONER
	for I = 1; I < INIT; I++ {
		SR = SR + TR*CWRKR[I]
		SI = SI + TR*CWRKI[I]
		TR = -TR
	}
	SUMR = SR
	SUMI = SI
	PHIR = CWRKR[16] * CON[2]
	PHII = CWRKI[16] * CON[2]
	return ZRR, ZRI, FNU, IKFLG, IPMTR, TOL, INIT, PHIR, PHII, ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI
}

// ZBUNI COMPUTES THE I BESSEL FUNCTION FOR LARGE CABS(Z) >  FNUL AND FNU+N-1 < FNUL.
// THE ORDER IS INCREASED FROM FNU+N-1 GREATER THAN FNUL BY ADDING NUI AND COMPUTING
// ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR I(FNU,Z) ON IFORM=1 AND
// THE EXPANSION FOR J(FNU,Z) ON IFORM=2
func ZBUNI(ZR float64, ZI float64, FNU float64, KODE int, N int, YR []float64, YI []float64, NZ int, NUI int, NLAST int,
	FNUL float64, TOL float64, ELIM float64, ALIM float64) (float64, float64, float64, int, int, []float64, []float64, int, int, int, float64, float64, float64, float64) {

	var AX, AY, CSCLR, CSCRR, DFNU, FNUI, GNU, RAZ, RZI, RZR, STI, STR, S1I, S1R, S2I, S2R, ASCLE, C1R, C1I, C1M float64

	var I, IFLAG, IFORM, K, NL, NW int
	var CYR, CYI []float64
	var BRY [4]float64

	NZ = 0
	AX = math.Abs(ZR) * 1.7321e0
	AY = math.Abs(ZI)
	IFORM = 1
	if AY > AX {
		IFORM = 2
	}
	if NUI == 0 {
		goto L60
	}
	FNUI = float64(float32(NUI))
	DFNU = FNU + float64(float32(N-1))
	GNU = DFNU + FNUI
	if IFORM == 2 {
		goto L10
	}

	// ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN -PI/3 <= ARG(Z) <= PI/3
	ZR, ZI, GNU, KODE, _, CYR, CYI, NW, NLAST, FNUL, TOL, ELIM, ALIM = ZUNI1(ZR, ZI, GNU, KODE, 2, CYR, CYI, NW, NLAST, FNUL, TOL, ELIM, ALIM)
	goto L20
L10:
	// ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU APPLIED IN
	// PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I AND HPI=PI/2
	ZR, ZI, GNU, KODE, _, CYR, CYI, NW, NLAST, FNUL, TOL, ELIM, ALIM = ZUNI2(ZR, ZI, GNU, KODE, 2, CYR, CYI, NW, NLAST, FNUL, TOL, ELIM, ALIM)
L20:
	if NW < 0 {
		goto L50
	}
	if NW != 0 {
		goto L90
	}
	STR = ZABS(CYR[1], CYI[1])

	// SCALE BACKWARD RECURRENCE, BRY(3) IS DEFINED BUT NEVER USED
	BRY[1] = 1000.0e0 * D1MACH[1] / TOL
	BRY[2] = 1.0e0 / BRY[1]
	BRY[3] = BRY[2]
	IFLAG = 2
	ASCLE = BRY[2]
	CSCLR = 1.0e0
	if STR > BRY[1] {
		goto L21
	}
	IFLAG = 1
	ASCLE = BRY[1]
	CSCLR = 1.0e0 / TOL
	goto L25
L21:
	if STR < BRY[2] {
		goto L25
	}
	IFLAG = 3
	ASCLE = BRY[3]
	CSCLR = TOL
L25:
	CSCRR = 1.0e0 / CSCLR
	S1R = CYR[2] * CSCLR
	S1I = CYI[2] * CSCLR
	S2R = CYR[1] * CSCLR
	S2I = CYI[1] * CSCLR
	RAZ = 1.0e0 / ZABS(ZR, ZI)
	STR = ZR * RAZ
	STI = -ZI * RAZ
	RZR = (STR + STR) * RAZ
	RZI = (STI + STI) * RAZ
	for I = 1; I <= NUI; I++ {
		STR = S2R
		STI = S2I
		S2R = (DFNU+FNUI)*(RZR*STR-RZI*STI) + S1R
		S2I = (DFNU+FNUI)*(RZR*STI+RZI*STR) + S1I
		S1R = STR
		S1I = STI
		FNUI = FNUI - 1.0e0
		if IFLAG >= 3 {
			continue
		}
		STR = S2R * CSCRR
		STI = S2I * CSCRR
		C1R = math.Abs(STR)
		C1I = math.Abs(STI)
		C1M = math.Max(C1R, C1I)
		if C1M <= ASCLE {
			continue
		}
		IFLAG = IFLAG + 1
		ASCLE = BRY[IFLAG]
		S1R = S1R * CSCRR
		S1I = S1I * CSCRR
		S2R = STR
		S2I = STI
		CSCLR = CSCLR * TOL
		CSCRR = 1.0e0 / CSCLR
		S1R = S1R * CSCLR
		S1I = S1I * CSCLR
		S2R = S2R * CSCLR
		S2I = S2I * CSCLR
	}
	YR[N] = S2R * CSCRR
	YI[N] = S2I * CSCRR
	if N == 1 {
		return ZR, ZI, FNU, KODE, N, YR, YI, NZ, NUI, NLAST, FNUL, TOL, ELIM, ALIM
	}
	NL = N - 1
	FNUI = float64(float32(NL))
	K = NL
	for I = 1; I <= NL; I++ {

		STR = S2R
		STI = S2I
		S2R = (FNU+FNUI)*(RZR*STR-RZI*STI) + S1R
		S2I = (FNU+FNUI)*(RZR*STI+RZI*STR) + S1I
		S1R = STR
		S1I = STI
		STR = S2R * CSCRR
		STI = S2I * CSCRR
		YR[K] = STR
		YI[K] = STI
		FNUI = FNUI - 1.0e0
		K = K - 1
		if IFLAG >= 3 {
			continue
		}
		C1R = math.Abs(STR)
		C1I = math.Abs(STI)
		C1M = math.Max(C1R, C1I)
		if C1M <= ASCLE {
			continue
		}
		IFLAG = IFLAG + 1
		ASCLE = BRY[IFLAG]
		S1R = S1R * CSCRR
		S1I = S1I * CSCRR
		S2R = STR
		S2I = STI
		CSCLR = CSCLR * TOL
		CSCRR = 1.0e0 / CSCLR
		S1R = S1R * CSCLR
		S1I = S1I * CSCLR
		S2R = S2R * CSCLR
		S2I = S2I * CSCLR
	}
	return ZR, ZI, FNU, KODE, N, YR, YI, NZ, NUI, NLAST, FNUL, TOL, ELIM, ALIM
L50:
	NZ = -1
	if NW == -2 {
		NZ = -2
	}
	return ZR, ZI, FNU, KODE, N, YR, YI, NZ, NUI, NLAST, FNUL, TOL, ELIM, ALIM
L60:
	if IFORM == 2 {
		goto L70
	}
	// ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN -PI/3 <= ARG(Z) <= PI/3
	ZR, ZI, FNU, KODE, N, YR, YI, NW, NLAST, FNUL, TOL, ELIM, ALIM = ZUNI1(ZR, ZI, FNU, KODE, N, YR, YI, NW, NLAST, FNUL, TOL, ELIM, ALIM)
	goto L80
L70:
	// ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU APPLIED IN PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I AND HPI=PI/2
	ZR, ZI, FNU, KODE, N, YR, YI, NW, NLAST, FNUL, TOL, ELIM, ALIM = ZUNI2(ZR, ZI, FNU, KODE, N, YR, YI, NW, NLAST, FNUL, TOL, ELIM, ALIM)
L80:
	if NW < 0 {
		goto L50
	}
	NZ = NW
	return ZR, ZI, FNU, KODE, N, YR, YI, NZ, NUI, NLAST, FNUL, TOL, ELIM, ALIM
L90:
	NLAST = N
	return ZR, ZI, FNU, KODE, N, YR, YI, NZ, NUI, NLAST, FNUL, TOL, ELIM, ALIM
}

// ZUNI1 COMPUTES I(FNU,Z)  BY MEANS OF THE UNIFORM ASYMPTOTIC
// EXPANSION FOR I(FNU,Z) IN -PI/3 <= ARG Z <= PI/3.
//
// FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
// EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
// NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
// FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1 < FNUL.
// Y(I)=CZERO FOR I=NLAST+1,N
func ZUNI1(ZR float64, ZI float64, FNU float64, KODE int, N int, YR []float64, YI []float64, NZ int, NLAST int, FNUL float64, TOL float64, ELIM float64, ALIM float64) (float64, float64, float64, int, int, []float64, []float64, int, int, float64, float64, float64, float64) {

	const (
		ZEROR = 0.0e0
		ZEROI = 0.0e0
		CONER = 1.0e0
	)
	var APHI, ASCLE, CRSC, CSCL, C1R, C2I, C2M, C2R, FN,
		PHII, PHIR, RAST, RS1, RZI, RZR, STI, STR, SUMI,
		SUMR, S1I, S1R, S2I, S2R, ZETA1I,
		ZETA1R, ZETA2I, ZETA2R float64
	var I, IFLAG, INIT, K, M, ND, NN, NUF, NW int

	var CYR, CYI [3]float64
	var CSRR, CSSR, BRY [4]float64
	var CWRKR, CWRKI []float64

	NZ = 0
	ND = N
	NLAST = 0
	// COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
	// NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
	// EXP(ALIM)=EXP(ELIM)*TOL
	CSCL = 1.0e0 / TOL
	CRSC = TOL
	CSSR[1] = CSCL
	CSSR[2] = CONER
	CSSR[3] = CRSC
	CSRR[1] = CRSC
	CSRR[2] = CONER
	CSRR[3] = CSCL
	BRY[1] = 1000.0e0 * D1MACH[1] / TOL
	// CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
	FN = math.Max(FNU, 1.0e0)
	INIT = 0
	ZR, ZI, FN, _, _, TOL, INIT, PHIR, PHII, ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI = ZUNIK(ZR, ZI, FN, 1, 1, TOL, INIT, PHIR, PHII, ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
	if KODE == 1 {
		goto L10
	}
	STR = ZR + ZETA2R
	STI = ZI + ZETA2I
	RAST = FN / ZABS(STR, STI)
	STR = STR * RAST * RAST
	STI = -STI * RAST * RAST
	S1R = -ZETA1R + STR
	S1I = -ZETA1I + STI
	goto L20
L10:
	S1R = -ZETA1R + ZETA2R
	S1I = -ZETA1I + ZETA2I
L20:
	RS1 = S1R
	if math.Abs(RS1) > ELIM {
		goto L130
	}
L30:
	NN = MIN(2, ND)
	for I = 1; I <= NN; I++ {
		FN = FNU + float64(float32(ND-I))
		INIT = 0
		ZR, ZI, FN, _, _, TOL, INIT, PHIR, PHII, ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI = ZUNIK(ZR, ZI, FN, 1, 0, TOL, INIT, PHIR, PHII, ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
		if KODE == 1 {
			continue
		}
		STR = ZR + ZETA2R
		STI = ZI + ZETA2I
		RAST = FN / ZABS(STR, STI)
		STR = STR * RAST * RAST
		STI = -STI * RAST * RAST
		S1R = -ZETA1R + STR
		S1I = -ZETA1I + STI + ZI
		goto L50

		S1R = -ZETA1R + ZETA2R
		S1I = -ZETA1I + ZETA2I
	L50:
		// TEST FOR UNDERFLOW AND OVERFLOW
		RS1 = S1R
		if math.Abs(RS1) > ELIM {
			goto L110
		}
		if I == 1 {
			IFLAG = 2
		}
		if math.Abs(RS1) < ALIM {
			goto L60
		}
		// REFINE  TEST AND SCALE
		APHI = ZABS(PHIR, PHII)
		RS1 = RS1 + math.Log(APHI)
		if math.Abs(RS1) > ELIM {
			goto L110
		}
		if I == 1 {
			IFLAG = 1
		}
		if RS1 < 0.0e0 {
			goto L60
		}
		if I == 1 {
			IFLAG = 3
		}
	L60:
		// SCALE S1 IF CABS(S1) < ASCLE
		S2R = PHIR*SUMR - PHII*SUMI
		S2I = PHIR*SUMI + PHII*SUMR
		STR = math.Exp(S1R) * CSSR[IFLAG]
		S1R = STR * math.Cos(S1I)
		S1I = STR * math.Sin(S1I)
		STR = S2R*S1R - S2I*S1I
		S2I = S2R*S1I + S2I*S1R
		S2R = STR
		if IFLAG != 1 {
			goto L70
		}
		S2R, S2I, NW, _, TOL = ZUCHK(S2R, S2I, NW, BRY[1], TOL)
		if NW != 0 {
			goto L110
		}
	L70:
		CYR[I] = S2R
		CYI[I] = S2I
		M = ND - I + 1
		YR[M] = S2R * CSRR[IFLAG]
		YI[M] = S2I * CSRR[IFLAG]
	}
	if ND <= 2 {
		goto L100
	}
	RAST = 1.0e0 / ZABS(ZR, ZI)
	STR = ZR * RAST
	STI = -ZI * RAST
	RZR = (STR + STR) * RAST
	RZI = (STI + STI) * RAST
	BRY[2] = 1.0e0 / BRY[1]
	BRY[3] = D1MACH[2]
	S1R = CYR[1]
	S1I = CYI[1]
	S2R = CYR[2]
	S2I = CYI[2]
	C1R = CSRR[IFLAG]
	ASCLE = BRY[IFLAG]
	K = ND - 2
	FN = float64(float32(K))
	for I = 3; I <= ND; I++ {
		C2R = S2R
		C2I = S2I
		S2R = S1R + (FNU+FN)*(RZR*C2R-RZI*C2I)
		S2I = S1I + (FNU+FN)*(RZR*C2I+RZI*C2R)
		S1R = C2R
		S1I = C2I
		C2R = S2R * C1R
		C2I = S2I * C1R
		YR[K] = C2R
		YI[K] = C2I
		K = K - 1
		FN = FN - 1.0e0
		if IFLAG >= 3 {
			continue
		}
		STR = math.Abs(C2R)
		STI = math.Abs(C2I)
		C2M = math.Max(STR, STI)
		if C2M <= ASCLE {
			continue
		}
		IFLAG = IFLAG + 1
		ASCLE = BRY[IFLAG]
		S1R = S1R * C1R
		S1I = S1I * C1R
		S2R = C2R
		S2I = C2I
		S1R = S1R * CSSR[IFLAG]
		S1I = S1I * CSSR[IFLAG]
		S2R = S2R * CSSR[IFLAG]
		S2I = S2I * CSSR[IFLAG]
		C1R = CSRR[IFLAG]
	}
L100:
	return ZR, ZI, FNU, KODE, N, YR, YI, NZ, NLAST, FNUL, TOL, ELIM, ALIM
	// SET UNDERFLOW AND UPDATE PARAMETERS
L110:
	if RS1 > 0.0e0 {
		goto L120
	}
	YR[ND] = ZEROR
	YI[ND] = ZEROI
	NZ = NZ + 1
	ND = ND - 1
	if ND == 0 {
		goto L100
	}
	ZR, ZI, FNU, KODE, _, ND, YR, YI, NUF, TOL, ELIM, ALIM = ZUOIK(ZR, ZI, FNU, KODE, 1, ND, YR, YI, NUF, TOL, ELIM, ALIM)
	if NUF < 0 {
		goto L120
	}
	ND = ND - NUF
	NZ = NZ + NUF
	if ND == 0 {
		goto L100
	}
	FN = FNU + float64(float32(ND-1))
	if FN >= FNUL {
		goto L30
	}
	NLAST = ND
	return ZR, ZI, FNU, KODE, N, YR, YI, NZ, NLAST, FNUL, TOL, ELIM, ALIM
L120:
	NZ = -1
	return ZR, ZI, FNU, KODE, N, YR, YI, NZ, NLAST, FNUL, TOL, ELIM, ALIM
L130:
	if RS1 > 0.0e0 {
		goto L120
	}
	NZ = N
	for I = 1; I <= N; I++ {
		YR[I] = ZEROR
		YI[I] = ZEROI
	}
	return ZR, ZI, FNU, KODE, N, YR, YI, NZ, NLAST, FNUL, TOL, ELIM, ALIM
}

// ZBINU COMPUTES THE I FUNCTION IN THE RIGHT HALF Z PLANE
func ZBINU(ZR float64, ZI float64, FNU float64, KODE int, N int, CYR []float64, CYI []float64, NZ int,
	RL float64, FNUL float64, TOL float64, ELIM float64, ALIM float64) (float64, float64, float64, int, int, []float64, []float64, int, float64, float64, float64, float64, float64) {

	const (
		ZEROR = 0.0e0
		ZEROI = 0.0e0
	)

	var AZ, DFNU float64
	var I, INW, NLAST, NN, NUI, NW int

	CWR := []float64{math.NaN(), 0, 0}
	CWI := []float64{math.NaN(), 0, 0}

	NZ = 0
	AZ = ZABS(ZR, ZI)
	NN = N
	DFNU = FNU + float64(float32(N-1))
	if AZ <= 2.0e0 {
		goto L10
	}
	if AZ*AZ*0.25e0 > DFNU+1.0e0 {
		goto L20
	}
L10:
	// POWER SERIES
	ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, TOL, ELIM, ALIM = ZSERI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, TOL, ELIM, ALIM)
	INW = ABS(NW)
	NZ = NZ + INW
	NN = NN - INW
	if NN == 0 {
		return ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, RL, FNUL, TOL, ELIM, ALIM
	}
	if NW >= 0 {
		goto L120
	}
	DFNU = FNU + float64(float32(NN-1))
L20:
	if AZ < RL {
		goto L40
	}
	if DFNU <= 1.0e0 {
		goto L30
	}
	if AZ+AZ < DFNU*DFNU {
		goto L50
	}

	// ASYMPTOTIC EXPANSION FOR LARGE Z
L30:
	ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, RL, TOL, ELIM, ALIM = ZASYI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, RL, TOL, ELIM, ALIM)
	if NW < 0 {
		goto L130
	}
	goto L120
L40:
	if DFNU <= 1.0e0 {
		goto L70
	}
L50:
	// OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM
	ZR, ZI, FNU, KODE, _, NN, CYR, CYI, NW, TOL, ELIM, ALIM = ZUOIK(ZR, ZI, FNU, KODE, 1, NN, CYR, CYI, NW, TOL, ELIM, ALIM)
	if NW < 0 {
		goto L130
	}
	NZ = NZ + NW
	NN = NN - NW
	if NN == 0 {
		return ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, RL, FNUL, TOL, ELIM, ALIM
	}
	DFNU = FNU + float64(float32(NN-1))
	if DFNU > FNUL {
		goto L110
	}
	if AZ > FNUL {
		goto L110
	}
L60:
	if AZ > RL {
		goto L80
	}
L70:
	// MILLER ALGORITHM NORMALIZED BY THE SERIES
	ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, TOL = ZMLRI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, TOL)
	if NW < 0 {
		goto L130
	}
	goto L120
L80:
	// MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN
	// OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN
	ZR, ZI, FNU, KODE, _, _, CWR, CWI, NW, TOL, ELIM, ALIM = ZUOIK(ZR, ZI, FNU, KODE, 2, 2, CWR, CWI, NW, TOL, ELIM, ALIM)
	if NW >= 0 {
		goto L100
	}
	NZ = NN
	for I = 1; I <= NN; I++ {
		CYR[I] = ZEROR
		CYI[I] = ZEROI
	}
	return ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, RL, FNUL, TOL, ELIM, ALIM
L100:
	if NW > 0 {
		goto L130
	}
	ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, CWR, CWI, TOL, ELIM, ALIM = ZWRSK(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, CWR, CWI, TOL, ELIM, ALIM)
	if NW < 0 {
		goto L130
	}
	goto L120
L110:
	// INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD
	NUI = int(float32(FNUL-DFNU)) + 1
	NUI = MAX(NUI, 0)
	ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, NUI, NLAST, FNUL, TOL, ELIM, ALIM = ZBUNI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, NUI, NLAST, FNUL, TOL, ELIM, ALIM)
	if NW < 0 {
		goto L130
	}
	NZ = NZ + NW
	if NLAST == 0 {
		goto L120
	}
	NN = NLAST
	goto L60
L120:
	return ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, RL, FNUL, TOL, ELIM, ALIM

L130:
	NZ = -1
	if NW == -2 {
		NZ = -2
	}
	return ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, RL, FNUL, TOL, ELIM, ALIM
}

// ZWRSK COMPUTES THE I BESSEL FUNCTION FOR RE(Z) >= 0.0 BY
// NORMALIZING THE I FUNCTION RATIOS FROM ZRATI BY THE WRONSKIAN
func ZWRSK(ZRR float64, ZRI float64, FNU float64, KODE int, N int, YR []float64, YI []float64, NZ int, CWR []float64, CWI []float64, TOL float64, ELIM float64, ALIM float64) (float64, float64, float64, int, int, []float64, []float64, int, []float64, []float64, float64, float64, float64) {

	var ACT, ACW, ASCLE, CINUI, CINUR, CSCLR, CTI, CTR, C1R, C1I, C2I, C2R, PTI, PTR, RACT, STI, STR float64
	var I, NW int

	// I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS
	// Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM CRATI NORMALIZED BY THE
	// WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM CBKNU.
	NZ = 0
	ZRR, ZRI, FNU, KODE, _, CWR, CWI, NW, TOL, ELIM, ALIM = ZBKNU(ZRR, ZRI, FNU, KODE, 2, CWR, CWI, NW, TOL, ELIM, ALIM)
	if NW != 0 {
		goto L50
	}
	ZRR, ZRI, FNU, N, YR, YI, TOL = ZRATI(ZRR, ZRI, FNU, N, YR, YI, TOL)

	// RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z),
	// R(FNU+J-1,Z)=Y(J),  J=1,...,N
	CINUR = 1.0e0
	CINUI = 0.0e0
	if KODE == 1 {
		goto L10
	}
	CINUR = math.Cos(ZRI)
	CINUI = math.Sin(ZRI)
L10:
	// ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH
	// THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE
	// SCALED TO PREVENT OVER OR UNDERFLOW. CUOIK HAS DETERMINED THAT
	// THE RESULT IS ON SCALE.
	ACW = ZABS(CWR[2], CWI[2])
	ASCLE = 1000.0e0 * D1MACH[1] / TOL
	CSCLR = 1.0e0
	if ACW > ASCLE {
		goto L20
	}
	CSCLR = 1.0e0 / TOL
	goto L30
L20:
	ASCLE = 1.0e0 / ASCLE
	if ACW < ASCLE {
		goto L30
	}
	CSCLR = TOL
L30:
	C1R = CWR[1] * CSCLR
	C1I = CWI[1] * CSCLR
	C2R = CWR[2] * CSCLR
	C2I = CWI[2] * CSCLR
	STR = YR[1]
	STI = YI[1]

	// CINU=CINU*(CONJG(CT)/CABS(CT))*(1.0e0/CABS(CT) PREVENTS
	// UNDER- OR OVERFLOW PREMATURELY BY SQUARING CABS(CT)
	PTR = STR*C1R - STI*C1I
	PTI = STR*C1I + STI*C1R
	PTR = PTR + C2R
	PTI = PTI + C2I
	CTR = ZRR*PTR - ZRI*PTI
	CTI = ZRR*PTI + ZRI*PTR
	ACT = ZABS(CTR, CTI)
	RACT = 1.0e0 / ACT
	CTR = CTR * RACT
	CTI = -CTI * RACT
	PTR = CINUR * RACT
	PTI = CINUI * RACT
	CINUR = PTR*CTR - PTI*CTI
	CINUI = PTR*CTI + PTI*CTR
	YR[1] = CINUR * CSCLR
	YI[1] = CINUI * CSCLR
	if N == 1 {
		return ZRR, ZRI, FNU, KODE, N, YR, YI, NZ, CWR, CWI, TOL, ELIM, ALIM
	}
	for I = 2; I <= N; I++ {

		PTR = STR*CINUR - STI*CINUI
		CINUI = STR*CINUI + STI*CINUR
		CINUR = PTR
		STR = YR[I]
		STI = YI[I]
		YR[I] = CINUR * CSCLR
		YI[I] = CINUI * CSCLR
	}
	return ZRR, ZRI, FNU, KODE, N, YR, YI, NZ, CWR, CWI, TOL, ELIM, ALIM
L50:
	NZ = -1
	if NW == -2 {
		NZ = -2
	}
	return ZRR, ZRI, FNU, KODE, N, YR, YI, NZ, CWR, CWI, TOL, ELIM, ALIM
}

// ZRATI COMPUTES RATIOS OF I BESSEL FUNCTIONS BY BACKWARD
// RECURRENCE.  THE STARTING INDEX IS DETERMINED BY FORWARD
// RECURRENCE AS DESCRIBED IN J. RES. OF NAT. BUR. OF STANDARDS-B,
// MATHEMATICAL SCIENCES, VOL 77B, P111-114, SEPTEMBER, 1973,
// BESSEL FUNCTIONS I AND J OF COMPLEX ARGUMENT AND INTEGER ORDER,
//   BY D. J. SOOKNE.
func ZRATI(ZR float64, ZI float64, FNU float64, N int, CYR []float64, CYI []float64, TOL float64) (float64, float64, float64, int, []float64, []float64, float64) {

	const (
		CZEROR = 0.0e0
		CZEROI = 0.0e0
		CONER  = 1.0e0
		CONEI  = 0.0e0
		RT2    = 1.41421356237309505e0
	)

	var AK, AMAGZ, AP1, AP2, ARG, AZ, CDFNUI, CDFNUR, DFNU, FDNU, FLAM, FNUP, PTI, PTR, P1I, P1R,
		P2I, P2R, RAK, RAP1, RHO, RZI, RZR, TEST, TEST1, TTI, TTR, T1I, T1R float64
	var I, ID, IDNU, INU, ITIME, K, KK, MAGZ int

	AZ = ZABS(ZR, ZI)
	INU = int(float32(FNU))
	IDNU = INU + N - 1
	MAGZ = int(float32(AZ))
	AMAGZ = float64(float32(MAGZ + 1))
	FDNU = float64(float32(IDNU))
	FNUP = math.Max(AMAGZ, FDNU)
	ID = IDNU - MAGZ - 1
	ITIME = 1
	K = 1
	PTR = 1.0e0 / AZ
	RZR = PTR * (ZR + ZR) * PTR
	RZI = -PTR * (ZI + ZI) * PTR
	T1R = RZR * FNUP
	T1I = RZI * FNUP
	P2R = -T1R
	P2I = -T1I
	P1R = CONER
	P1I = CONEI
	T1R = T1R + RZR
	T1I = T1I + RZI
	if ID > 0 {
		ID = 0
	}
	AP2 = ZABS(P2R, P2I)
	AP1 = ZABS(P1R, P1I)
	// THE OVERFLOW TEST ON K(FNU+I-1,Z) BEFORE THE CALL TO CBKNU
	// GUARANTEES THAT P2 IS ON SCALE. SCALE TEST1 AND ALL SUBSEQUENT
	// P2 VALUES BY AP1 TO ENSURE THAT AN OVERFLOW DOES NOT OCCUR
	// PREMATURELY.
	ARG = (AP2 + AP2) / (AP1 * TOL)
	TEST1 = math.Sqrt(ARG)
	TEST = TEST1
	RAP1 = 1.0e0 / AP1
	P1R = P1R * RAP1
	P1I = P1I * RAP1
	P2R = P2R * RAP1
	P2I = P2I * RAP1
	AP2 = AP2 * RAP1
L10:
	K = K + 1
	AP1 = AP2
	PTR = P2R
	PTI = P2I
	P2R = P1R - (T1R*PTR - T1I*PTI)
	P2I = P1I - (T1R*PTI + T1I*PTR)
	P1R = PTR
	P1I = PTI
	T1R = T1R + RZR
	T1I = T1I + RZI
	AP2 = ZABS(P2R, P2I)
	if AP1 <= TEST {
		goto L10
	}
	if ITIME == 2 {
		goto L20
	}
	AK = ZABS(T1R, T1I) * 0.5e0
	FLAM = AK + math.Sqrt(AK*AK-1.0e0)
	RHO = math.Min(AP2/AP1, FLAM)
	TEST = TEST1 * math.Sqrt(RHO/(RHO*RHO-1.0e0))
	ITIME = 2
	goto L10
L20:
	KK = K + 1 - ID
	AK = float64(float32(KK))
	T1R = AK
	T1I = CZEROI
	DFNU = FNU + float64(float32(N-1))
	P1R = 1.0e0 / AP2
	P1I = CZEROI
	P2R = CZEROR
	P2I = CZEROI
	for I = 1; I <= KK; I++ {
		PTR = P1R
		PTI = P1I
		RAP1 = DFNU + T1R
		TTR = RZR * RAP1
		TTI = RZI * RAP1
		P1R = (PTR*TTR - PTI*TTI) + P2R
		P1I = (PTR*TTI + PTI*TTR) + P2I
		P2R = PTR
		P2I = PTI
		T1R = T1R - CONER
	}
	if P1R != CZEROR || P1I != CZEROI {
		goto L40
	}
	P1R = TOL
	P1I = TOL
L40:
	CYR[N], CYI[N] = ZDIV(P2R, P2I, P1R, P1I)
	if N == 1 {
		return ZR, ZI, FNU, N, CYR, CYI, TOL
	}
	K = N - 1
	AK = float64(float32(K))
	T1R = AK
	T1I = CZEROI
	CDFNUR = FNU * RZR
	CDFNUI = FNU * RZI
	for I = 2; I <= N; I++ {
		PTR = CDFNUR + (T1R*RZR - T1I*RZI) + CYR[K+1]
		PTI = CDFNUI + (T1R*RZI + T1I*RZR) + CYI[K+1]
		AK = ZABS(PTR, PTI)
		if AK != CZEROR {
			goto L50
		}
		PTR = TOL
		PTI = TOL
		AK = TOL * RT2
	L50:
		RAK = CONER / AK
		CYR[K] = RAK * PTR * RAK
		CYI[K] = -RAK * PTI * RAK
		T1R = T1R - CONER
		K = K - 1
	}
	return ZR, ZI, FNU, N, CYR, CYI, TOL
}

// ZUNI2 COMPUTES I(FNU,Z) IN THE RIGHT HALF PLANE BY MEANS OF
// UNIFORM ASYMPTOTIC EXPANSION FOR J(FNU,ZN) WHERE ZN IS Z*I
// OR -Z*I AND ZN IS IN THE RIGHT HALF PLANE ALSO.
//
// FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
// EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
// NLAST != 0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
// FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1 < FNUL.
// Y(I)=CZERO FOR I=NLAST+1,N
func ZUNI2(ZR float64, ZI float64, FNU float64, KODE int, N int, YR []float64, YI []float64, NZ int, NLAST int, FNUL float64, TOL float64, ELIM float64, ALIM float64) (float64, float64, float64, int, int, []float64, []float64, int, int, float64, float64, float64, float64) {

	const (
		ZEROR = 0.0e0
		ZEROI = 0.0e0
		CONER = 1.0e0
		RT2   = 1.41421356237309505e0
		HPI   = 1.57079632679489662e+00
		AIC   = 1.265512123484645396e+00
	)

	var AARG, AII, AIR, ANG, APHI, ARGI, ARGR, ASCLE, ASUMI, ASUMR, BSUMI, BSUMR, CIDI,
		CRSC, CSCL, C1R, C2I, C2M, C2R, DAII, DAIR, FN, PHII, PHIR, RAST, RAZ, RS1, RZI, RZR,
		STI, STR, S1I, S1R, S2I, S2R, ZBI, ZBR, ZETA1I, ZETA1R, ZETA2I, ZETA2R, ZNI, ZNR, CAR, SAR float64

	var I, IFLAG, IN, INU, J, K, ND, NN, NUF, NW int

	var BRY, CSSR, CSRR [4]float64
	var CYR, CYI [4]float64
	var CIPR = [5]float64{math.NaN(), 1.0e0, 0.0e0, -1.0e0, 0.0e0}
	var CIPI = [5]float64{math.NaN(), 0.0e0, 1.0e0, 0.0e0, -1.0e0}

	NZ = 0
	ND = N
	NLAST = 0

	// COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
	// NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
	// EXP(ALIM)=EXP(ELIM)*TOL
	CSCL = 1.0e0 / TOL
	CRSC = TOL
	CSSR[1] = CSCL
	CSSR[2] = CONER
	CSSR[3] = CRSC
	CSRR[1] = CRSC
	CSRR[2] = CONER
	CSRR[3] = CSCL
	BRY[1] = 1000.0e0 * D1MACH[1] / TOL

	// ZN IS IN THE RIGHT HALF PLANE AFTER ROTATION BY CI OR -CI
	ZNR = ZI
	ZNI = -ZR
	ZBR = ZR
	ZBI = ZI
	CIDI = -CONER
	INU = int(float32(FNU))
	ANG = HPI * (FNU - float64(float32(INU)))
	C2R = math.Cos(ANG)
	C2I = math.Sin(ANG)
	CAR = C2R
	SAR = C2I
	IN = INU + N - 1
	IN = IN%4 + 1
	STR = C2R*CIPR[IN] - C2I*CIPI[IN]
	C2I = C2R*CIPI[IN] + C2I*CIPR[IN]
	C2R = STR
	if ZI > 0.0e0 {
		goto L10
	}
	ZNR = -ZNR
	ZBI = -ZBI
	CIDI = -CIDI
	C2I = -C2I
L10:
	// CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
	FN = math.Max(FNU, 1.0e0)
	ZNR, ZNI, FN, _, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R, ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI = ZUNHJ(ZNR, ZNI, FN, 1, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R, ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
	if KODE == 1 {
		goto L20
	}
	STR = ZBR + ZETA2R
	STI = ZBI + ZETA2I
	RAST = FN / ZABS(STR, STI)
	STR = STR * RAST * RAST
	STI = -STI * RAST * RAST
	S1R = -ZETA1R + STR
	S1I = -ZETA1I + STI
	goto L30
L20:
	S1R = -ZETA1R + ZETA2R
	S1I = -ZETA1I + ZETA2I
L30:
	RS1 = S1R
	if math.Abs(RS1) > ELIM {
		goto L150
	}
L40:
	NN = MIN(2, ND)
	for I = 1; I <= NN; I++ {
		FN = FNU + float64(float32(ND-I))
		ZNR, ZNI, FN, _, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R, ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI = ZUNHJ(ZNR, ZNI, FN, 0, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R, ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
		if KODE == 1 {
			goto L50
		}
		STR = ZBR + ZETA2R
		STI = ZBI + ZETA2I
		RAST = FN / ZABS(STR, STI)
		STR = STR * RAST * RAST
		STI = -STI * RAST * RAST
		S1R = -ZETA1R + STR
		S1I = -ZETA1I + STI + math.Abs(ZI)
		goto L60
	L50:
		S1R = -ZETA1R + ZETA2R
		S1I = -ZETA1I + ZETA2I
	L60:
		// TEST FOR UNDERFLOW AND OVERFLOW
		RS1 = S1R
		if math.Abs(RS1) > ELIM {
			goto L120
		}
		if I == 1 {
			IFLAG = 2
		}
		if math.Abs(RS1) < ALIM {
			goto L70
		}

		// REFINE  TEST AND SCALE
		APHI = ZABS(PHIR, PHII)
		AARG = ZABS(ARGR, ARGI)
		RS1 = RS1 + math.Log(APHI) - 0.25e0*math.Log(AARG) - AIC
		if math.Abs(RS1) > ELIM {
			goto L120
		}
		if I == 1 {
			IFLAG = 1
		}
		if RS1 < 0.0e0 {
			goto L70
		}
		if I == 1 {
			IFLAG = 3
		}
	L70:

		// SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
		// EXPONENT EXTREMES
		AIR, AII, _, _ = ZAIRY(ARGR, ARGI, 0, 2)
		DAIR, DAII, _, _ = ZAIRY(ARGR, ARGI, 1, 2)

		STR = DAIR*BSUMR - DAII*BSUMI
		STI = DAIR*BSUMI + DAII*BSUMR
		STR = STR + (AIR*ASUMR - AII*ASUMI)
		STI = STI + (AIR*ASUMI + AII*ASUMR)
		S2R = PHIR*STR - PHII*STI
		S2I = PHIR*STI + PHII*STR
		STR = math.Exp(S1R) * CSSR[IFLAG]
		S1R = STR * math.Cos(S1I)
		S1I = STR * math.Sin(S1I)
		STR = S2R*S1R - S2I*S1I
		S2I = S2R*S1I + S2I*S1R
		S2R = STR
		if IFLAG != 1 {
			goto L80
		}
		S2R, S2I, NW, BRY[1], TOL = ZUCHK(S2R, S2I, NW, BRY[1], TOL)
		if NW != 0 {
			goto L120
		}
	L80:
		if ZI <= 0.0e0 {
			S2I = -S2I
		}
		STR = S2R*C2R - S2I*C2I
		S2I = S2R*C2I + S2I*C2R
		S2R = STR
		CYR[I] = S2R
		CYI[I] = S2I
		J = ND - I + 1
		YR[J] = S2R * CSRR[IFLAG]
		YI[J] = S2I * CSRR[IFLAG]
		STR = -C2I * CIDI
		C2I = C2R * CIDI
		C2R = STR
	}
	if ND <= 2 {
		goto L110
	}
	RAZ = 1.0e0 / ZABS(ZR, ZI)
	STR = ZR * RAZ
	STI = -ZI * RAZ
	RZR = (STR + STR) * RAZ
	RZI = (STI + STI) * RAZ
	BRY[2] = 1.0e0 / BRY[1]
	BRY[3] = D1MACH[2]
	S1R = CYR[1]
	S1I = CYI[1]
	S2R = CYR[2]
	S2I = CYI[2]
	C1R = CSRR[IFLAG]
	ASCLE = BRY[IFLAG]
	K = ND - 2
	FN = float64(float32(K))
	for I = 3; I <= ND; I++ {
		C2R = S2R
		C2I = S2I
		S2R = S1R + (FNU+FN)*(RZR*C2R-RZI*C2I)
		S2I = S1I + (FNU+FN)*(RZR*C2I+RZI*C2R)
		S1R = C2R
		S1I = C2I
		C2R = S2R * C1R
		C2I = S2I * C1R
		YR[K] = C2R
		YI[K] = C2I
		K = K - 1
		FN = FN - 1.0e0
		if IFLAG >= 3 {
			continue
		}
		STR = math.Abs(C2R)
		STI = math.Abs(C2I)
		C2M = math.Max(STR, STI)
		if C2M <= ASCLE {
			continue
		}
		IFLAG = IFLAG + 1
		ASCLE = BRY[IFLAG]
		S1R = S1R * C1R
		S1I = S1I * C1R
		S2R = C2R
		S2I = C2I
		S1R = S1R * CSSR[IFLAG]
		S1I = S1I * CSSR[IFLAG]
		S2R = S2R * CSSR[IFLAG]
		S2I = S2I * CSSR[IFLAG]
		C1R = CSRR[IFLAG]
	}
L110:
	return ZR, ZI, FNU, KODE, N, YR, YI, NZ, NLAST, FNUL, TOL, ELIM, ALIM
L120:
	if RS1 > 0.0e0 {
		goto L140
	}

	// SET UNDERFLOW AND UPDATE PARAMETERS
	YR[ND] = ZEROR
	YI[ND] = ZEROI
	NZ = NZ + 1
	ND = ND - 1
	if ND == 0 {
		goto L110
	}
	ZR, ZI, FNU, KODE, _, ND, YR, YI, NUF, TOL, ELIM, ALIM = ZUOIK(ZR, ZI, FNU, KODE, 1, ND, YR, YI, NUF, TOL, ELIM, ALIM)
	if NUF < 0 {
		goto L140
	}
	ND = ND - NUF
	NZ = NZ + NUF
	if ND == 0 {
		goto L110
	}
	FN = FNU + float64(float32(ND-1))
	if FN < FNUL {
		goto L130
	}
	IN = INU + ND - 1
	IN = IN%4 + 1
	C2R = CAR*CIPR[IN] - SAR*CIPI[IN]
	C2I = CAR*CIPI[IN] + SAR*CIPR[IN]
	if ZI <= 0.0e0 {
		C2I = -C2I
	}
	goto L40
L130:
	NLAST = ND
	return ZR, ZI, FNU, KODE, N, YR, YI, NZ, NLAST, FNUL, TOL, ELIM, ALIM
L140:
	NZ = -1
	return ZR, ZI, FNU, KODE, N, YR, YI, NZ, NLAST, FNUL, TOL, ELIM, ALIM
L150:
	if RS1 > 0.0e0 {
		goto L140
	}
	NZ = N
	for I = 1; I <= N; I++ {
		YR[I] = ZEROR
		YI[I] = ZEROI
	}
	return ZR, ZI, FNU, KODE, N, YR, YI, NZ, NLAST, FNUL, TOL, ELIM, ALIM

}

// ZUNHJ COMPUTES PARAMETERS FOR BESSEL FUNCTIONS C(FNU,Z) =
// J(FNU,Z), Y(FNU,Z) OR H(I,FNU,Z) I=1,2 FOR LARGE ORDERS FNU
// BY MEANS OF THE UNIFORM ASYMPTOTIC EXPANSION
//
// C(FNU,Z)=C1*PHI*( ASUM*AIRY(ARG) + C2*BSUM*DAIRY(ARG) )
//
// FOR PROPER CHOICES OF C1, C2, AIRY AND DAIRY WHERE AIRY IS
// AN AIRY FUNCTION AND DAIRY IS ITS DERIVATIVE.
//
//       (2/3)*FNU*ZETA**1.5 = ZETA1-ZETA2,
//
// ZETA1=0.5*FNU*CLOG((1+W)/(1-W)), ZETA2=FNU*W FOR SCALING
// PURPOSES IN AIRY FUNCTIONS FROM CAIRY OR CBIRY.
//
// MCONJ=SIGN OF AIMAG(Z), BUT IS AMBIGUOUS WHEN Z IS REAL AND
// MUST BE SPECIFIED. IPMTR=0 RETURNS ALL PARAMETERS. IPMTR=
// 1 COMPUTES ALL EXCEPT ASUM AND BSUM.
//
// HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ AND I.A.
// STEGUN, AMS55, NATIONAL BUREAU OF STANDARDS, 1965, CHAPTER 9.
// ASYMPTOTICS AND SPECIAL FUNCTIONS BY F.W.J. OLVER, ACADEMIC
// PRESS, N.Y., 1974, PAGE 420
func ZUNHJ(ZR float64, ZI float64, FNU float64, IPMTR int, TOL float64, PHIR float64, PHII float64,
	ARGR float64, ARGI float64, ZETA1R float64, ZETA1I float64, ZETA2R float64, ZETA2I float64,
	ASUMR float64, ASUMI float64, BSUMR float64, BSUMI float64) (float64, float64, float64, int, float64, float64, float64,
	float64, float64, float64, float64, float64, float64, float64, float64, float64, float64) {

	const (
		ZEROR = 0.0e0
		ZEROI = 0.0e0
		CONER = 1.0e0
		CONEI = 1.0e0
		EX1   = 3.33333333333333333e-01
		EX2   = 6.66666666666666667e-01
		HPI   = 1.57079632679489662e+00
		GPI   = 3.14159265358979324e+00
		THPI  = 4.71238898038468986e+00
	)

	var ANG, ATOL, AW2, AZTH, BTOL, FN13, FN23,
		PP, PRZTHI, PRZTHR, PTFNI, PTFNR, RAW, RAW2, RAZTH, RFNU, RFNU2, RFN13, RTZTI, RTZTR, RZTHI, RZTHR, STI, STR,
		SUMAI, SUMAR, SUMBI, SUMBR, TEST, TFNI, TFNR, TZAI, TZAR, T2I, T2R, WI, WR, W2I, W2R, ZAI, ZAR, ZBI, ZBR,
		ZCI, ZCR, ZETAI, ZETAR, ZTHI, ZTHR, AC float64
	var IAS, IBS, IS, J, JR, JU, K, KMAX, KP1, KS, L, LR, LRP1, L1, L2, M int

	var AP, PR, PI [31]float64
	var UPR, UPI, CRR, CRI, DRR, DRI [15]float64

	var AR = [15]float64{math.NaN(),
		1.00000000000000000e+00, 1.04166666666666667e-01,
		8.35503472222222222e-02, 1.28226574556327160e-01,
		2.91849026464140464e-01, 8.81627267443757652e-01,
		3.32140828186276754e+00, 1.49957629868625547e+01,
		7.89230130115865181e+01, 4.74451538868264323e+02,
		3.20749009089066193e+03, 2.40865496408740049e+04,
		1.98923119169509794e+05, 1.79190200777534383e+06}

	var BR = [15]float64{math.NaN(),
		1.00000000000000000e+00, -1.45833333333333333e-01,
		-9.87413194444444444e-02, -1.43312053915895062e-01,
		-3.17227202678413548e-01, -9.42429147957120249e-01,
		-3.51120304082635426e+00, -1.57272636203680451e+01,
		-8.22814390971859444e+01, -4.92355370523670524e+02,
		-3.31621856854797251e+03, -2.48276742452085896e+04,
		-2.04526587315129788e+05, -1.83844491706820990e+06}

	var C = [106]float64{math.NaN(),
		1.00000000000000000e+00, -2.08333333333333333e-01,
		1.25000000000000000e-01, 3.34201388888888889e-01,
		-4.01041666666666667e-01, 7.03125000000000000e-02,
		-1.02581259645061728e+00, 1.84646267361111111e+00,
		-8.91210937500000000e-01, 7.32421875000000000e-02,
		4.66958442342624743e+00, -1.12070026162229938e+01,
		8.78912353515625000e+00, -2.36408691406250000e+00,
		1.12152099609375000e-01, -2.82120725582002449e+01,
		8.46362176746007346e+01, -9.18182415432400174e+01,
		4.25349987453884549e+01, -7.36879435947963170e+00,
		2.27108001708984375e-01, 2.12570130039217123e+02,
		-7.65252468141181642e+02, 1.05999045252799988e+03,
		-6.99579627376132541e+02, 2.18190511744211590e+02,
		-2.64914304869515555e+01, 5.72501420974731445e-01,
		-1.91945766231840700e+03, 8.06172218173730938e+03,
		-1.35865500064341374e+04, 1.16553933368645332e+04,
		-5.30564697861340311e+03, 1.20090291321635246e+03,
		-1.08090919788394656e+02, 1.72772750258445740e+00,
		2.02042913309661486e+04, -9.69805983886375135e+04,
		1.92547001232531532e+05, -2.03400177280415534e+05,
		1.22200464983017460e+05, -4.11926549688975513e+04,
		7.10951430248936372e+03, -4.93915304773088012e+02,
		6.07404200127348304e+00, -2.42919187900551333e+05,
		1.31176361466297720e+06, -2.99801591853810675e+06,
		3.76327129765640400e+06, -2.81356322658653411e+06,
		1.26836527332162478e+06, -3.31645172484563578e+05,
		4.52187689813627263e+04, -2.49983048181120962e+03,
		2.43805296995560639e+01, 3.28446985307203782e+06,
		-1.97068191184322269e+07, 5.09526024926646422e+07,
		-7.41051482115326577e+07, 6.63445122747290267e+07,
		-3.75671766607633513e+07, 1.32887671664218183e+07,
		-2.78561812808645469e+06, 3.08186404612662398e+05,
		-1.38860897537170405e+04, 1.10017140269246738e+02,
		-4.93292536645099620e+07, 3.25573074185765749e+08,
		-9.39462359681578403e+08, 1.55359689957058006e+09,
		-1.62108055210833708e+09, 1.10684281682301447e+09,
		-4.95889784275030309e+08, 1.42062907797533095e+08,
		-2.44740627257387285e+07, 2.24376817792244943e+06,
		-8.40054336030240853e+04, 5.51335896122020586e+02,
		8.14789096118312115e+08, -5.86648149205184723e+09,
		1.86882075092958249e+10, -3.46320433881587779e+10,
		4.12801855797539740e+10, -3.30265997498007231e+10,
		1.79542137311556001e+10, -6.56329379261928433e+09,
		1.55927986487925751e+09, -2.25105661889415278e+08,
		1.73951075539781645e+07, -5.49842327572288687e+05,
		3.03809051092238427e+03, -1.46792612476956167e+10,
		1.14498237732025810e+11, -3.99096175224466498e+11,
		8.19218669548577329e+11, -1.09837515608122331e+12,
		1.00815810686538209e+12, -6.45364869245376503e+11,
		2.87900649906150589e+11, -8.78670721780232657e+10,
		1.76347306068349694e+10, -2.16716498322379509e+09,
		1.43157876718888981e+08, -3.87183344257261262e+06,
		1.82577554742931747e+04}

	var ALFA = [181]float64{math.NaN(),
		-4.44444444444444444e-03, -9.22077922077922078e-04,
		-8.84892884892884893e-05, 1.65927687832449737e-04,
		2.46691372741792910e-04, 2.65995589346254780e-04,
		2.61824297061500945e-04, 2.48730437344655609e-04,
		2.32721040083232098e-04, 2.16362485712365082e-04,
		2.00738858762752355e-04, 1.86267636637545172e-04,
		1.73060775917876493e-04, 1.61091705929015752e-04,
		1.50274774160908134e-04, 1.40503497391269794e-04,
		1.31668816545922806e-04, 1.23667445598253261e-04,
		1.16405271474737902e-04, 1.09798298372713369e-04,
		1.03772410422992823e-04, 9.82626078369363448e-05,
		9.32120517249503256e-05, 8.85710852478711718e-05,
		8.42963105715700223e-05, 8.03497548407791151e-05,
		7.66981345359207388e-05, 7.33122157481777809e-05,
		7.01662625163141333e-05, 6.72375633790160292e-05,
		6.93735541354588974e-04, 2.32241745182921654e-04,
		-1.41986273556691197e-05, -1.16444931672048640e-04,
		-1.50803558053048762e-04, -1.55121924918096223e-04,
		-1.46809756646465549e-04, -1.33815503867491367e-04,
		-1.19744975684254051e-04, -1.06184319207974020e-04,
		-9.37699549891194492e-05, -8.26923045588193274e-05,
		-7.29374348155221211e-05, -6.44042357721016283e-05,
		-5.69611566009369048e-05, -5.04731044303561628e-05,
		-4.48134868008882786e-05, -3.98688727717598864e-05,
		-3.55400532972042498e-05, -3.17414256609022480e-05,
		-2.83996793904174811e-05, -2.54522720634870566e-05,
		-2.28459297164724555e-05, -2.05352753106480604e-05,
		-1.84816217627666085e-05, -1.66519330021393806e-05,
		-1.50179412980119482e-05, -1.35554031379040526e-05,
		-1.22434746473858131e-05, -1.10641884811308169e-05,
		-3.54211971457743841e-04, -1.56161263945159416e-04,
		3.04465503594936410e-05, 1.30198655773242693e-04,
		1.67471106699712269e-04, 1.70222587683592569e-04,
		1.56501427608594704e-04, 1.36339170977445120e-04,
		1.14886692029825128e-04, 9.45869093034688111e-05,
		7.64498419250898258e-05, 6.07570334965197354e-05,
		4.74394299290508799e-05, 3.62757512005344297e-05,
		2.69939714979224901e-05, 1.93210938247939253e-05,
		1.30056674793963203e-05, 7.82620866744496661e-06,
		3.59257485819351583e-06, 1.44040049814251817e-07,
		-2.65396769697939116e-06, -4.91346867098485910e-06,
		-6.72739296091248287e-06, -8.17269379678657923e-06,
		-9.31304715093561232e-06, -1.02011418798016441e-05,
		-1.08805962510592880e-05, -1.13875481509603555e-05,
		-1.17519675674556414e-05, -1.19987364870944141e-05,
		3.78194199201772914e-04, 2.02471952761816167e-04,
		-6.37938506318862408e-05, -2.38598230603005903e-04,
		-3.10916256027361568e-04, -3.13680115247576316e-04,
		-2.78950273791323387e-04, -2.28564082619141374e-04,
		-1.75245280340846749e-04, -1.25544063060690348e-04,
		-8.22982872820208365e-05, -4.62860730588116458e-05,
		-1.72334302366962267e-05, 5.60690482304602267e-06,
		2.31395443148286800e-05, 3.62642745856793957e-05,
		4.58006124490188752e-05, 5.24595294959114050e-05,
		5.68396208545815266e-05, 5.94349820393104052e-05,
		6.06478527578421742e-05, 6.08023907788436497e-05,
		6.01577894539460388e-05, 5.89199657344698500e-05,
		5.72515823777593053e-05, 5.52804375585852577e-05,
		5.31063773802880170e-05, 5.08069302012325706e-05,
		4.84418647620094842e-05, 4.60568581607475370e-05,
		-6.91141397288294174e-04, -4.29976633058871912e-04,
		1.83067735980039018e-04, 6.60088147542014144e-04,
		8.75964969951185931e-04, 8.77335235958235514e-04,
		7.49369585378990637e-04, 5.63832329756980918e-04,
		3.68059319971443156e-04, 1.88464535514455599e-04,
		3.70663057664904149e-05, -8.28520220232137023e-05,
		-1.72751952869172998e-04, -2.36314873605872983e-04,
		-2.77966150694906658e-04, -3.02079514155456919e-04,
		-3.12594712643820127e-04, -3.12872558758067163e-04,
		-3.05678038466324377e-04, -2.93226470614557331e-04,
		-2.77255655582934777e-04, -2.59103928467031709e-04,
		-2.39784014396480342e-04, -2.20048260045422848e-04,
		-2.00443911094971498e-04, -1.81358692210970687e-04,
		-1.63057674478657464e-04, -1.45712672175205844e-04,
		-1.29425421983924587e-04, -1.14245691942445952e-04,
		1.92821964248775885e-03, 1.35592576302022234e-03,
		-7.17858090421302995e-04, -2.58084802575270346e-03,
		-3.49271130826168475e-03, -3.46986299340960628e-03,
		-2.82285233351310182e-03, -1.88103076404891354e-03,
		-8.89531718383947600e-04, 3.87912102631035228e-06,
		7.28688540119691412e-04, 1.26566373053457758e-03,
		1.62518158372674427e-03, 1.83203153216373172e-03,
		1.91588388990527909e-03, 1.90588846755546138e-03,
		1.82798982421825727e-03, 1.70389506421121530e-03,
		1.55097127171097686e-03, 1.38261421852276159e-03,
		1.20881424230064774e-03, 1.03676532638344962e-03,
		8.71437918068619115e-04, 7.16080155297701002e-04,
		5.72637002558129372e-04, 4.42089819465802277e-04,
		3.24724948503090564e-04, 2.20342042730246599e-04,
		1.28412898401353882e-04, 4.82005924552095464e-05}

	var BETA = [211]float64{math.NaN(),
		1.79988721413553309e-02, 5.59964911064388073e-03,
		2.88501402231132779e-03, 1.80096606761053941e-03,
		1.24753110589199202e-03, 9.22878876572938311e-04,
		7.14430421727287357e-04, 5.71787281789704872e-04,
		4.69431007606481533e-04, 3.93232835462916638e-04,
		3.34818889318297664e-04, 2.88952148495751517e-04,
		2.52211615549573284e-04, 2.22280580798883327e-04,
		1.97541838033062524e-04, 1.76836855019718004e-04,
		1.59316899661821081e-04, 1.44347930197333986e-04,
		1.31448068119965379e-04, 1.20245444949302884e-04,
		1.10449144504599392e-04, 1.01828770740567258e-04,
		9.41998224204237509e-05, 8.74130545753834437e-05,
		8.13466262162801467e-05, 7.59002269646219339e-05,
		7.09906300634153481e-05, 6.65482874842468183e-05,
		6.25146958969275078e-05, 5.88403394426251749e-05,
		-1.49282953213429172e-03, -8.78204709546389328e-04,
		-5.02916549572034614e-04, -2.94822138512746025e-04,
		-1.75463996970782828e-04, -1.04008550460816434e-04,
		-5.96141953046457895e-05, -3.12038929076098340e-05,
		-1.26089735980230047e-05, -2.42892608575730389e-07,
		8.05996165414273571e-06, 1.36507009262147391e-05,
		1.73964125472926261e-05, 1.98672978842133780e-05,
		2.14463263790822639e-05, 2.23954659232456514e-05,
		2.28967783814712629e-05, 2.30785389811177817e-05,
		2.30321976080909144e-05, 2.28236073720348722e-05,
		2.25005881105292418e-05, 2.20981015361991429e-05,
		2.16418427448103905e-05, 2.11507649256220843e-05,
		2.06388749782170737e-05, 2.01165241997081666e-05,
		1.95913450141179244e-05, 1.90689367910436740e-05,
		1.85533719641636667e-05, 1.80475722259674218e-05,
		5.52213076721292790e-04, 4.47932581552384646e-04,
		2.79520653992020589e-04, 1.52468156198446602e-04,
		6.93271105657043598e-05, 1.76258683069991397e-05,
		-1.35744996343269136e-05, -3.17972413350427135e-05,
		-4.18861861696693365e-05, -4.69004889379141029e-05,
		-4.87665447413787352e-05, -4.87010031186735069e-05,
		-4.74755620890086638e-05, -4.55813058138628452e-05,
		-4.33309644511266036e-05, -4.09230193157750364e-05,
		-3.84822638603221274e-05, -3.60857167535410501e-05,
		-3.37793306123367417e-05, -3.15888560772109621e-05,
		-2.95269561750807315e-05, -2.75978914828335759e-05,
		-2.58006174666883713e-05, -2.41308356761280200e-05,
		-2.25823509518346033e-05, -2.11479656768912971e-05,
		-1.98200638885294927e-05, -1.85909870801065077e-05,
		-1.74532699844210224e-05, -1.63997823854497997e-05,
		-4.74617796559959808e-04, -4.77864567147321487e-04,
		-3.20390228067037603e-04, -1.61105016119962282e-04,
		-4.25778101285435204e-05, 3.44571294294967503e-05,
		7.97092684075674924e-05, 1.03138236708272200e-04,
		1.12466775262204158e-04, 1.13103642108481389e-04,
		1.08651634848774268e-04, 1.01437951597661973e-04,
		9.29298396593363896e-05, 8.40293133016089978e-05,
		7.52727991349134062e-05, 6.69632521975730872e-05,
		5.92564547323194704e-05, 5.22169308826975567e-05,
		4.58539485165360646e-05, 4.01445513891486808e-05,
		3.50481730031328081e-05, 3.05157995034346659e-05,
		2.64956119950516039e-05, 2.29363633690998152e-05,
		1.97893056664021636e-05, 1.70091984636412623e-05,
		1.45547428261524004e-05, 1.23886640995878413e-05,
		1.04775876076583236e-05, 8.79179954978479373e-06,
		7.36465810572578444e-04, 8.72790805146193976e-04,
		6.22614862573135066e-04, 2.85998154194304147e-04,
		3.84737672879366102e-06, -1.87906003636971558e-04,
		-2.97603646594554535e-04, -3.45998126832656348e-04,
		-3.53382470916037712e-04, -3.35715635775048757e-04,
		-3.04321124789039809e-04, -2.66722723047612821e-04,
		-2.27654214122819527e-04, -1.89922611854562356e-04,
		-1.55058918599093870e-04, -1.23778240761873630e-04,
		-9.62926147717644187e-05, -7.25178327714425337e-05,
		-5.22070028895633801e-05, -3.50347750511900522e-05,
		-2.06489761035551757e-05, -8.70106096849767054e-06,
		1.13698686675100290e-06, 9.16426474122778849e-06,
		1.56477785428872620e-05, 2.08223629482466847e-05,
		2.48923381004595156e-05, 2.80340509574146325e-05,
		3.03987774629861915e-05, 3.21156731406700616e-05,
		-1.80182191963885708e-03, -2.43402962938042533e-03,
		-1.83422663549856802e-03, -7.62204596354009765e-04,
		2.39079475256927218e-04, 9.49266117176881141e-04,
		1.34467449701540359e-03, 1.48457495259449178e-03,
		1.44732339830617591e-03, 1.30268261285657186e-03,
		1.10351597375642682e-03, 8.86047440419791759e-04,
		6.73073208165665473e-04, 4.77603872856582378e-04,
		3.05991926358789362e-04, 1.60315694594721630e-04,
		4.00749555270613286e-05, -5.66607461635251611e-05,
		-1.32506186772982638e-04, -1.90296187989614057e-04,
		-2.32811450376937408e-04, -2.62628811464668841e-04,
		-2.82050469867598672e-04, -2.93081563192861167e-04,
		-2.97435962176316616e-04, -2.96557334239348078e-04,
		-2.91647363312090861e-04, -2.83696203837734166e-04,
		-2.73512317095673346e-04, -2.61750155806768580e-04,
		6.38585891212050914e-03, 9.62374215806377941e-03,
		7.61878061207001043e-03, 2.83219055545628054e-03,
		-2.09841352012720090e-03, -5.73826764216626498e-03,
		-7.70804244495414620e-03, -8.21011692264844401e-03,
		-7.65824520346905413e-03, -6.47209729391045177e-03,
		-4.99132412004966473e-03, -3.45612289713133280e-03,
		-2.01785580014170775e-03, -7.59430686781961401e-04,
		2.84173631523859138e-04, 1.10891667586337403e-03,
		1.72901493872728771e-03, 2.16812590802684701e-03,
		2.45357710494539735e-03, 2.61281821058334862e-03,
		2.67141039656276912e-03, 2.65203073395980430e-03,
		2.57411652877287315e-03, 2.45389126236094427e-03,
		2.30460058071795494e-03, 2.13684837686712662e-03,
		1.95896528478870911e-03, 1.77737008679454412e-03,
		1.59690280765839059e-03, 1.42111975664438546e-03}

	var GAMA = [31]float64{math.NaN(),
		6.29960524947436582e-01, 2.51984209978974633e-01,
		1.54790300415655846e-01, 1.10713062416159013e-01,
		8.57309395527394825e-02, 6.97161316958684292e-02,
		5.86085671893713576e-02, 5.04698873536310685e-02,
		4.42600580689154809e-02, 3.93720661543509966e-02,
		3.54283195924455368e-02, 3.21818857502098231e-02,
		2.94646240791157679e-02, 2.71581677112934479e-02,
		2.51768272973861779e-02, 2.34570755306078891e-02,
		2.19508390134907203e-02, 2.06210828235646240e-02,
		1.94388240897880846e-02, 1.83810633800683158e-02,
		1.74293213231963172e-02, 1.65685837786612353e-02,
		1.57865285987918445e-02, 1.50729501494095594e-02,
		1.44193250839954639e-02, 1.38184805735341786e-02,
		1.32643378994276568e-02, 1.27517121970498651e-02,
		1.22761545318762767e-02, 1.18338262398482403e-02}

	RFNU = 1.0e0 / FNU

	// OVERFLOW TEST (Z/FNU TOO SMALL)
	TEST = D1MACH[1] * 1.0e+3
	AC = FNU * TEST
	if math.Abs(ZR) > AC || math.Abs(ZI) > AC {
		goto L15
	}

	ZETA1R = 2.0e0*math.Abs(math.Log(TEST)) + FNU
	ZETA1I = 0.0e0
	ZETA2R = FNU
	ZETA2I = 0.0e0
	PHIR = 1.0e0
	PHII = 0.0e0
	ARGR = 1.0e0
	ARGI = 0.0e0
	return ZR, ZI, FNU, IPMTR, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R, ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI
L15:
	ZBR = ZR * RFNU
	ZBI = ZI * RFNU
	RFNU2 = RFNU * RFNU

	// COMPUTE IN THE FOURTH QUADRANT
	FN13 = math.Pow(FNU, EX1)
	FN23 = FN13 * FN13
	RFN13 = 1.0e0 / FN13
	W2R = CONER - ZBR*ZBR + ZBI*ZBI
	W2I = CONEI - ZBR*ZBI - ZBR*ZBI
	AW2 = ZABS(W2R, W2I)
	if AW2 > 0.25e0 {
		goto L130
	}

	// POWER SERIES FOR CABS(W2) <= 0.25e0
	K = 1
	PR[1] = CONER
	PI[1] = CONEI
	SUMAR = GAMA[1]
	SUMAI = ZEROI
	AP[1] = 1.0e0
	if AW2 < TOL {
		goto L20
	}
	for K = 2; K <= 30; K++ {
		PR[K] = PR[K-1]*W2R - PI[K-1]*W2I
		PI[K] = PR[K-1]*W2I + PI[K-1]*W2R
		SUMAR = SUMAR + PR[K]*GAMA[K]
		SUMAI = SUMAI + PI[K]*GAMA[K]
		AP[K] = AP[K-1] * AW2
		if AP[K] < TOL {
			goto L20
		}
	}
	K = 30
L20:
	KMAX = K
	ZETAR = W2R*SUMAR - W2I*SUMAI
	ZETAI = W2R*SUMAI + W2I*SUMAR
	ARGR = ZETAR * FN23
	ARGI = ZETAI * FN23
	ZAR, ZAI = ZSQRT(SUMAR, SUMAI)
	STR, STI = ZSQRT(W2R, W2I)
	ZETA2R = STR * FNU
	ZETA2I = STI * FNU
	STR = CONER + EX2*(ZETAR*ZAR-ZETAI*ZAI)
	STI = CONEI + EX2*(ZETAR*ZAI+ZETAI*ZAR)
	ZETA1R = STR*ZETA2R - STI*ZETA2I
	ZETA1I = STR*ZETA2I + STI*ZETA2R
	ZAR = ZAR + ZAR
	ZAI = ZAI + ZAI
	STR, STI = ZSQRT(ZAR, ZAI)
	PHIR = STR * RFN13
	PHII = STI * RFN13
	if IPMTR == 1 {
		goto L120
	}

	// SUM SERIES FOR ASUM AND BSUM
	SUMBR = ZEROR
	SUMBI = ZEROI
	for K = 1; K <= KMAX; K++ {
		SUMBR = SUMBR + PR[K]*BETA[K]
		SUMBI = SUMBI + PI[K]*BETA[K]
	}
	ASUMR = ZEROR
	ASUMI = ZEROI
	BSUMR = SUMBR
	BSUMI = SUMBI
	L1 = 0
	L2 = 30
	BTOL = TOL * (math.Abs(BSUMR) + math.Abs(BSUMI))
	ATOL = TOL
	PP = 1.0e0
	IAS = 0
	IBS = 0
	if RFNU2 < TOL {
		goto L110
	}
	for IS = 2; IS <= 7; IS++ {
		ATOL = ATOL / RFNU2
		PP = PP * RFNU2
		if IAS == 1 {
			goto L60
		}
		SUMAR = ZEROR
		SUMAI = ZEROI
		for K = 1; K <= KMAX; K++ {
			M = L1 + K
			SUMAR = SUMAR + PR[K]*ALFA[M]
			SUMAI = SUMAI + PI[K]*ALFA[M]
			if AP[K] < ATOL {
				goto L50
			}
		}
	L50:
		ASUMR = ASUMR + SUMAR*PP
		ASUMI = ASUMI + SUMAI*PP
		if PP < TOL {
			IAS = 1
		}
	L60:
		if IBS == 1 {
			goto L90
		}
		SUMBR = ZEROR
		SUMBI = ZEROI
		for K = 1; K <= KMAX; K++ {
			M = L2 + K
			SUMBR = SUMBR + PR[K]*BETA[M]
			SUMBI = SUMBI + PI[K]*BETA[M]
			if AP[K] < ATOL {
				goto L80
			}
		}
	L80:
		BSUMR = BSUMR + SUMBR*PP
		BSUMI = BSUMI + SUMBI*PP
		if PP < BTOL {
			IBS = 1
		}
	L90:
		if IAS == 1 && IBS == 1 {
			goto L110
		}
		L1 = L1 + 30
		L2 = L2 + 30
	}
L110:
	ASUMR = ASUMR + CONER
	PP = RFNU * RFN13
	BSUMR = BSUMR * PP
	BSUMI = BSUMI * PP
L120:
	return ZR, ZI, FNU, IPMTR, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R, ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI

	// CABS(W2) > 0.25e0
L130:
	WR, WI = ZSQRT(W2R, W2I)
	if WR < 0.0e0 {
		WR = 0.0e0
	}
	if WI < 0.0e0 {
		WI = 0.0e0
	}
	STR = CONER + WR
	STI = WI
	ZAR, ZAI = ZDIV(STR, STI, ZBR, ZBI)
	ZCR, ZCI = ZLOG(ZAR, ZAI)
	if ZCI < 0.0e0 {
		ZCI = 0.0e0
	}
	if ZCI > HPI {
		ZCI = HPI
	}
	if ZCR < 0.0e0 {
		ZCR = 0.0e0
	}
	ZTHR = (ZCR - WR) * 1.5e0
	ZTHI = (ZCI - WI) * 1.5e0
	ZETA1R = ZCR * FNU
	ZETA1I = ZCI * FNU
	ZETA2R = WR * FNU
	ZETA2I = WI * FNU
	AZTH = ZABS(ZTHR, ZTHI)
	ANG = THPI
	if ZTHR >= 0.0e0 && ZTHI < 0.0e0 {
		goto L140
	}
	ANG = HPI
	if ZTHR == 0.0e0 {
		goto L140
	}
	ANG = math.Atan(ZTHI / ZTHR)
	if ZTHR < 0.0e0 {
		ANG = ANG + GPI
	}
L140:
	PP = math.Pow(AZTH, EX2)
	ANG = ANG * EX2
	ZETAR = PP * math.Cos(ANG)
	ZETAI = PP * math.Sin(ANG)
	if ZETAI < 0.0e0 {
		ZETAI = 0.0e0
	}
	ARGR = ZETAR * FN23
	ARGI = ZETAI * FN23
	RTZTR, RTZTI = ZDIV(ZTHR, ZTHI, ZETAR, ZETAI)
	ZAR, ZAI = ZDIV(RTZTR, RTZTI, WR, WI)
	TZAR = ZAR + ZAR
	TZAI = ZAI + ZAI
	STR, STI = ZSQRT(TZAR, TZAI)
	PHIR = STR * RFN13
	PHII = STI * RFN13
	if IPMTR == 1 {
		goto L120
	}
	RAW = 1.0e0 / math.Sqrt(AW2)
	STR = WR * RAW
	STI = -WI * RAW
	TFNR = STR * RFNU * RAW
	TFNI = STI * RFNU * RAW
	RAZTH = 1.0e0 / AZTH
	STR = ZTHR * RAZTH
	STI = -ZTHI * RAZTH
	RZTHR = STR * RAZTH * RFNU
	RZTHI = STI * RAZTH * RFNU
	ZCR = RZTHR * AR[2]
	ZCI = RZTHI * AR[2]
	RAW2 = 1.0e0 / AW2
	STR = W2R * RAW2
	STI = -W2I * RAW2
	T2R = STR * RAW2
	T2I = STI * RAW2
	STR = T2R*C[2] + C[3]
	STI = T2I * C[2]
	UPR[2] = STR*TFNR - STI*TFNI
	UPI[2] = STR*TFNI + STI*TFNR
	BSUMR = UPR[2] + ZCR
	BSUMI = UPI[2] + ZCI
	ASUMR = ZEROR
	ASUMI = ZEROI
	if RFNU < TOL {
		goto L220
	}
	PRZTHR = RZTHR
	PRZTHI = RZTHI
	PTFNR = TFNR
	PTFNI = TFNI
	UPR[1] = CONER
	UPI[1] = CONEI
	PP = 1.0e0
	BTOL = TOL * (math.Abs(BSUMR) + math.Abs(BSUMI))
	KS = 0
	KP1 = 2
	L = 3
	IAS = 0
	IBS = 0
	for LR = 2; LR < 12; LR = LR + 2 {
		LRP1 = LR + 1

		// COMPUTE TWO ADDITIONAL CR, DR, AND UP FOR TWO MORE TERMS IN
		// NEXT SUMA AND SUMB
		for K = LR; K <= LRP1; K++ {
			KS = KS + 1
			KP1 = KP1 + 1
			L = L + 1
			ZAR = C[L]
			ZAI = ZEROI
			for J = 2; J <= KP1; J++ {
				L = L + 1
				STR = ZAR*T2R - T2I*ZAI + C[L]
				ZAI = ZAR*T2I + ZAI*T2R
				ZAR = STR
			}
			STR = PTFNR*TFNR - PTFNI*TFNI
			PTFNI = PTFNR*TFNI + PTFNI*TFNR
			PTFNR = STR
			UPR[KP1] = PTFNR*ZAR - PTFNI*ZAI
			UPI[KP1] = PTFNI*ZAR + PTFNR*ZAI
			CRR[KS] = PRZTHR * BR[KS+1]
			CRI[KS] = PRZTHI * BR[KS+1]
			STR = PRZTHR*RZTHR - PRZTHI*RZTHI
			PRZTHI = PRZTHR*RZTHI + PRZTHI*RZTHR
			PRZTHR = STR
			DRR[KS] = PRZTHR * AR[KS+2]
			DRI[KS] = PRZTHI * AR[KS+2]
		}
		PP = PP * RFNU2
		if IAS == 1 {
			goto L180
		}
		SUMAR = UPR[LRP1]
		SUMAI = UPI[LRP1]
		JU = LRP1
		for JR = 1; JR <= LR; JR++ {
			JU = JU - 1
			SUMAR = SUMAR + CRR[JR]*UPR[JU] - CRI[JR]*UPI[JU]
			SUMAI = SUMAI + CRR[JR]*UPI[JU] + CRI[JR]*UPR[JU]
		}
		ASUMR = ASUMR + SUMAR
		ASUMI = ASUMI + SUMAI
		TEST = math.Abs(SUMAR) + math.Abs(SUMAI)
		if PP < TOL && TEST < TOL {
			IAS = 1
		}
	L180:
		if IBS == 1 {
			goto L200
		}
		SUMBR = UPR[LR+2] + UPR[LRP1]*ZCR - UPI[LRP1]*ZCI
		SUMBI = UPI[LR+2] + UPR[LRP1]*ZCI + UPI[LRP1]*ZCR
		JU = LRP1
		for JR = 1; JR <= LR; JR++ {
			JU = JU - 1
			SUMBR = SUMBR + DRR[JR]*UPR[JU] - DRI[JR]*UPI[JU]
			SUMBI = SUMBI + DRR[JR]*UPI[JU] + DRI[JR]*UPR[JU]
		}
		BSUMR = BSUMR + SUMBR
		BSUMI = BSUMI + SUMBI
		TEST = math.Abs(SUMBR) + math.Abs(SUMBI)
		if PP < BTOL && TEST < BTOL {
			IBS = 1
		}
	L200:
		if IAS == 1 && IBS == 1 {
			goto L220
		}
	}
L220:
	ASUMR = ASUMR + CONER
	STR = -BSUMR * RFN13
	STI = -BSUMI * RFN13
	BSUMR, BSUMI = ZDIV(STR, STI, RTZTR, RTZTI)
	goto L120
}

// ZUOIK COMPUTES THE LEADING TERMS OF THE UNIFORM ASYMPTOTIC
// EXPANSIONS FOR THE I AND K FUNCTIONS AND COMPARES THEM
// (IN LOGARITHMIC FORM) TO ALIM AND ELIM FOR OVER AND UNDERFLOW
// WHERE ALIM < ELIM. if THE MAGNITUDE, BASED ON THE LEADING
// EXPONENTIAL, IS LESS THAN ALIM OR GREATER THAN -ALIM, THEN
// THE RESULT IS ON SCALE. if NOT, THEN A REFINED TEST USING OTHER
// MULTIPLIERS (IN LOGARITHMIC FORM) IS MADE BASED ON ELIM. HERE
// EXP(-ELIM)=SMALLEST MACHINE NUMBER*1.0E+3 AND EXP(-ALIM)=
// EXP(-ELIM)/TOL
//
// IKFLG=1 MEANS THE I SEQUENCE IS TESTED
//   =2 MEANS THE K SEQUENCE IS TESTED
// NUF = 0 MEANS THE LAST MEMBER OF THE SEQUENCE IS ON SCALE
//  =-1 MEANS AN OVERFLOW WOULD OCCUR
// IKFLG=1 AND NUF > 0 MEANS THE LAST NUF Y VALUES WERE SET TO ZERO
//      THE FIRST N-NUF VALUES MUST BE SET BY ANOTHER ROUTINE
// IKFLG=2 AND NUF == N MEANS ALL Y VALUES WERE SET TO ZERO
// IKFLG=2 AND 0 < NUF < N NOT CONSIDERED. Y MUST BE SET BY
//      ANOTHER ROUTINE
func ZUOIK(ZR float64, ZI float64, FNU float64, KODE int, IKFLG int, N int, YR []float64, YI []float64, NUF int, TOL float64, ELIM float64, ALIM float64) (
	float64, float64, float64, int, int, int, []float64, []float64, int, float64, float64, float64) {

	const (
		ZEROR = 0.0e0
		ZEROI = 0.0e0
		AIC   = 1.265512123484645396e+00
	)

	var AARG, APHI, ARGI, ARGR, ASUMI, ASUMR, ASCLE, AX, AY, BSUMI, BSUMR, CZI, CZR, FNN,
		GNN, GNU, PHII, PHIR, RCZ, STR, STI, SUMI, SUMR, ZBI, ZBR, ZETA1I, ZETA1R, ZETA2I, ZETA2R,
		ZNI, ZNR, ZRI, ZRR float64
	var I, IFORM, INIT, NN, NW int

	var CWRKR = []float64{math.NaN(), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
	var CWRKI = []float64{math.NaN(), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}

	NUF = 0
	NN = N
	ZRR = ZR
	ZRI = ZI
	if ZR >= 0.0e0 {
		goto L10
	}
	ZRR = -ZR
	ZRI = -ZI
L10:
	ZBR = ZRR
	ZBI = ZRI
	AX = math.Abs(ZR) * 1.7321e0
	AY = math.Abs(ZI)
	IFORM = 1
	if AY > AX {
		IFORM = 2
	}
	GNU = math.Max(FNU, 1.0e0)
	if IKFLG == 1 {
		goto L20
	}
	FNN = float64(float32(NN))
	GNN = FNU + FNN - 1.0e0
	GNU = math.Max(GNN, FNN)
L20:
	// ONLY THE MAGNITUDE OF ARG AND PHI ARE NEEDED ALONG WITH THE
	// REAL PARTS OF ZETA1, ZETA2 AND ZB. NO ATTEMPT IS MADE TO GET
	// THE SIGN OF THE IMAGINARY PART CORRECT.
	if IFORM == 2 {
		goto L30
	}
	INIT = 0
	ZRR, ZRI, GNU, IKFLG, _, TOL, INIT, PHIR, PHII, ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI = ZUNIK(ZRR, ZRI, GNU, IKFLG, 1, TOL, INIT, PHIR, PHII, ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
	CZR = -ZETA1R + ZETA2R
	CZI = -ZETA1I + ZETA2I
	goto L50
L30:
	ZNR = ZRI
	ZNI = -ZRR
	if ZI > 0.0e0 {
		goto L40
	}
	ZNR = -ZNR
L40:
	ZNR, ZNI, GNU, _, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R, ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI = ZUNHJ(ZNR, ZNI, GNU, 1, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R, ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
	CZR = -ZETA1R + ZETA2R
	CZI = -ZETA1I + ZETA2I
	AARG = ZABS(ARGR, ARGI)
L50:
	if KODE == 1 {
		goto L60
	}
	CZR = CZR - ZBR
	CZI = CZI - ZBI
L60:
	if IKFLG == 1 {
		goto L70
	}
	CZR = -CZR
	CZI = -CZI
L70:
	APHI = ZABS(PHIR, PHII)
	RCZ = CZR

	// OVERFLOW TEST
	if RCZ > ELIM {
		goto L210
	}
	if RCZ < ALIM {
		goto L80
	}
	RCZ = RCZ + math.Log(APHI)
	if IFORM == 2 {
		RCZ = RCZ - 0.25e0*math.Log(AARG) - AIC
	}
	if RCZ > ELIM {
		goto L210
	}
	goto L130
L80:

	// UNDERFLOW TEST
	if RCZ < -ELIM {
		goto L90
	}
	if RCZ > -ALIM {
		goto L130
	}
	RCZ = RCZ + math.Log(APHI)
	if IFORM == 2 {
		RCZ = RCZ - 0.25e0*math.Log(AARG) - AIC
	}
	if RCZ > -ELIM {
		goto L110
	}
L90:
	for I = 1; I <= NN; I++ {
		YR[I] = ZEROR
		YI[I] = ZEROI
	}
	NUF = NN
	return ZR, ZI, FNU, KODE, IKFLG, N, YR, YI, NUF, TOL, ELIM, ALIM
L110:
	ASCLE = 1000.0e0 * D1MACH[1] / TOL
	STR, STI = ZLOG(PHIR, PHII)
	CZR = CZR + STR
	CZI = CZI + STI
	if IFORM == 1 {
		goto L120
	}
	STR, STI = ZLOG(ARGR, ARGI)
	CZR = CZR - 0.25e0*STR - AIC
	CZI = CZI - 0.25e0*STI
L120:
	AX = math.Exp(RCZ) / TOL
	AY = CZI
	CZR = AX * math.Cos(AY)
	CZI = AX * math.Sin(AY)
	CZR, CZI, NW, ASCLE, TOL = ZUCHK(CZR, CZI, NW, ASCLE, TOL)
	if NW != 0 {
		goto L90
	}
L130:
	if IKFLG == 2 {
		return ZR, ZI, FNU, KODE, IKFLG, N, YR, YI, NUF, TOL, ELIM, ALIM
	}
	if N == 1 {
		return ZR, ZI, FNU, KODE, IKFLG, N, YR, YI, NUF, TOL, ELIM, ALIM
	}

	// SET UNDERFLOWS ON I SEQUENCE
L140:
	GNU = FNU + float64(float32(NN-1))
	if IFORM == 2 {
		goto L150
	}
	INIT = 0
	ZRR, ZRI, GNU, IKFLG, _, TOL, INIT, PHIR, PHII, ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI = ZUNIK(ZRR, ZRI, GNU, IKFLG, 1, TOL, INIT, PHIR, PHII, ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
	CZR = -ZETA1R + ZETA2R
	CZI = -ZETA1I + ZETA2I
	goto L160
L150:
	ZNR, ZNI, GNU, _, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R, ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI = ZUNHJ(ZNR, ZNI, GNU, 1, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R, ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
	CZR = -ZETA1R + ZETA2R
	CZI = -ZETA1I + ZETA2I
	AARG = ZABS(ARGR, ARGI)
L160:
	if KODE == 1 {
		goto L170
	}
	CZR = CZR - ZBR
	CZI = CZI - ZBI
L170:
	APHI = ZABS(PHIR, PHII)
	RCZ = CZR
	if RCZ < -ELIM {
		goto L180
	}
	if RCZ > -ALIM {
		return ZR, ZI, FNU, KODE, IKFLG, N, YR, YI, NUF, TOL, ELIM, ALIM
	}
	RCZ = RCZ + math.Log(APHI)
	if IFORM == 2 {
		RCZ = RCZ - 0.25e0*math.Log(AARG) - AIC
	}
	if RCZ > -ELIM {
		goto L190
	}
L180:
	YR[NN] = ZEROR
	YI[NN] = ZEROI
	NN = NN - 1
	NUF = NUF + 1
	if NN == 0 {
		return ZR, ZI, FNU, KODE, IKFLG, N, YR, YI, NUF, TOL, ELIM, ALIM
	}
	goto L140
L190:
	ASCLE = 1000.0e0 * D1MACH[1] / TOL
	STR, STI = ZLOG(PHIR, PHII)
	CZR = CZR + STR
	CZI = CZI + STI
	if IFORM == 1 {
		goto L200
	}
	STR, STI = ZLOG(ARGR, ARGI)
	CZR = CZR - 0.25e0*STR - AIC
	CZI = CZI - 0.25e0*STI
L200:
	AX = math.Exp(RCZ) / TOL
	AY = CZI
	CZR = AX * math.Cos(AY)
	CZI = AX * math.Sin(AY)
	CZR, CZI, NW, ASCLE, TOL = ZUCHK(CZR, CZI, NW, ASCLE, TOL)
	if NW != 0 {
		goto L180
	}
	return ZR, ZI, FNU, KODE, IKFLG, N, YR, YI, NUF, TOL, ELIM, ALIM
L210:
	NUF = -1
	return ZR, ZI, FNU, KODE, IKFLG, N, YR, YI, NUF, TOL, ELIM, ALIM
}

//  ZBIRY COMPUTES THE COMPLEX AIRY FUNCTION BI(Z) OR
//  ITS DERIVATIVE DBI(Z)/DZ ON ID=0 OR ID=1 RESPECTIVELY FOR KODE=1
//
//  FOR KODE=2, A SCALING OPTION CEXP(-AXZTA)*BI(Z) OR CEXP(-AXZTA)*
//  DBI(Z)/DZ IS PROVIDED TO REMOVE THE EXPONENTIAL BEHAVIOR IN
//  BOTH THE LEFT AND RIGHT HALF PLANES WHERE
//  ZTA=(2/3)*Z*CSQRT(Z)=CMPLX(XZTA,YZTA) AND AXZTA=ABS(XZTA).
//  DEFINTIONS AND NOTATION ARE FOUND IN THE NBS HANDBOOK OF
//  MATHEMATICAL FUNCTIONS (REF. 1).
//
//  INPUT      ZR,ZI ARE DOUBLE PRECISION
//    ZR,ZI  - Z=CMPLX(ZR,ZI)
//    ID     - ORDER OF DERIVATIVE, ID=0 OR ID=1
//    KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
//             KODE= 1  return S
//                      BI=BI(Z)                 ON ID=0 OR
//                      BI=DBI(Z)/DZ             ON ID=1
//                 = 2  return S
//                      BI=CEXP(-AXZTA)*BI(Z)     ON ID=0 OR
//                      BI=CEXP(-AXZTA)*DBI(Z)/DZ ON ID=1 WHERE
//                      ZTA=(2/3)*Z*CSQRT(Z)=CMPLX(XZTA,YZTA)
//                      AND AXZTA=ABS(XZTA)
//
//  OUTPUT     BIR,BII ARE DOUBLE PRECISION
//    BIR,BII- COMPLEX ANSWER DEPENDING ON THE CHOICES FOR ID AND
//             KODE
//    IERR   - ERROR FLAG
//             IERR=0, NORMAL return  - COMPUTATION COMPLETED
//             IERR=1, INPUT ERROR   - NO COMPUTATION
//             IERR=2, OVERFLOW      - NO COMPUTATION, REAL(Z)
//                     TOO LARGE ON KODE=1
//             IERR=3, CABS(Z) LARGE      - COMPUTATION COMPLETED
//                     LOSSES OF SIGNIFCANCE BY ARGUMENT REDUCTION
//                     PRODUCE LESS THAN HALF OF MACHINE ACCURACY
//             IERR=4, CABS(Z) TOO LARGE  - NO COMPUTATION
//                     COMPLETE LOSS OF ACCURACY BY ARGUMENT
//                     REDUCTION
//             IERR=5, ERROR              - NO COMPUTATION,
//                     ALGORITHM TERMINATION CONDITION NOT MET
func ZBIRY(ZR float64, ZI float64, ID int, KODE int) (float64, float64, int) {

	//  BI AND DBI ARE COMPUTED FOR CABS(Z) > 1.0 FROM THE I BESSEL
	//  FUNCTIONS BY
	//
	//         BI(Z)=C*SQRT(Z)*( I(-1/3,ZTA) + I(1/3,ZTA) )
	//        DBI(Z)=C *  Z  * ( I(-2/3,ZTA) + I(2/3,ZTA) )
	//                        C=1.0/SQRT(3.0)
	//                      ZTA=(2/3)*Z**(3/2)
	//
	//  WITH THE POWER SERIES FOR CABS(Z) <= 1.0.
	//
	//  IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
	//  MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z IS LARGE, LOSSES
	//  OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR. CONSEQUENTLY, IF
	//  THE MAGNITUDE OF ZETA=(2/3)*Z**1.5 EXCEEDS U1=SQRT(0.5/UR),
	//  THEN LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR
	//  FLAG IERR=3 IS TRIGGERED WHERE UR=math.Max(D1MACH(4),1.0e-18) IS
	//  DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
	//  ALSO, if THE MAGNITUDE OF ZETA IS LARGER THAN U2=0.5/UR, THEN
	//  ALL SIGNIFICANCE IS LOST AND IERR=4. IN ORDER TO USE THE INT
	//  FUNCTION, ZETA MUST BE FURTHER RESTRICTED NOT TO EXCEED THE
	//  LARGEST INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF ZETA
	//  MUST BE RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2,
	//  AND U3 ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE
	//  PRECISION ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE
	//  PRECISION ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMIT-
	//  ING IN THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT THE MAG-
	//  NITUDE OF Z CANNOT EXCEED 3.1E+4 IN SINGLE AND 2.1E+6 IN
	//  DOUBLE PRECISION ARITHMETIC. THIS ALSO MEANS THAT ONE CAN
	//  EXPECT TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES,
	//  NO DIGITS IN SINGLE PRECISION AND ONLY 7 DIGITS IN DOUBLE
	//  PRECISION ARITHMETIC. SIMILAR CONSIDERATIONS HOLD FOR OTHER
	//  MACHINES.
	//
	//  THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
	//  BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
	//  ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
	//  SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
	//  ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
	//  ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
	//  CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
	//  HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
	//  ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
	//  SEVERAL ORDERS OF MAGNITUDE. if ONE COMPONENT IS 10**K LARGER
	//  THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
	//  0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
	//  THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
	//  COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
	//  BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
	//  COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
	//  MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
	//  THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
	//  OR -PI/2+P.
	//
	//        HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
	//          AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
	//          COMMERCE, 1955.
	//
	//        COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
	//          AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
	//
	//        A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
	//          ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
	//          1018, MAY, 1985
	//
	//        A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
	//          ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, TRANS.
	//          MATH. SOFTWARE, 1986

	const (
		CONER = 1.0e0
		CONEI = 0.0e0
		TTH   = 6.66666666666666667e-01
		C1    = 6.14926627446000736e-01
		C2    = 4.48288357353826359e-01
		COEF  = 5.77350269189625765e-01
	)
	var AA, AD, AK, ALIM, ATRM, AZ, AZ3, BB, BII, BIR, BK, CC, CK, CSQI, CSQR, DIG, DK,
		D1, D2, EAA, ELIM, FID, FMR, FNU, FNUL, RL, R1M5, SFAC, STI, STR, S1I, S1R, S2I, S2R,
		TOL, TRM1I, TRM1R, TRM2I, TRM2R, ZTAI, ZTAR, Z3I, Z3R float64
	var IERR, K, K1, K2, NZ int

	CYR := []float64{math.NaN(), 0, 0}
	CYI := []float64{math.NaN(), 0, 0}

	IERR = 0
	NZ = 0
	if ID < 0 || ID > 1 {
		IERR = 1
	}
	if KODE < 1 || KODE > 2 {
		IERR = 1
	}
	if IERR != 0 {
		return BIR, BII, IERR
	}
	AZ = ZABS(ZR, ZI)
	TOL = math.Max(D1MACH[4], 1.0e-18)
	FID = float64(float32(ID))
	if AZ > 1.0e0 {
		goto L70
	}

	// POWER SERIES FOR CABS(Z) <= 1.
	S1R = CONER
	S1I = CONEI
	S2R = CONER
	S2I = CONEI
	if AZ < TOL {
		goto L130
	}
	AA = AZ * AZ
	if AA < TOL/AZ {
		goto L40
	}
	TRM1R = CONER
	TRM1I = CONEI
	TRM2R = CONER
	TRM2I = CONEI
	ATRM = 1.0e0
	STR = ZR*ZR - ZI*ZI
	STI = ZR*ZI + ZI*ZR
	Z3R = STR*ZR - STI*ZI
	Z3I = STR*ZI + STI*ZR
	AZ3 = AZ * AA
	AK = 2.0e0 + FID
	BK = 3.0e0 - FID - FID
	CK = 4.0e0 - FID
	DK = 3.0e0 + FID + FID
	D1 = AK * DK
	D2 = BK * CK
	AD = math.Min(D1, D2)
	AK = 24.0e0 + 9.0e0*FID
	BK = 30.0e0 - 9.0e0*FID
	for K = 1; K <= 25; K++ {
		STR = (TRM1R*Z3R - TRM1I*Z3I) / D1
		TRM1I = (TRM1R*Z3I + TRM1I*Z3R) / D1
		TRM1R = STR
		S1R = S1R + TRM1R
		S1I = S1I + TRM1I
		STR = (TRM2R*Z3R - TRM2I*Z3I) / D2
		TRM2I = (TRM2R*Z3I + TRM2I*Z3R) / D2
		TRM2R = STR
		S2R = S2R + TRM2R
		S2I = S2I + TRM2I
		ATRM = ATRM * AZ3 / AD
		D1 = D1 + AK
		D2 = D2 + BK
		AD = math.Min(D1, D2)
		if ATRM < TOL*AD {
			goto L40
		}
		AK = AK + 18.0e0
		BK = BK + 18.0e0
	}
L40:
	if ID == 1 {
		goto L50
	}
	BIR = C1*S1R + C2*(ZR*S2R-ZI*S2I)
	BII = C1*S1I + C2*(ZR*S2I+ZI*S2R)
	if KODE == 1 {
		return BIR, BII, IERR
	}
	STR, STI = ZSQRT(ZR, ZI)
	ZTAR = TTH * (ZR*STR - ZI*STI)
	ZTAI = TTH * (ZR*STI + ZI*STR)
	AA = ZTAR
	AA = -math.Abs(AA)
	EAA = math.Exp(AA)
	BIR = BIR * EAA
	BII = BII * EAA
	return BIR, BII, IERR
L50:
	BIR = S2R * C2
	BII = S2I * C2
	if AZ <= TOL {
		goto L60
	}
	CC = C1 / (1.0e0 + FID)
	STR = S1R*ZR - S1I*ZI
	STI = S1R*ZI + S1I*ZR
	BIR = BIR + CC*(STR*ZR-STI*ZI)
	BII = BII + CC*(STR*ZI+STI*ZR)
L60:
	if KODE == 1 {
		return BIR, BII, IERR
	}
	STR, STI = ZSQRT(ZR, ZI)
	ZTAR = TTH * (ZR*STR - ZI*STI)
	ZTAI = TTH * (ZR*STI + ZI*STR)
	AA = ZTAR
	AA = -math.Abs(AA)
	EAA = math.Exp(AA)
	BIR = BIR * EAA
	BII = BII * EAA
	return BIR, BII, IERR

	// CASE FOR CABS(Z) > 1.0
L70:
	FNU = (1.0e0 + FID) / 3.0e0

	// SET PARAMETERS RELATED TO MACHINE CONSTANTS.
	// TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
	// ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
	// EXP(-ELIM) < EXP(-ALIM)=EXP(-ELIM)/TOL    AND
	// EXP(ELIM) > EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
	// UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
	// RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
	// DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
	// FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU.
	K1 = I1MACH[15]
	K2 = I1MACH[16]
	R1M5 = D1MACH[5]
	K = MIN(ABS(K1), ABS(K2))
	ELIM = 2.303e0 * (float64(float32(K))*R1M5 - 3.0e0)
	K1 = I1MACH[14] - 1
	AA = R1M5 * float64(float32(K1))
	DIG = math.Min(AA, 18.0e0)
	AA = AA * 2.303e0
	ALIM = ELIM + math.Max(-AA, -41.45e0)
	RL = 1.2e0*DIG + 3.0e0
	FNUL = 10.0e0 + 6.0e0*(DIG-3.0e0)

	// TEST FOR RANGE
	AA = 0.5e0 / TOL
	BB = float64(float32(I1MACH[9])) * 0.5e0
	AA = math.Min(AA, BB)
	AA = math.Pow(AA, TTH)
	if AZ > AA {
		goto L260
	}
	AA = math.Sqrt(AA)
	if AZ > AA {
		IERR = 3
	}
	CSQR, CSQI = ZSQRT(ZR, ZI)
	ZTAR = TTH * (ZR*CSQR - ZI*CSQI)
	ZTAI = TTH * (ZR*CSQI + ZI*CSQR)

	// RE(ZTA) <= 0 WHEN RE(Z) < 0, ESPECIALLY WHEN IM(Z) IS SMALL
	SFAC = 1.0e0
	AK = ZTAI
	if ZR >= 0.0e0 {
		goto L80
	}
	BK = ZTAR
	CK = -math.Abs(BK)
	ZTAR = CK
	ZTAI = AK
L80:
	if ZI != 0.0e0 || ZR > 0.0e0 {
		goto L90
	}
	ZTAR = 0.0e0
	ZTAI = AK
L90:
	AA = ZTAR
	if KODE == 2 {
		goto L100
	}

	// OVERFLOW TEST
	BB = math.Abs(AA)
	if BB < ALIM {
		goto L100
	}
	BB = BB + 0.25e0*math.Log(AZ)
	SFAC = TOL
	if BB > ELIM {
		goto L190
	}
L100:
	FMR = 0.0e0
	if AA >= 0.0e0 && ZR > 0.0e0 {
		goto L110
	}
	FMR = math.Pi
	if ZI < 0.0e0 {
		FMR = -math.Pi
	}
	ZTAR = -ZTAR
	ZTAI = -ZTAI

L110:
	// AA=FACTOR FOR ANALYTIC CONTINUATION OF I(FNU,ZTA)
	// KODE=2 return S EXP(-ABS(XZTA))*I(FNU,ZTA) FROM CBESI
	ZTAR, ZTAI, FNU, KODE, _, CYR, CYI, NZ, RL, FNUL, TOL, ELIM, ALIM = ZBINU(ZTAR, ZTAI, FNU, KODE, 1, CYR, CYI, NZ, RL, FNUL, TOL, ELIM, ALIM)
	if NZ < 0 {
		goto L200
	}
	AA = FMR * FNU
	Z3R = SFAC
	STR = math.Cos(AA)
	STI = math.Sin(AA)
	S1R = (STR*CYR[1] - STI*CYI[1]) * Z3R
	S1I = (STR*CYI[1] + STI*CYR[1]) * Z3R
	FNU = (2.0e0 - FID) / 3.0e0
	ZTAR, ZTAI, FNU, KODE, _, CYR, CYI, NZ, RL, FNUL, TOL, ELIM, ALIM = ZBINU(ZTAR, ZTAI, FNU, KODE, 2, CYR, CYI, NZ, RL, FNUL, TOL, ELIM, ALIM)
	CYR[1] = CYR[1] * Z3R
	CYI[1] = CYI[1] * Z3R
	CYR[2] = CYR[2] * Z3R
	CYI[2] = CYI[2] * Z3R

	// BACKWARD RECUR ONE STEP FOR ORDERS -1/3 OR -2/3
	STR, STI = ZDIV(CYR[1], CYI[1], ZTAR, ZTAI)
	S2R = (FNU+FNU)*STR + CYR[2]
	S2I = (FNU+FNU)*STI + CYI[2]
	AA = FMR * (FNU - 1.0e0)
	STR = math.Cos(AA)
	STI = math.Sin(AA)
	S1R = COEF * (S1R + S2R*STR - S2I*STI)
	S1I = COEF * (S1I + S2R*STI + S2I*STR)
	if ID == 1 {
		goto L120
	}
	STR = CSQR*S1R - CSQI*S1I
	S1I = CSQR*S1I + CSQI*S1R
	S1R = STR
	BIR = S1R / SFAC
	BII = S1I / SFAC
	return BIR, BII, IERR
L120:
	STR = ZR*S1R - ZI*S1I
	S1I = ZR*S1I + ZI*S1R
	S1R = STR
	BIR = S1R / SFAC
	BII = S1I / SFAC
	return BIR, BII, IERR
L130:
	AA = C1*(1.0e0-FID) + FID*C2
	BIR = AA
	BII = 0.0e0
	return BIR, BII, IERR
L190:
	IERR = 2
	NZ = 0
	return BIR, BII, IERR
L200:
	if NZ == -1 {
		goto L190
	}
	NZ = 0
	IERR = 5
	return BIR, BII, IERR
L260:
	IERR = 4
	NZ = 0
	return BIR, BII, IERR
}
