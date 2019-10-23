// Copyright 2019 Infin IT Pty Ltd. All rights reserved.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
//
// The original Fortran code are from
//  ALGORITHM 644, TRANSACTIONS ON MATHEMATICAL SOFTWARE,
//	VOL. 21, NO. 4, December, 1995, P.  388--393.
//  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES

package amos

import (
	"github.com/dreading/gospecfunc/machine"
	"math"
	"math/cmplx"
)

// DGAMLN computes the natural log of the gamma function
func DGAMLN(z float64) float64 {
	r, _ := math.Lgamma(z)
	return r
}

// ZABS computes the absolute value or magnitude of a complex number real and imag parts
func ZABS(ZR float64, ZI float64) float64 {
	return cmplx.Abs(complex(ZR, ZI))
}

// ZFLOATS splits a complex number into real and imag parts
func ZFLOATS(Z complex128) (float64, float64) {
	return real(Z), imag(Z)
}

// ZLOG returns complex logarithm of a complex number real and imag parts
func ZLOG(AR float64, AI float64) (float64, float64) {
	return ZFLOATS(cmplx.Log(complex(AR, AI)))
}

// ZMLT returns complex multipliTION of complex numberSs real and imag parts
func ZMLT(AR float64, AI float64, BR float64, BI float64) (float64, float64) {
	return ZFLOATS(complex(AR, AI) * complex(BR, BI))
}

// ZSQRT returns complex square root of a complex number real and imag parts
func ZSQRT(AR float64, AI float64) (float64, float64) {
	return ZFLOATS(cmplx.Sqrt(complex(AR, AI)))
}

// ZDIV returns complex divide of a complex number real and imag parts
func ZDIV(AR float64, AI float64, BR float64, BI float64) (float64, float64) {
	return ZFLOATS(complex(AR, AI) / complex(BR, BI))
}

// ZEXP returns complex exponential function of a complex number real and imag parts
func ZEXP(AR float64, AI float64) (float64, float64) {
	return ZFLOATS(cmplx.Exp(complex(AR, AI)))
}

// ABS returns the absolute value of a
func ABS(a int) int {
	if a >= 0 {
		return a
	}
	return -a
}

// MIN returns the minimum value of a and b
func MIN(a int, b int) int {
	if a < b {
		return a
	}
	return b
}

// MAX returns the maximum value of a and b
func MAX(a int, b int) int {
	if a > b {
		return a
	}
	return b
}

// ZSHCH computes the complex hyperbolic functions SINH(x+i*y) and COSH(x+i*y)
// returns real and imag parts of sinh and cosh
func ZSHCH(ZR float64, ZI float64) (float64, float64, float64, float64) {
	SH := math.Sinh(ZR)
	CH := math.Cosh(ZR)
	SN, CN := math.Sincos(ZI)
	return SH * CN, CH * SN, CH * CN, SH * SN
}

// ZS1S2 tests for a possible underflow resulting from the addition of the I and K functions in the analytic
// continuation formula where S1=K function and S2=I function.
// For KODE=1 the I and K functions are different orders of magnitude, but for KODE=2 they can be of the same order
// of magnitude and the maximum must be at least one precision above the underflow limit.
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

// ZUCHK takes a scaled quantity Y whose magnitude is greater than
//
//		EXP(-ALIM)=ASCLE=1000*D1MACH[1]/TOL
//
// and tests if the magnitude of the real or imaginary part would underflow when y is scaled
// by TOL to its proper value. Y is accepted if the underflow is at least one precision
// below the magnitude of the largest component; otherwise the phase angle does not have
// absolute accuracy and an underflow is assumed.
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

// ZACAI applies the analytic continuation formula
//
// K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
// MP=PI*MR*CMPLX(0.0,1.0)
//
// to continue the K function from the right half to the left half z plane for use with ZAIRY
// where FNU=1/3 or 2/3 and N=1.
// ZACAI is the same as ZACON with the parts for larger orders an recurrence removed.
// A recursive call to ZACON can result if ZACON is called from ZAIRY.
func ZACAI(ZR float64, ZI float64, FNU float64, KODE int, MR int, N int, YR []float64, YI []float64, NZ int, RL float64, TOL float64, ELIM float64, ALIM float64) (float64, float64, float64, int, int, int, []float64, []float64, int, float64, float64, float64, float64) {

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
	ASCLE = 1000.0e0 * machine.D1MACH[1] / TOL
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

// ZAIRY computes the Airy functions A(z) and DAI(z) for complex z
// For KODE=1, ZAIRY computes the complex airy function Ai(z) or
// its derivative dAi(z)/dz on ID=0 or ID=1 respectively.
// For KODE=2, a scaling option Exp(zta)*Ai(z) or Exp(zta)*dAi(z)/dz
// is provided to remove the exponential decay in -𝛑/3 < arg(z) < 𝛑/3
// and the exponential growth in pi/3 < abs(arg(z)) < pi where
// zta=(2/3)*z*csqrt(z).
//
// While the airy functions Ai(z) and dAi(z)/dz are analytic in
// the whole z plane, the corresponding scaled functions defined
// for KODE=2 have a cut along the negative real axis.
// Defintions and notation are found in the nbs handbook of
// mathematical functions (ref. 1).
//
// 	 	INPUT    ZR,ZI are double precision
//   	ZR,ZI  - Z=CMPLX(ZR,ZI)
//   	ID     - order of derivative, ID=0 OR ID=1
//   	KODE   - a parameter to indicate the scaling option
//            	KODE= 1 returns
//                     AI=AI(Z)                on ID=0 or
//                     AI=DAI(Z)/DZ            on ID=1
//                = 2  returns
//                     AI=CEXP(ZTA)*AI(Z)       on ID=0 or
//                     AI=CEXP(ZTA)*DAI(Z)/DZ   on ID=1 where
//                     ZTA=(2/3)*Z*CSQRT(Z)
//
// 		OUTPUT AIR,AII are double precision
//   	AIR,AII- complex answer depending on the choices for ID and KODE
//   	NZ     - underflow indicator
//            	NZ= 0   , normal return
//            	NZ= 1   , AI=CMPLX(0.0e0,0.0e0) due to underflow in
//                      	-𝛑//3 < ARG(Z) < 𝛑//3 on KODE=1
//   	IERR   - error flag
//            	IERR=0, normal return - computation completed
//            	IERR=1, input error   - no computation
//            	IERR=2, overflow      - no computation, real(ZTA) too large on KODE=1
//            	IERR=3, ABS(Z) large      - computation completed
//                    losses of signifcance by argument reduction produce less than half of machine accuracy
//            	IERR=4, CABS(Z) TOO LARGE  - no computation complete loss of accuracy by argument reduction
//            	IERR=5, ERROR              - no computation, algorithm termination condition not met
//
// AI and DAI are computed for ABS(z) > 1.0 from the K Bessel functions by
//
//    	AI(Z)=C*SQRT(Z)*K(1/3,ZTA) , DAI(Z)=-C*Z*K(2/3,ZTA)
//      C=1.0/(PI*SQRT(3.0)) ZTA=(2/3)*Z**(3/2)
//
// with the power series for abs(z) <= 1.0
func ZAIRY(ZR float64, ZI float64, ID int, KODE int) (float64, float64, int, int) {

	// In most complex variable computation, one must evaluate elementary functions.
	// when the magnitude of Z is large, losses of significance by argument reduction
	// occur. consequently, if the magnitude of ZETA=(2/3)*Z^1.5 exceeds U1=SQRT(0.5/UR),
	// then losses exceeding half precision are likely and an error flag IERR=3 is
	// triggered where UR=Max(D1MACH(4),1.0e-18) is double precision unit roundoff
	// limited to 18 digits precision.
	// Also, if the magnitude of zeta is larger than U2=0.5/UR, then all significance is lost
	// and IERR=4. In order to use the int function, zeta must be further restricted not to
	// exceed the largest integer, U3=I1MACH(9). Thus, the magnitude of ZETA must be restricted
	// by min(U2,U3). On 32 bit machines, U1,U2, and U3 are approximately 2.0e+3, 4.2e+6, 2.1e+9.
	// In singLE precision arithmetic and 1.3e+8, 1.8e+16, 2.1e+9 in double precision arithmetic
	// respectively. This makes U2 and U3 limiting in their respective arithmetics. This means that
	// the magnitude of Z cannot exceed 3.1e+4 in single and 2.1e+6 in double precision arithmetic.
	// This also means that one can expect to retain, in the worst cases on 32 bit machines,
	// no digits in single precision and only 7 digits in double precision arithmetic. Similar
	// considerations hold for other machines.
	//
	// The approximate relative error in the magnitude of a complex bessel function can be
	// expressed by p*10^s where p=max(unit roundoff,1.0e-18) is the nominal precision and
	// 10^s represents the increase in error due to argument reduction in the elementary functions.
	// Here, s=max(1,abs(log10(abs(z))), abs(log10(fnu))) approximately (i.e. s=max(1,abs(exponent of
	// abs(z),abs(exponent of fnu)) ). However, the phase angle may have only absolute accuracy.
	// This is most likely to occur when one component (in absolute value) is larger than the other by
	// several orders of magnitude. If one component is 10**k larger than the other, then one can expect
	// only max(abs(log10(p))-k,0) significant digits; or, stated another way, when k exceeds
	// the exponent of p, no significant digits remain in the smaller component. however, the phase
	// angle retains absolute accuracy because, in complex arithmetic with precision p, the smaller
	// component will not (as a rule) decrease below p times the magnitude of the larger component.
	// In these extreme cases, the principal phase angle is on the order of +p, -p, 𝛑/2-p,
	// or -𝛑/2+p.
	//
	// References
	// Handbook of Mathematical Functions by M. Abramowitz and I. A. Stegun,
	// NBS AMS Series 55, U.S. Dept. of Commerce, 1955.
	//
	// Computation of Bessel functions of complex argument and large order
	// by D. E. Amos, sand83-0643, May, 1983
	//
	// A subroutine package for Bessel functions of a complex argument and
	// nonnegative order by D. E. Amos, sand85-1018, May, 1985
	//
	// A portable package for Bessel functions of a complex argument and
	// nonnegative order by D. E. Amos, Trans. Math. Software, 1986

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
	TOL = math.Max(machine.D1MACH[4], 1.0e-18)
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

	K1 = machine.I1MACH[15]
	K2 = machine.I1MACH[16]
	R1M5 = machine.D1MACH[5]
	K = MIN(ABS(K1), ABS(K2))
	ELIM = 2.303e0 * (float64(float32(K))*R1M5 - 3.0e0)
	K1 = machine.I1MACH[14] - 1
	AA = R1M5 * float64(float32(K1))
	DIG = math.Min(AA, 18.0e0)
	AA = AA * 2.303e0
	ALIM = ELIM + math.Max(-AA, -41.45e0)
	RL = 1.2e0*DIG + 3.0e0
	ALAZ = math.Log(AZ)

	// TEST FOR PROPER RANGE
	AA = 0.5e0 / TOL
	BB = float64(float32(machine.I1MACH[9])) * 0.5e0
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
	AA = 1000 * machine.D1MACH[1]
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
	ARM = 1000 * machine.D1MACH[1]
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
	BRY[1] = 1000.0 * machine.D1MACH[1] / TOL
	BRY[2] = 1.0e0 / BRY[1]
	BRY[3] = machine.D1MACH[2]
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
	T1 = float64(machine.I1MACH[14] - 1)
	T1 = T1 * machine.D1MACH[5] * 3.321928094e0
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

// ZKSCL SET K FUNCTIONS TO ZERO ON UNDERFLOW, CONTINUE RECURRENCE ON SCALED FUNCTIONS
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

	SCLE = machine.D1MACH[1] / TOL
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
	ARM = 1000.0 * machine.D1MACH[1]
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
