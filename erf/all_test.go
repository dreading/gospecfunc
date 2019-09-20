// Copyright 2019 Infin IT Pty Ltd. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package erf_test

import (
	. "github.com/dreading/gospecfunc/erf"
	"testing"
)

func TestErf(t *testing.T) {
	testCases := []struct {
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{0.5 + 0.5i, 0.64261291485482052831942135847199581467891706129770675079 + 0.45788139443519221584208890063522927648352418879876980026i},
	}

	for _, tc := range testCases {
		y := Erf(tc.x)
		if veryclose(real(y), real(tc.y)) == false {
			t.Fatalf("real(Feaddeyeva(%v)): expected %v, got %v", tc.x, real(tc.y), real(y))
		}
		if veryclose(imag(y), imag(tc.y)) == false {
			t.Fatalf("imag(Feaddeyeva(%v)): expected %v, got %v", tc.x, imag(tc.y), imag(y))
		}
	}
}

func TestErfi(t *testing.T) {
	testCases := []struct {
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{0.5 + 0.5i, 0.45788139443519221584208890063522927648352418879876980026 + 0.64261291485482052831942135847199581467891706129770675079i},
	}

	for _, tc := range testCases {
		y := Erfi(tc.x)
		if veryclose(real(y), real(tc.y)) == false {
			t.Fatalf("real(Feaddeyeva(%v)): expected %v, got %v", tc.x, real(tc.y), real(y))
		}
		if veryclose(imag(y), imag(tc.y)) == false {
			t.Fatalf("imag(Feaddeyeva(%v)): expected %v, got %v", tc.x, imag(tc.y), imag(y))
		}
	}
}

func TestErfc(t *testing.T) {
	testCases := []struct {
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{0.5 + 0.5i, 0.35738708514517947168057864152800418532108293870229324920 - 0.45788139443519221584208890063522927648352418879876980026i},
	}

	for _, tc := range testCases {
		y := Erfc(tc.x)
		if veryclose(real(y), real(tc.y)) == false {
			t.Fatalf("real(Feaddeyeva(%v)): expected %v, got %v", tc.x, real(tc.y), real(y))
		}
		if veryclose(imag(y), imag(tc.y)) == false {
			t.Fatalf("imag(Feaddeyeva(%v)): expected %v, got %v", tc.x, imag(tc.y), imag(y))
		}
	}
}

func TestDawson(t *testing.T) {
	testCases := []struct {
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{0.5 + 0.5i, 0.62914469771362783370956774548491213744194694694642085904 + 0.30523946561753882091323070230974006125244937735067903082i},
	}

	for _, tc := range testCases {
		y := Dawson(tc.x)
		if veryclose(real(y), real(tc.y)) == false {
			t.Fatalf("real(Feaddeyeva(%v)): expected %v, got %v", tc.x, real(tc.y), real(y))
		}
		if veryclose(imag(y), imag(tc.y)) == false {
			t.Fatalf("imag(Feaddeyeva(%v)): expected %v, got %v", tc.x, imag(tc.y), imag(y))
		}
	}
}

func TestFaddeyeva(t *testing.T) {
	testCases := []struct {
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{0.5 + 0.5i, 0.53315670791217491376822891204271112100489475433534083731 + 0.23048823138445840870767807113455955862989744369028872730i},
		{0.5 - 0.5i, 1.22200841586857051846433425316494818297839563988414726868 + 1.18933930859286440925425394156570233479350417957149007767i},
	}

	for _, tc := range testCases {
		y := Faddeyeva(tc.x)
		if veryclose(real(y), real(tc.y)) == false {
			t.Fatalf("real(Feaddeyeva(%v)): expected %v, got %v", tc.x, real(tc.y), real(y))
		}
		if veryclose(imag(y), imag(tc.y)) == false {
			t.Fatalf("imag(Feaddeyeva(%v)): expected %v, got %v", tc.x, imag(tc.y), imag(y))
		}
	}
}

func TestFrenselC(t *testing.T) {
	testCases := []struct {
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{0.5 + 0.5i, 0.53173595500717020917826034727737322043857351236142204171 + 0.53173595500717020917826034727737322043857351236142204171i},
	}

	for _, tc := range testCases {
		y, _ := Fresnel(tc.x)
		if veryclose(real(y), real(tc.y)) == false {
			t.Fatalf("real(FresnelC(%v)): expected %v, got %v", tc.x, real(tc.y), real(y))
		}
		if veryclose(imag(y), imag(tc.y)) == false {
			t.Fatalf("imag(FresnelC(%v)): expected %v, got %v", tc.x, imag(tc.y), imag(y))
		}
	}
}

func TestFrenselS(t *testing.T) {
	testCases := []struct {
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{0.5 + 0.5i, -0.1367816577291388468638554200956170098637129887108742800 + 0.1367816577291388468638554200956170098637129887108742800i},
	}

	for _, tc := range testCases {
		_, y := Fresnel(tc.x)
		if veryclose(real(y), real(tc.y)) == false {
			t.Fatalf("real(FresnelS(%v)): expected %v, got %v", tc.x, real(tc.y), real(y))
		}
		if veryclose(imag(y), imag(tc.y)) == false {
			t.Fatalf("imag(FresnelS(%v)): expected %v, got %v", tc.x, imag(tc.y), imag(y))
		}
	}
}

func TestVoigt(t *testing.T) {
	testCases := []struct {
		x, t, vr, vi float64
	}{
		{0.5, 0.5, 0.6181707462326899, 0.16454229282575275},
		{1, 1e-34, 0.5, 0.5},
		{0, 0, 1, 0},
	}

	for _, tc := range testCases {
		vr, vi := Voigt(tc.x, tc.t)
		if veryclose(vr, tc.vr) == false {
			t.Fatalf("Voigt Real(%v, %v)): expected %v, got %v", tc.x, tc.t, tc.vr, vr)
		}
		if veryclose(vi, tc.vi) == false {
			t.Fatalf("Voigt Imag(%v, %v)): expected %v, got %v", tc.x, tc.t, tc.vi, vi)
		}
	}
}

func TestErfcxLargeNegative(t *testing.T) {

	testCases := []struct {
		x, y float64
	}{
		// extended precision valeus computed using Mathematica
		{-24.6438055098077, 1.1363627998302537408199038021823478437489991303892970e264},
		{-18.3606202026281, 5.0941755016494381919909749370853183247108146946632678e+146},
		{-9.72124040525614, 2.2026768805395279137105540037081856761793654217012477e+41},
		{-6.57964775166634, 1.2659147073374019834469401791473064878069367955720208e+19},
		{-5.7942495882689, 7.6162503596169089643842718139903640837346133123846494e14},
		{-3.43805509807655, 271948.8992883505292945176350662819469190070229287071247260},
	}

	for _, tc := range testCases {
		if y := Erfcx(complex(tc.x, 0)); soclose(real(y), tc.y, 1e-13) == false {
			t.Fatalf("Erfcx(%v): expected %v, got %v", tc.x, tc.y, y)
		}
	}
}

func TestErfcx(t *testing.T) {

	testCases := []struct {
		x, y float64
	}{
		// extended precision valeus computed using Mathematica
		{-2.1, 164.2938078850581786033535756010373886889407300324149511047},
		{-2, 108.9409043899779724123554338248132140422788747719728953862},
		{-1.86725877128166, 65.08259822571413624464960190629088108079118255820727492512},
		{-1.08186060788421, 6.040496719420559069042337844922667141609211183706177952338},
		{-0.296462444486757, 1.446697524293911950679978694884888005694830774823915998041},
		{0.488935718910689, 0.621407170543139014621201513672065153966191496240882991744},
		{1.27433388230814, 0.362803991242338377291710601242918924966030832510719387918},
		{2.05973204570559, 0.249162426872206342991325437388674726308559083733372555711},
		{2.84513020910303, 0.187820125319619076448476381651127670262889920476862189097},
		{99.4491043069892, 0.005672862203514460617907910565771973977540996637069157954},
	}

	for _, tc := range testCases {
		if y := Erfcx(complex(tc.x, 0)); veryclose(real(y), tc.y) == false {
			t.Fatalf("Erfcx(%v): expected %v, got %v", tc.x, tc.y, y)
		}
	}
}

// The floating point comparison tests are copied from from math/all_test.go.
func tolerance(a, b, e float64) bool {
	// Multiplying by e here can underflow denormal values to zero.
	// Check a==b so that at least if a and b are small and identical
	// we say they match.
	if a == b {
		return true
	}
	d := a - b
	if d < 0 {
		d = -d
	}
	// note: b is correct (expected) value, a is actual value.
	// make error tolerance a fraction of b, not a.
	if b != 0 {
		e = e * b
		if e < 0 {
			e = -e
		}
	}
	return d < e
}

func close(a, b float64) bool      { return tolerance(a, b, 1e-14) }
func veryclose(a, b float64) bool  { return tolerance(a, b, 5e-16) }
func soclose(a, b, e float64) bool { return tolerance(a, b, e) }
