// Copyright 2019 Infin IT Pty Ltd. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package bessel_test

import (
	. "github.com/dreading/gospecfunc/bessel"
	"testing"
)

func TestBesselI(t *testing.T) {
	testCases := []struct {
		α    float64
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{0, 5, 27.23987182360444689454423207588441928247906183222098381517},
		{0, 0.4 + 0.1i, 1.037751775187953907718589364295309545291505675072069744 + 0.02037712677480810300172139986765507486164198939561617240i},
		{0.25, 0.4 + 0.1i, 0.7636677889558093934973812584746228980011733000823678786 + 0.05894567689760461852351886060104048793278243522220770649i},
		{1, 0.4 + 0.1i, 0.2032605049177812708944870481129207354193884289297667737 + 0.05296680165856376304875068492944952721943227874651157176i},
		{0, 333 - 876i, -2.42466995791572885758160970181684018139584605311695e142 - 4.86259747781671556106159806889776193425074468393636e142i},
	}

	for _, tc := range testCases {
		ζ := I(tc.α, tc.x)
		if veryclose(real(ζ), real(tc.y)) == false {
			t.Fatalf("real(I(%v, %v)): expected %v, got %v", tc.α, tc.x, real(tc.y), real(ζ))
		}
		if veryclose(imag(ζ), imag(tc.y)) == false {
			t.Fatalf("imag(I(%v,%v)): expected %v, got %v", tc.α, tc.x, imag(tc.y), imag(ζ))
		}
	}
}

func TestBesselJ(t *testing.T) {
	testCases := []struct {
		α    float64
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{0, 1, 0.765197686557966551449717526102663220909274289755325241861},
		{0, 0.4 + 0.1i, 0.9627513455152844523412642515338186980389471552835888233 - 0.01962711629303271170774431370757696151493347606833454781i},
	}

	for _, tc := range testCases {
		ζ := J(tc.α, tc.x)
		if veryclose(real(ζ), real(tc.y)) == false {
			t.Fatalf("real(J(%v, %v)): expected %v, got %v", tc.α, tc.x, real(tc.y), real(ζ))
		}
		if veryclose(imag(ζ), imag(tc.y)) == false {
			t.Fatalf("imag(J(%v,%v)): expected %v, got %v", tc.α, tc.x, imag(tc.y), imag(ζ))
		}
	}
}

func TestBesselK(t *testing.T) {
	testCases := []struct {
		α    float64
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{0, 1, 0.421024438240708333335627379212609036136219748226660472298},
		{0, 0.4 + 0.1i, 1.0826035097235081563303065645729426351352961828491217647 - 0.21324459634740555852107854608689805665579172058695006929i},
	}

	for _, tc := range testCases {
		ζ := K(tc.α, tc.x)
		if veryclose(real(ζ), real(tc.y)) == false {
			t.Fatalf("real(K(%v, %v)): expected %v, got %v", tc.α, tc.x, real(tc.y), real(ζ))
		}
		if veryclose(imag(ζ), imag(tc.y)) == false {
			t.Fatalf("imag(K(%v,%v)): expected %v, got %v", tc.α, tc.x, imag(tc.y), imag(ζ))
		}
	}
}

func TestBesselY(t *testing.T) {
	testCases := []struct {
		α    float64
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{0, 1, 0.088256964215676957982926766023515162827817523090675546711},
		{0, 0.4 + 0.1i, -0.5873828735984329430094809243367886365339042810690287965 + 0.1750446663657107385528975496319911274395118490119615635i},
	}

	for _, tc := range testCases {
		ζ := Y(tc.α, tc.x)
		if close(real(ζ), real(tc.y)) == false {
			t.Fatalf("real(K(%v, %v)): expected %v, got %v", tc.α, tc.x, real(tc.y), real(ζ))
		}
		if close(imag(ζ), imag(tc.y)) == false {
			t.Fatalf("imag(K(%v,%v)): expected %v, got %v", tc.α, tc.x, imag(tc.y), imag(ζ))
		}
	}
}

func TestBesselH1(t *testing.T) {
	testCases := []struct {
		α    float64
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{0, 1, 0.7651976865579665514497175261026632209092742897553252418 + 0.08825696421567695798292676602351516282781752309067554671i},
	}

	for _, tc := range testCases {
		ζ := H1(tc.α, tc.x)
		if veryclose(real(ζ), real(tc.y)) == false {
			t.Fatalf("real(H(%v, %v)): expected %v, got %v", tc.α, tc.x, real(tc.y), real(ζ))
		}
		if veryclose(imag(ζ), imag(tc.y)) == false {
			t.Fatalf("imag(H(%v,%v)): expected %v, got %v", tc.α, tc.x, imag(tc.y), imag(ζ))
		}
	}
}

func TestBesselH2(t *testing.T) {
	testCases := []struct {
		α    float64
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{0, 1, 0.7651976865579665514497175261026632209092742897553252418 - 0.08825696421567695798292676602351516282781752309067554671i},
	}

	for _, tc := range testCases {
		ζ := H2(tc.α, tc.x)
		if veryclose(real(ζ), real(tc.y)) == false {
			t.Fatalf("real(H(%v, %v)): expected %v, got %v", tc.α, tc.x, real(tc.y), real(ζ))
		}
		if veryclose(imag(ζ), imag(tc.y)) == false {
			t.Fatalf("imag(H(%v,%v)): expected %v, got %v", tc.α, tc.x, imag(tc.y), imag(ζ))
		}
	}
}

func TestAiryAi(t *testing.T) {
	testCases := []struct {
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{5, 1.08344428136074417349865025033459804795777834796889391e-04},
		{5i, 29.9014823980071663793368297194577996107966955412306646481 + 21.6778315987836365558366270522814009526154445079208768325i},
		{-0.8 - 5i, -40.502439645308446470999918851421417821301812631645554975 - 117.58608694901610862936384033209723226387843162614560540i},
		{-80 - 50i, 1.28264717210434580202856048898366187480222379931511e196 + 4.86799742198478384026489758810414265077624185873133e195i},
		{69 + 27.345345i, -4.1680951604712214279525229780753604917334277179166e-158 - 3.1646991120150008455176064001093550612164673982528e-158i},
	}

	for _, tc := range testCases {
		ζ := Ai(tc.x)
		if soclose(real(ζ), real(tc.y), 1e-13) == false {
			t.Fatalf("real(AiryAi(%v)): expected %v, got %v", tc.x, real(tc.y), real(ζ))
		}
		if soclose(imag(ζ), imag(tc.y), 1e-13) == false {
			t.Fatalf("imag(AiryAi(%v)): expected %v, got %v", tc.x, imag(tc.y), imag(ζ))
		}
	}
}

func TestAiryBi(t *testing.T) {
	testCases := []struct {
		x, y complex128
	}{
		// extended precision values computed using Mathematica
		{5, 657.7920441711711824410805788744438785563124063292868570873},
	}

	for _, tc := range testCases {
		ζ := Bi(tc.x)
		if soclose(real(ζ), real(tc.y), 1e-13) == false {
			t.Fatalf("real(AiryBi(%v)): expected %v, got %v", tc.x, real(tc.y), real(ζ))
		}
		if soclose(imag(ζ), imag(tc.y), 1e-13) == false {
			t.Fatalf("imag(AiryBi(%v)): expected %v, got %v", tc.x, imag(tc.y), imag(ζ))
		}
	}
}

func TestGi(t *testing.T) {
	testCases := []struct {
		num, den, res float64
	}{
		{-1.0e0, 512.0e0, 0.20468308070040542435e0},
		{-1.0e0, 8.0e0, 0.18374662832557904078e0},
		{-1.0e0, 1.0e0, -0.11667221729601528265e0},
		{-4.0e0, 1.0e0, 0.31466934902729557596e0},
		{-8.0e0, 1.0e0, -0.37089040722426257729e0},
		{-33.0e0, 4.0e0, -0.25293059772424019694e0},
		{-9.0e0, 1.0e0, 0.28967410658692701936e0},
		{-10.0e0, 1.0e0, -0.34644836492634090590e0},
		{-11.0e0, 1.0e0, 0.28076035913873049496e0},
		{-13.0e0, 1.0e0, 0.21814994508094865815e0},
		{1.0e0, 512.0e0, 0.20526679000810503329e0},
		{1.0e0, 8.0e0, 0.22123695363784773258e0},
		{1.0e0, 1.0e0, 0.23521843981043793760e0},
		{4.0e0, 1.0e0, 0.82834303363768729338e-1},
		{7.0e0, 1.0e0, 0.45757385490989281893e-1},
		{29.0e0, 4.0e0, 0.44150012014605159922e-1},
		{8.0e0, 1.0e0, 0.39951133719508907541e-1},
		{9.0e0, 1.0e0, 0.35467706833949671483e-1},
		{10.0e0, 1.0e0, 0.31896005100679587981e-1},
		{12.0e0, 1.0e0, 0.26556892713512410405e-1},
	}

	for _, tc := range testCases {
		ζ := Gi(tc.num / tc.den)
		if close(ζ, tc.res) == false {
			t.Fatalf("Gi(%v): expected %v, got %v", tc.num/tc.den, tc.res, ζ)
		}

	}
}

func TestHi(t *testing.T) {
	testCases := []struct {
		num, den, res float64
	}{
		{-1.0e0, 512.0e0, 0.40936798278458884024e0},
		{-1.0e0, 8.0e0, 0.37495291608048868619e0},
		{-1.0e0, 1.0e0, 0.22066960679295989454e0},
		{-4.0e0, 1.0e0, 0.77565356679703713590e-1},
		{-8.0e0, 1.0e0, 0.39638826473124717315e-1},
		{-33.0e0, 4.0e0, 0.38450072575004151871e-1},
		{-9.0e0, 1.0e0, 0.35273216868317898556e-1},
		{-10.0e0, 1.0e0, 0.31768535282502272742e-1},
		{-11.0e0, 1.0e0, 0.28894408288051391369e-1},
		{-13.0e0, 1.0e0, 0.24463284011678541180e-1},
		{1.0e0, 512.0e0, 0.41053540139998941517e0},
		{1.0e0, 8.0e0, 0.44993502381204990817e0},
		{1.0e0, 1.0e0, 0.97220515514243332184e0},
		{4.0e0, 1.0e0, 0.83764237105104371193e2},
		{7.0e0, 1.0e0, 0.80327744952044756016e5},
		{29.0e0, 4.0e0, 0.15514138847749108298e6},
		{8.0e0, 1.0e0, 0.11995859641733262114e7},
		{9.0e0, 1.0e0, 0.21472868855967642259e8},
		{10.0e0, 1.0e0, 0.45564115351632913590e9},
		{12.0e0, 1.0e0, 0.32980722582904761929e12},
	}

	for _, tc := range testCases {
		ζ := Hi(tc.num / tc.den)
		if close(ζ, tc.res) == false {
			t.Fatalf("Hi(%v): expected %v, got %v", tc.num/tc.den, tc.res, ζ)
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
