// Copyright 2019 Infin IT Pty Ltd. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package airy_test

import (
	. "github.com/dreading/gospecfunc/airy"
	"testing"
)

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

func TestAiInt(t *testing.T) {
	testCases := []struct {
		num, den, res float64
	}{
		{-12.0e0, 1.0e0, -0.75228838916610124300e0},
		{-11.0e0, 1.0e0, -0.57348350185854889466e0},
		{-10.0e0, 1.0e0, -0.76569840313421291743e0},
		{-19.0e0, 2.0e0, -0.65181015505382467421e0},
		{-9.0e0, 1.0e0, -0.55881974894471876922e0},
		{-13.0e0, 2.0e0, -0.56902352870716815309e0},
		{-4.0e0, 1.0e0, -0.47800749642926168100e0},
		{-1.0e0, 1.0e0, -0.46567398346706861416e0},
		{-1.0e0, 4.0e0, -0.96783140945618013679e-1},
		{-1.0e0, 1024.0e0, -0.34683049857035607494e-3},
		{1.0e0, 1024.0e0, 0.34658366917927930790e-3},
		{1.0e0, 128.0e0, 0.27657581846051227124e-2},
		{1.0e0, 2.0e0, 0.14595330491185717833e0},
		{1.0e0, 1.0e0, 0.23631734191710977960e0},
		{4.0e0, 1.0e0, 0.33289264538612212697e0},
		{9.0e0, 2.0e0, 0.33318759129779422976e0},
		{6.0e0, 1.0e0, 0.33332945170523851439e0},
		{8.0e0, 1.0e0, 0.33333331724248357420e0},
		{10.0e0, 1.0e0, 0.33333333329916901594e0},
		{12.0e0, 1.0e0, 0.33333333333329380187e0},
	}

	for _, tc := range testCases {
		ζ := AiInt(tc.num / tc.den)
		if close(ζ, tc.res) == false {
			t.Fatalf("AiInt(%v): expected %v, got %v", tc.num/tc.den, tc.res, ζ)
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
