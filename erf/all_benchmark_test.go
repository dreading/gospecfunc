// Copyright 2019 Infin IT Pty Ltd. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package erf_test

import (
	. "github.com/dreading/gospecfunc/erf"
	"testing"
)
  
// Global exported variables are used to store the
// return values of functions measured in the benchmarks.
// Storing the results in these variables prevents the compiler
// from completely optimizing the benchmarked functions away.
var (
	GlobalC complex128
	GlobalF float64
)

func BenchmarkErfcxReal(b *testing.B) {
	var r complex128
	for n := 0; n < b.N; n++ {
		r = Erfcx(complex(0.9,0))
	}
	GlobalC = r
}

func BenchmarkErfcxRealLargeNegative(b *testing.B) {
	var r complex128
	for n := 0; n < b.N; n++ {
		r = Erfcx(complex(-4,0))
	}
	GlobalC = r
}

func BenchmarkErfcxRealLargePositive(b *testing.B) {
	var r complex128
	for n := 0; n < b.N; n++ {
		r = Erfcx(complex(99,0))
	}
	GlobalC = r
}


func BenchmarkErfcxComplex(b *testing.B) {
	var r complex128
	for n := 0; n < b.N; n++ {
		r = Erfcx(complex(0.9,0.5))
	}
	GlobalC = r
}

func BenchmarkErfcxComplexLargeNegative(b *testing.B) {
	var r complex128
	for n := 0; n < b.N; n++ {
		r = Erfcx(complex(-4,-4))
	}
	GlobalC = r
}

func BenchmarkErfcxComplexLargePositive(b *testing.B) {
	var r complex128
	for n := 0; n < b.N; n++ {
		r = Erfcx(complex(99,123))
	}
	GlobalC = r
}

func BenchmarkErf(b *testing.B) {
	var r complex128
	for n := 0; n < b.N; n++ {
		r = Erf(0.9 + 0.4i)
	}
	GlobalC = r
}

func BenchmarkErfi(b *testing.B) {
	var r complex128
	for n := 0; n < b.N; n++ {
		r = Erfi(0.9 + 0.4i)
	}
	GlobalC = r
}

func BenchmarkErfc(b *testing.B) {
	var r complex128
	for n := 0; n < b.N; n++ {
		r = Erfc(0.9 + 0.4i)
	}
	GlobalC = r
}

func BenchmarkDawson(b *testing.B) {
	var r complex128
	for n := 0; n < b.N; n++ {
		r = Dawson(0.9 + 0.4i)
	}
	GlobalC = r
}

func BenchmarkFaddeyeva(b *testing.B) {
	var r complex128
	for n := 0; n < b.N; n++ {
		r = Faddeyeva(0.9 + 0.4i)
	}
	GlobalC = r
}

func BenchmarkFresnel(b *testing.B) {
	var c complex128  
	for n := 0; n < b.N; n++ {
		c,_ = Fresnel(0.9 + 0.4i)
	}
	GlobalC = c
} 


func BenchmarkVoigt(b *testing.B) {
	var r float64  
	for n := 0; n < b.N; n++ {
		r,_ = Voigt(0.9, 0.4)
	}
	GlobalF = r
} 

