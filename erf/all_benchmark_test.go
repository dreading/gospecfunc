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

func BenchmarkErfcx(b *testing.B) {
	var r complex128
	for n := 0; n < b.N; n++ {
		r = Erfcx(complex(0.9,0))
	}
	GlobalC = r
}

func BenchmarkErfcxLargeNegative(b *testing.B) {
	var r complex128
	for n := 0; n < b.N; n++ {
		r = Erfcx(complex(-4,0))
	}
	GlobalC = r
}

func BenchmarkErfcxLargePositive(b *testing.B) {
	var r complex128
	for n := 0; n < b.N; n++ {
		r = Erfcx(complex(99,0))
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