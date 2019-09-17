package erf

import "github.com/dreading/gospecfunc/erf/internal/libcerf"

func Erfcx(x float64) float64 {
	return libcerf.Erfcx(x)
}
