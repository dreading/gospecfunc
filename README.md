[![Build Status](https://travis-ci.org/dreading/gospecfunc.svg?branch=master)](https://travis-ci.org/dreading/gospecfunc) 
[![Go Report Card](https://goreportcard.com/badge/github.com/dreading/gospecfunc)](https://goreportcard.com/report/github.com/dreading/gospecfunc)
[![stability-unstable](https://img.shields.io/badge/stability-unstable-yellow.svg)](https://github.com/emersion/stability-badges#unstable)

# Installation
The core packages of the gospecfunc suite are written in pure Go. Installation is done using go get.
```
go get -u github.com/dreading/gospecfunc
```

# Special Functions

## Bessel

Bessel and related functions:

Functions  | Domain | Description |
---------- | ------ | ----------- |
Ai     |  ℂ  | Airy Ai  function |
Aix    |  ℂ   | Exponentially scaled Airy Ai function|
Aid    |  ℂ   | First derivative of the Airy Ai function|
Aidx    |  ℂ   | Exponentially scaled first derivative of the Airy Ai function|
Bi     |  ℂ  | Biry Bi function |
Bix    |  ℂ   | Exponentially scaled Biry Bi function|
Bid    |  ℂ   | First derivative of the Biry Bi function|
Bidx    |  ℂ   | Exponentially scaled first derivative of the Biry Bi function|
Gi     |  ℝ  | Modified Airy Gi  function |
Hi     |  ℝ  | Modified Airy Hi  function |
I      | ℂ  | Modified Bessel function of the first kind  |
J      | ℂ  | Bessel function of the first kind |
K      | ℂ  | Modified Bessel function of the second kind  |
Y      | ℂ  | Bessel function of the second kind  |
H1      | ℂ  | Hankel fucntion of of the first kind  |
H2      | ℂ  | Hankel fucntion of of the second kind  |

## Erf

The error function and related functions:

Functions  | Domain | Description |
---------- | ------ | ----------- |
Erf    |  ℂ  | Error function |
Erfc    | ℂ  | Complementary error function  1 - Erf(ζ)
Erfcx    | ℂ  | Scaled complementary error function   exp(ζ²) Erfc(ζ) |
Erfi    |  ℂ  | Imaginary error function   -i erf(iζ) |
Dawson    |  ℂ  | Dawson's function. The one-sided Fourier–Laplace sine transform of the Gaussian function |
Fresnel |  ℂ  | Cos and Sin Fresnel integrals  |
Voigt |  ℝ  | Real and imaginary Voigt functions  𝖴(x,t) and 𝖵(x,t) |
Faddeyeva |  ℂ  | Plasma dispersion Faddeyeva function exp(-ζ²) Erfc(-iζ) |
  

## Integrals


# Testing 
```
 go test ./*/. 
```
# Benchmarking
```
 go test -bench=. ./*/.
```
# License
Original code is licensed under the Gonum License found in the LICENSE file. Portions of the code are subject to the additional licenses found in THIRD_PARTY_LICENSES.  
