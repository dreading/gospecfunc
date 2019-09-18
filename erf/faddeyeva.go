// The code below are based on Algorithm 916: Computing the Faddeyeva and Voigt Functions
// by Zaghloul, Mofreh R. and Ali, Ahmed N. in ACM Trans. Math. Softw. December 2011 Vol 38 No 2 Jan 2012
// http://doi.acm.org/10.1145/2049673.2049679. The go code is a port of the original matlab code.
package erf

import (
	"github.com/dreading/gospecfunc/erf/internal/libcerf"
	"math"
	"math/cmplx"
)

func Faddeyeva(z complex128) complex128 {
	// Faddeyeva or the plasma dispersion function, w(z) = exp(-z^2) * erfc(-i*z)
	// where z=x+iy and erfc(z) is the complex complementary error function of z;

	//For purely imaginary input values, use the built-in function erfcx.
	if IsReal(1i * z) {
		return complex(libcerf.Erfcx(imag(z)), 0)
	}

	var Rmin float64 = math.Nextafter(0, 1)
	var sqrt_log_Rmin float64 = math.Sqrt(-math.Log(Rmin))

	//Faddeyeva cannot calculate w for y negative & exp(y^2-x^2)>=the largest
	//floating point number due to unavoidable over flow problems
	if imag(z) < 0 && (imag(z)*imag(z)-real(z)*real(z)) >= sqrt_log_Rmin*sqrt_log_Rmin {
		return cmplx.Inf()
	}

	var eps float64 = 1e-16
	var tiny float64 = 0.06447 * eps
	var a float64 = math.Sqrt(-1 * math.Pi * math.Pi / math.Log(tiny/2))

	var half_a float64 = 0.5 * a
	var a_sqr float64 = a * a
	var two_a float64 = 2 * a
	var two_a_sqr float64 = 2 * a_sqr
	var four_a_sqr float64 = 2 * two_a_sqr
	var a_pi float64 = a / math.Pi
	var two_a_pi float64 = 2 * a_pi
	var one_sqrt_pi float64 = 1 / math.Sqrt(math.Pi)
	var erfcx_y float64 = libcerf.Erfcx(math.Abs(imag(z)))

	var x float64 = real(z)
	var y float64 = imag(z)
	var xsign float64 = Sign(x)
	var ysign float64 = Sign(y)
	x = xsign * x
	y = math.Max(Rmin, ysign*y)

	var x_sqr float64 = x * x
	var y_sqr float64 = y * y
	var two_yx float64 = 2 * y * x
	var two_a_x float64 = two_a * x
	var exp_x_sqr float64 = math.Exp(-x_sqr)
	var cos_2yx float64 = math.Cos(two_yx)
	var sin_2yx float64 = math.Sin(two_yx)

	var w complex128

	if x == 0 {
		//% For x=0, use the asymptotic expressions for x--->0 from Eq. (6) in the manuscript
		w = complex(erfcx_y, 0) * (complex(1, xsign*(x/y)))
	} else {
		var V_old float64 = exp_x_sqr * (erfcx_y*cos_2yx + (two_a_pi/y)*(math.Sin(two_yx/2)*math.Sin(two_yx/2)))
		var L_old float64 = -erfcx_y + a_pi/y
		var Sigma3 float64 = Rmin
		var Sigma5 float64 = Rmin
		var Sigma4_5 float64 = 0
		var delta3 float64 = 1
		var delta5 float64 = 1
		var n float64 = 0
		var n3 float64 = math.Ceil(x / a)
		var n3_3 = n3 - 1

		if (sqrt_log_Rmin - x) > 0 {
			var Sigma1 float64 = Rmin
			var Sigma2 float64 = Rmin
			var Sigma4 float64 = Rmin
			var exp1 float64 = math.Exp(-two_a_x)
			var exp2 float64 = math.Exp((four_a_sqr*n3 - 2*two_a_x) - two_a_sqr)
			var exp3 float64 = math.Exp(-((two_a_sqr*n3 - two_a_x) - two_a_sqr))
			var del2_tmp float64 = 1.0
			var del3_tmp float64 = math.Exp(-((a_sqr*n3*n3 - two_a_x*n3 - two_a_sqr*n3) + x_sqr + two_a_x + a_sqr))
			var del3_3_tmp float64 = math.Exp(a_sqr - (two_a_sqr*n3 - two_a_x))
			var del3 float64
			var del5 float64
			var den1 float64
			var exp_del1 float64
			var exp3_den float64
			var exp3_3_den float64

			for delta3 >= tiny || (delta5 >= tiny && n <= 50) {
				n = n + 1
				den1 = a_sqr*n*n + y_sqr
				exp_del1 = math.Exp(-(a_sqr * n * n)) / den1
				del2_tmp = del2_tmp * exp1
				if n3_3 >= 1 {
					del3_tmp = del3_tmp * exp3
					exp3_den = del3_tmp * exp_del1 * (den1 / (a_sqr*n3*n3 + y_sqr))
					del3_3_tmp = del3_3_tmp * exp2
					exp3_3_den = exp3_den * del3_3_tmp * ((a_sqr*n3*n3 + y_sqr) / (a_sqr*n3_3*n3_3 + y_sqr))
					del5 = (n3_3*exp3_3_den + n3*exp3_den)
					del3 = (exp3_3_den + exp3_den)
				} else {
					del3_tmp = del3_tmp * exp3
					del3 = del3_tmp * exp_del1 * (den1 / (a_sqr*n3*n3 + y_sqr))
					del5 = (n3 * del3)
				}
				delta3 = del3 / Sigma3
				delta5 = del5 / Sigma5

				Sigma1 = Sigma1 + exp_del1
				Sigma2 = Sigma2 + del2_tmp*exp_x_sqr*exp_del1
				Sigma3 = Sigma3 + del3
				Sigma4 = Sigma4 + n*del2_tmp*exp_x_sqr*exp_del1
				Sigma5 = Sigma5 + del5

				if x >= 5e-4 {
					Sigma4_5 = -Sigma4 + Sigma5
				} else {
					var two_a_xn_2 float64 = two_a_x * n * two_a_x * n
					Sigma4_5 = Sigma4_5 + 2*n*n*two_a_x*exp_x_sqr*exp_del1*(1+1.666666666666667e-001*(two_a_xn_2)+8.333333333333333e-003*(two_a_xn_2)*(two_a_xn_2))
				}
				n3 = n3 + 1
				n3_3 = n3_3 - 1
			}

			if y <= 5e0 && two_yx > Rmin {
				var w_real = V_old + y*two_a_pi*(-cos_2yx*exp_x_sqr*Sigma1+0.5*(Sigma2+Sigma3))
				var w_imag = xsign * (sin_2yx*exp_x_sqr*(L_old+two_a_pi*y*Sigma1) + two_a_pi*half_a*Sigma4_5)
				w = complex(w_real, w_imag)
			} else if y <= 5e0 && two_yx <= Rmin {
				var w_real = V_old + y*two_a_pi*(-cos_2yx*exp_x_sqr*Sigma1+.5*(Sigma2+Sigma3))
				var w_imag = xsign * (2*y*exp_x_sqr*(x*L_old+x*two_a_pi*y*Sigma1) + two_a_pi*half_a*Sigma4_5)
				w = complex(w_real, w_imag)
			} else {
				var w_real = V_old + y*two_a_pi*(-cos_2yx*exp_x_sqr*Sigma1+0.5*(Sigma2+Sigma3))
				var w_imag = xsign * (sin_2yx*exp_x_sqr*math.Min(0, math.Abs(L_old+(two_a_pi*y*Sigma1))) + two_a_pi*half_a*Sigma4_5)
				w = complex(w_real, w_imag)
			}

		} else if x >= sqrt_log_Rmin && x < 1e15 {
			var exp2 = math.Exp((four_a_sqr*n3 - 2*two_a_x) - two_a_sqr)
			var del3_3_tmp = math.Exp(a_sqr + (two_a_x - two_a_sqr*n3))

			for delta3 >= tiny || delta5 >= tiny && n <= 50 {
				n = n + 1
				var del3 float64
				var del5 float64
				if n3_3 >= 1 {
					var exp3_den = math.Exp(-(a*n3-x)*(a*n3-x)) / (a_sqr*n3*n3 + y_sqr)
					del3_3_tmp = del3_3_tmp * exp2
					var exp3_3_den = exp3_den * del3_3_tmp * ((a_sqr*n3*n3 + y_sqr) / (a_sqr*n3_3*n3_3 + y_sqr))
					del5 = (n3_3*exp3_3_den + n3*exp3_den)
					del3 = (exp3_3_den + exp3_den)
				} else {
					del3 = math.Exp(-(a*n3-x)*(a*n3-x)) / (a_sqr*n3*n3 + y_sqr)
					del5 = n3 * del3
				}
				delta3 = del3 / Sigma3
				delta5 = del5 / Sigma5
				Sigma3 = Sigma3 + del3
				Sigma5 = Sigma5 + del5
				n3 = n3 + 1
				n3_3 = n3_3 - 1
			}

			var w_real = V_old + y*a_pi*Sigma3
			var w_imag = xsign * (sin_2yx*exp_x_sqr*L_old + two_a_pi*half_a*Sigma5)
			w = complex(w_real, w_imag)

		} else {
			var w_real = one_sqrt_pi * (y / (x_sqr + y_sqr))
			var w_imag = one_sqrt_pi * (xsign * x / (x_sqr + y_sqr))
			w = complex(w_real, w_imag)
		}
	}

	if ysign < 0 {
		var two_exp_x_sqr_ysqr = 2 * math.Exp(-x_sqr+y_sqr)
		var w_real = two_exp_x_sqr_ysqr*cos_2yx - real(w)
		var w_imag = xsign*two_exp_x_sqr_ysqr*sin_2yx + imag(w)
		w = complex(w_real, w_imag)
	}

	return w
}

func IsReal(cmplx complex128) bool {
	return imag(cmplx) == 0
}

func Sign(x float64) float64 {
	switch {
	case math.IsNaN(x):
		return math.NaN()
	case math.IsInf(x, 1):
		return 1
	case math.IsInf(x, -1):
		return -1
	case x == 0:
		return 0
	case x < 0:
		return -1
	}
	return 1
}

/*

function w = Faddeyeva(z,tiny)
%----------
%   w = Faddeyeva(z,tiny) is the Faddeyeva or the plasma dispersion
%   function, w(z) = exp(-z^2) * erfc(-i*z) where z=x+iy and erfc(z) is the
%   complex  complementary error function of z; z is usually an array (with
%   one or two dimensions) but can be a single scalar as well.
%----------
%   The parameter "tiny" can be assigned values between the default value
%   tiny=1e-4 (for lowest accuracy & shortest run time) and tiny_min<eps
%   (for highest accuracy with longer run time) where eps is the floating
%   point relative accuracy in the computational platform. However, it is
%   not recommended for tiny_min to be << eps. A default value of tiny_min=
%   0.06447*eps is used herein and any smaller input value for
%   tiny will be automatically changed to this value with a warning
%   message.
%-------------
%   The accompanying driver code Faddeyeva_driver.m can be run for
%   computation of the Faddeyeva function and the partial derivatives
%   of its real part, V(x,y) on an array of the complex variable z.The
%   partial derivatives of the imaginary part, L(x,y) are given simply
%   by Eq. (23) in the manuscript.
%   An example of generating an array of z is included in the driver code.
%-------------
%   Authors: M. R. Zaghloul & A. N. Ali,
%   United Arab Emirates University, May 28,2011
%----------

% Machine-related parameters
Rmin=realmin;                     % Smallest positive floating point number
eps0=eps;                         % Floating point relative accuracy
sqrt_log_Rmin=sqrt(-log(Rmin));



a=realsqrt(-(pi*pi)/reallog(tiny/2));          % parameter a in Eq. (8)
tiny=max(tiny,eps0);

%-------------------

% Repeatedly used parameters & constants
half_a=.5*a;
a_sqr=a*a;
two_a=2*a;
two_a_sqr=2*a_sqr;
four_a_sqr=2*two_a_sqr;
a_pi=a/pi;
two_a_pi=2*a_pi;
one_sqrt_pi=(1/sqrt(pi));

% Reset outputs
sizein=size(z);
w=repmat(NaN,sizein);
idxmax=prod(sizein);
%--------------------

% Calculate erfcx(y).
erfcx_y=erfcx(abs(imag(z)));
%--------------------

x=real(z(idx));
y=imag(z(idx));
xsign=sign(x);
ysign=sign(y);
x=xsign*x;
y=max(Rmin,ysign*y);

x_sqr=x*x;
y_sqr=y*y;
two_yx=2*y*x;
two_a_x=two_a*x;
exp_x_sqr=exp(-x_sqr);
cos_2yx=cos(two_yx);
sin_2yx=sin(two_yx);

% For x=0, use the asymptotic expressions for x--->0 from Eq. (6) in
% the manuscript
if x==0
    w(idx)=erfcx_y(idx)*(1+1i*xsign*(x/y));
else

    V_old=exp_x_sqr*(erfcx_y(idx)*cos_2yx+(two_a_pi/y)*sin(two_yx/2)^2);
    L_old=(-erfcx_y(idx)+a_pi/y);

    % Initialization of the sums Sigma3, Sigma5 & Sigma4_5
    Sigma3=Rmin;
    Sigma5=Rmin;
    Sigma4_5=0;

    delta3=1;
    delta5=1;

    n=0;
    n3=ceil(x/a);
    n3_3=n3-1;

    if (sqrt_log_Rmin-x)>0

        % Initialization of the sums Sigma1, Sigma2 & Sigma4
        Sigma1=Rmin;
        Sigma2=Rmin;
        Sigma4=Rmin;

        exp1=exp(-two_a_x);
        exp2=exp((four_a_sqr*n3-2*two_a_x)-two_a_sqr);
        exp3=exp(-((two_a_sqr*n3-two_a_x)-two_a_sqr));
        del2_tmp=1.0;
        del3_tmp=exp(-((a_sqr*n3*n3-two_a_x*n3-two_a_sqr*n3)+x_sqr+two_a_x+a_sqr));
        del3_3_tmp=exp(a_sqr-(two_a_sqr*n3-two_a_x));

        while (delta3>=tiny || delta5>=tiny && n<=50)
            n=n+1;
            den1=a_sqr*n*n+y_sqr;
            exp_del1=exp(-(a_sqr*n*n))/den1;
            del2_tmp=del2_tmp*exp1;

            if n3_3>=1
                del3_tmp=del3_tmp*exp3;
                exp3_den=del3_tmp*exp_del1*(den1/(a_sqr*n3*n3+y_sqr));
                del3_3_tmp=del3_3_tmp*exp2;
                exp3_3_den=exp3_den*del3_3_tmp*((a_sqr*n3*n3+y_sqr)/(a_sqr*n3_3*n3_3+y_sqr));
                del5=(n3_3*exp3_3_den+n3*exp3_den);
                del3=(exp3_3_den+exp3_den);
            else
                del3_tmp=del3_tmp*exp3;
                del3=del3_tmp*exp_del1*(den1/(a_sqr*n3*n3+y_sqr));
                del5=(n3*del3);
            end

            delta3=del3/Sigma3;
            delta5=del5/Sigma5;

            Sigma1=Sigma1+exp_del1;
            Sigma2=Sigma2+del2_tmp*exp_x_sqr*exp_del1;
            Sigma3=Sigma3+del3;
            Sigma4=Sigma4+n*del2_tmp*exp_x_sqr*exp_del1;
            Sigma5=Sigma5+del5 ;

            if x>=5e-4
                Sigma4_5=-Sigma4+Sigma5;
            else
                Sigma4_5=Sigma4_5+2*n*n*two_a_x*exp_x_sqr*exp_del1*(1+1.666666666666667e-001*(two_a_x*n)^2+...
                    8.333333333333333e-003*(two_a_x*n)^4);
            end

            n3=n3+1;
            n3_3=n3_3-1;

        end

        if y <= 5e0 && two_yx>Rmin
            w(idx)=V_old+y*two_a_pi*(-cos_2yx*exp_x_sqr*Sigma1+.5*(Sigma2+Sigma3))+...
                1i* xsign*(sin_2yx*exp_x_sqr*(L_old+two_a_pi*y*Sigma1)+two_a_pi*half_a*Sigma4_5);
        elseif y <= 5e0 && two_yx<=Rmin
            w(idx)=V_old+y*two_a_pi*(-cos_2yx*exp_x_sqr*Sigma1+.5*(Sigma2+Sigma3))+...
                1i* xsign*(2*y*exp_x_sqr*(x*L_old+x*two_a_pi*y*Sigma1)+two_a_pi*half_a*Sigma4_5);
        else
            w(idx)=V_old+y*two_a_pi*(-cos_2yx*exp_x_sqr*Sigma1+.5*(Sigma2+Sigma3))+...
                1i* xsign*(sin_2yx*exp_x_sqr*min(0,abs(L_old+(two_a_pi*y*Sigma1)))+two_a_pi*half_a*Sigma4_5);
        end

    elseif x>=sqrt_log_Rmin && x<1e15
        exp2=exp((four_a_sqr*n3-2*two_a_x)-two_a_sqr);
        del3_3_tmp=exp(a_sqr+(two_a_x-two_a_sqr*n3));

        while (delta3>=tiny || delta5>=tiny && n<=50)
            n=n+1;

            if n3_3>=1
                exp3_den=exp(-(a*n3-x)*(a*n3-x))/(a_sqr*n3*n3+y_sqr);
                del3_3_tmp=del3_3_tmp*exp2;
                exp3_3_den=exp3_den*del3_3_tmp*((a_sqr*n3*n3+y_sqr)/(a_sqr*n3_3*n3_3+y_sqr));
                del5=(n3_3*exp3_3_den+n3*exp3_den);
                del3=(exp3_3_den+exp3_den);
            else
                del3=exp(-(a*n3-x)^2)/(a_sqr*n3*n3+y_sqr);
                del5=n3*del3;
            end

            delta3=del3/Sigma3;
            delta5=del5/Sigma5;
            Sigma3=Sigma3+del3 ;
            Sigma5=Sigma5+del5;
            n3=n3+1;
            n3_3=n3_3-1;

        end
        w(idx)=V_old+y*a_pi*Sigma3+1i*xsign*(sin_2yx*exp_x_sqr*L_old+two_a_pi*half_a*Sigma5);
    else
        w(idx)=one_sqrt_pi*((y+1i*xsign*x)/(x_sqr+y_sqr)) ;
    end
end

if ysign<0
    two_exp_x_sqr_ysqr=2*exp(-x_sqr+y_sqr);
    w(idx)=two_exp_x_sqr_ysqr*cos_2yx-real(w(idx))-1i*(-xsign*two_exp_x_sqr_ysqr*sin_2yx-imag(w(idx)));
end
*/
//(0.9849331541661086+0.05937031524864555i)
//   9849424862549036e-001 6.965909657459020e-002
