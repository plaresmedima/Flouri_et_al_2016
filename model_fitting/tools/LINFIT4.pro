;Description
;-----------
;
;Calculates the values [a,b,c,d] for the coefficients
;in Eq.[16] of Flouri et al [submitted paper]
;
;
;Syntax
;------
;
;coeffs = LINFIT4(v, x, y, z, t, FIT=fit, CHISQ=chisq)
;
;Arguments
;---------
;
;v: system of equations to be solved
;x, y, z, t: n-element arrays
;
;
;Keywords
;--------
;
;FIT: optional, named variable which returns the best fit to the measured concentrations
;CHISQ: optional, named variable which returns the value of the summed, squared residuals
;       for the returned parameter values
;
;
;Returns
;-------
;
;coefficients: A 4-element array which contains
;              the values [a,b,c,d] for the coefficients of Eq.[16]
;Example
;-------
;
;IDL> x = findgen(100)/99.
;IDL> y = x^3
;IDL> z = -sqrt(x)
;IDL> t = exp(-x)
;IDL> v = 3*x + 4*y + 5*z + 6*t
;IDL> print, LINFIT4(v, x, y, z, t, FIT=vfit1)
;   2.99995      3.99999      4.99994      5.99998
;IDL> print, LINFIT4(v+0.1, x, y, z, t, FIT=vfit2)
;   3.05959      3.98666      4.98302      6.09845
;IDL> print, LINFIT4(v+1, x, y, z, t, FIT=vfit3)
;   3.59933      3.86565      4.83290      6.98502
;IDL> plot, v, linestyle=2
;IDL> oplot, vfit3
;IDL> oplot, vfit2
;IDL> oplot, vfit1
;

;----------------------------------------------------------------------------------------
;    Copyright (C) 2015, Dimitra Flouri and Steven Sourbron
;
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.
;
;    You should have received a copy of the GNU General Public License along
;    with this program; if not, write to the Free Software Foundation, Inc.,
;    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
;-----------------------------------------------------------------------------------------


FUNCTION LINFIT4, v, x, y, z, t, FIT=fit, CHISQ=chisq

    ;v = ax + by + cz + dt

    xx = TOTAL(x*x)
    yx = TOTAL(y*x)
    zx = TOTAL(z*x)
    tx = TOTAL(t*x)
    vx = TOTAL(v*x)

    yy = TOTAL(y*y)
    zy = TOTAL(z*y)
    ty = TOTAL(t*y)
    vy = TOTAL(v*y)

    zz = TOTAL(z*z)
    tz = TOTAL(t*z)
    vz = TOTAL(v*z)

    tt = TOTAL(t*t)
    vt = TOTAL(v*t)

    m = [[xx,yx,zx,tx, vx], $
         [yx,yy,zy,ty, vy], $
         [zx,zy,zz,tz, vz], $
         [tx,ty,tz,tt, vt]]

    a = m[1:*,0]/m[0,0]
    m = m[1:*,1:*] - m[0,1:*] ## a

    b = m[1:*,0]/m[0,0]
    m = m[1:*,1:*] - m[0,1:*] ## b

    c = m[1:*,0]/m[0,0]
    m = m[1:*,1:*] - m[0,1:*] ## c

    d = m[1:*,0]/m[0,0]

    d = + d[0]
    c = - c[0]*d + c[1]
    b = - b[0]*c - b[1]*d + b[2]
    a = - a[0]*b - a[1]*c - a[2]*d + a[3]

    IF ARG_PRESENT(fit) THEN fit = a*x + b*y + c*z + d*t
    IF ARG_PRESENT(chisq) THEN chisq = TOTAL((v - a*x - b*y - c*z - d*t)^2)

    RETURN, [a,b,c,d]

END