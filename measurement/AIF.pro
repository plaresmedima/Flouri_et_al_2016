;Description
;-----------
;
;Returns the experimentally derived
;population-averaged AIF as measured by
;Parker et al.(2006) Magn Reson Med 56: 993-1000
;-- extended with a baseline of 20s
;
;
;Syntax
;------
;
;ca = AIF(t)
;
;
;Arguments
;--------
;
;t: an array of time points, in seconds
;   where the function values are to be calculated
;
;
;Returns
;-------
;
;ca: an arry of the same length as t,
;    containing the plasma concentration in mM
;    at the corresponding times
;
;
;Example
;-------
;
;Plot the AIF for 60s in steps of 2s
;
;IDL> tacq = 60.0 & dt = 2.0
;IDL> nt = 1+floor(tacq/dt)
;IDL> t = findgen(nt)*dt
;IDL> plot, t, AIF(t)
;

;---------------------------------------------------------------------------
;    Copyright (C) 2010 *Anonimised*
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
;---------------------------------------------------------------------------



function AIF, t

;   parameter values defined in table 1 (Parker 2006)
	A1 = 0.809
	A2 = 0.330
	T1 = 0.17046
	T2 = 0.365
	sigma1 = 0.0563
	sigma2 = 0.132
	alpha = 1.050
	beta = 0.1685
	s = 38.078
	tau = 0.483

    ;Eq. 1 (Parker 2006)

	AIF = exp(-(t/60.-T1)^2/(2*sigma1^2))*A1/(sigma1*sqrt(2*!PI)) $
	    + exp(-(t/60.-T2)^2/(2*sigma2^2))*A2/(sigma2*sqrt(2*!PI)) $
	    + alpha*exp(-beta*t/60.)/(1+exp(-s*(t/60.-tau)))

;   baseline in sec
    t0=20.0

;   baseline shift
	n0 = total(t LT t0)
	if n0 gt 0 then begin
	    AIF = Shift(AIF,n0)
	    AIF[0:n0-1]=0
	endif

	return, AIF
end