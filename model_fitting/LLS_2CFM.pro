;Description
;-----------
;
;Calculates the parameters of the 2CFM via a
;linear least squares method.
;
;For more details see Flouri et al [submitted paper]
;
;
;Syntax
;------
;
;pars = LLS_2CFM(t, ct, ca, FIT=fit, WEIGHTS=w)
;
;
;Arguments
;---------
;
;t: measured time points
;ct: tissue concentrations at the measured time points
;ca: arterial concentrations at the measured time points
;
;
;Keywords
;--------
;
;FIT: optional, named variable which returns the best fit to the measured concentrations
;WEIGHTS: weigthing function which performs weighthed linear least squares method
;         if set to 1, then the linear least squares method will be performed
;
;For more details see Flouri et al [submitted paper]
;
;
;Returns
;-------
;
;pars: 4-element floating point array with the values [FP, TP, PS, TE]
;    if a physiological solution does not exist, the scalar value 0B is returned
;
;
;Example
;-------
;
;Reconstruct 2CFM parameters for patient 3
;
;IDL> ct = EXACT_CONC(Model='2CFM', Tacq=300.0, TIME=t, AIF=ca)
;IDL> ct = ct[3,*]
;IDL> print, 'Exact parameters: ', PARS(3)
;Exact parameters:      0.026134801       7.2700000    0.0022222222       117.00000
;IDL> print, 'Reconstruction: ', LLS_2CFM(t, ct, ca, FIT=LLS_FIT)
;Reconstruction:       0.026134907       7.2699548    0.0022222305       116.99965


;-------------------------------------------------------------------------------------
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
;-------------------------------------------------------------------------------------



FUNCTION LLS_2CFM,  t, ct, ca, FIT=fit, WEIGHTS=w

    IF n_elements(w) EQ 0 THEN w=1+0*ct

    ct1 = INT_TRAP(t, ct)
    ct2 = INT_TRAP(t, ct1)
    ca1 = INT_TRAP(t, ca)
    ca2 = INT_TRAP(t, ca1)

    X = LINFIT4(ct*w, -ct2*w, -ct1*w, ca1*w, ca2*w)

    ;Extract physical parameters

    FP = X[2]
    TT = X[3]/(X[0]*X[2])

    prod = X[0]
    sum = X[1]
    det = sum^2 - 4*prod
    IF det LE 0 THEN root=0 ELSE root = sqrt(det)
    TP = (sum-root)/(2*prod)
    TE = (sum+root)/(2*prod)

    PS = FP*(TT-TP)/TE

    params = [FP, TP, PS, TE]

    IF ARG_PRESENT(fit) THEN fit = - X[0]*ct2 - X[1]*ct1 + X[2]*ca1 + X[3]*ca2

    RETURN, params
END
