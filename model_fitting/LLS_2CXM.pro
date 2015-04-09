;Description
;-----------
;
;Calculates the parameters of the 2CFM via a
;linear least squares method.
;
;For more details see *Anonymous* et al (submitted)
;
;
;Syntax
;------
;
;pars = LLS_2CXM(t, ct, ca, FIT=fit)
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
;
;
;Returns
;-------
;
;pars: 4-element floating point array with the values [TP, TE, VP, VE]
;
;
;Example
;-------
;
;Reconstruct 2CFM parameters for patient 3
;
;IDL> ct = EXACT_CONC(3, '2CXM', TIME=t, AIF=ca)
;IDL> print, 'Exact parameters: ', PARS(3)
;Exact parameters:        7.2700000       117.00000      0.19000000      0.26000000
;IDL> print, 'Reconstruction: ', LLS_2CXM(t, ct, ca)
;Reconstruction:        7.2699549       116.99965      0.18999959      0.26000018
;

;---------------------------------------------------------------------------
;    Copyright (C) 2014 *Anonimised*
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

FUNCTION LLS_2CXM,  t, ct, ca, FIT=fit, WEIGHTS=w

    if n_elements(w) EQ 0 then w=1+0*ct

    ct1 = INT_TRAP(t, ct)
    ct2 = INT_TRAP(t, ct1)
    ca1 = INT_TRAP(t, ca)
    ca2 = INT_TRAP(t, ca1)

    X = LINFIT4(ct*w, -ct2*w, -ct1*w, ca1*w, ca2*w)

    ;Extract physical parameters

    FP = X[2]
    TT = X[3]/(X[0]*X[2])
    TE = X[1]/X[0] - TT
    TP = 1/(X[0]*TE)
    PS = FP*(TT-TP)/TE

    params = [FP, TP, PS, TE]

   IF ARG_PRESENT(fit) THEN fit = - X[0]*ct2 - X[1]*ct1 + X[2]*ca1 + X[3]*ca2

    RETURN, params

END