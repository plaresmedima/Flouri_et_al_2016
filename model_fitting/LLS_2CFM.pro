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
;pars = LLS_2CFM(t, ct, ca, FIT=fit)
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
;    if a physiological solution does not exist, the scalar value 0B is returned
;
;
;Example
;-------
;
;Reconstruct 2CFM parameters for patient 3
;
;IDL> ct = EXACT_CONC(3, '2CFM', TIME=t, AIF=ca)
;IDL> print, 'Exact parameters: ', PARS(3)
;Exact parameters:      0.026134801      0.45000000     0.085029240      0.57777778
;IDL> print, 'Reconstruction: ', LLS_2CFM(t, ct, ca)
;Reconstruction:      0.026134802      0.45000000     0.085029239      0.57777778

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




FUNCTION LLS_2CFM,  t, ct, ca, FIT=fit, WEIGHTS=w

    if n_elements(w) EQ 0 then w=1+0*ct

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
    if det LE 0 then begin
        TP = 2/sum
        TE = 2/sum
    endif else begin
        root = sqrt(det)
        TP = (sum-root)/(2*prod)
        TE = (sum+root)/(2*prod)
    endelse

    PS = FP*(TT-TP)/TE

    params = [FP, TP, PS, TE]

    IF ARG_PRESENT(fit) THEN fit = - X[0]*ct2 - X[1]*ct1 + X[2]*ca1 + X[3]*ca2

    RETURN, params
END
