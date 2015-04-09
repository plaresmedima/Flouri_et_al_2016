;Description
;-----------
;
;Calculates the pseudo-continuous concentration-time curve
;for a given tissue (0-4) and model (2CXM or 2CFM))
;following the formula in Eq. 6 of [submitted paper]
;and using a population-average input function (see also AIF())
;
;
;Syntax
;------
;
;concentrations = EXACT_CONC(Tissue, Model, TIME=t, AIF=ca)
;
;
;Arguments
;---------
;
;tissue: number from 0 to 4 defining the tissue type to be simulated (see also PARS())

;Model: string, either '2CXM' or '2CFM'
;

;Kewords
;-------
;
;TIME: named variable which upon return contains the pseudo-continuous times
;AIF: named variable which upon return contains the input function
;
;Example
;-------
;
;Plot 2CFM concentrations for patient 3
;
;IDL> c = EXACT_CONC(3, '2CFM', TIME=t, AIF=ca)
;IDL> plot, t, ca
;IDL> oplot, t, c
;

;---------------------------------------------------------------------------
;    Copyright (C) 2014 *Anonymised*
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


FUNCTION EXACT_CONC, Model=Model, Tacq=Tacq, TIME=t, AIF=ca

    dt = 0.1 ;sec
    nt = 1 + floor(Tacq/dt)
    t = dt*dindgen(nt) ;Pseudo-continuous time points (sec)
    ca = AIF(t) ;Arterial Input Function (mM)
    ct = fltarr(5,nt) ;Concentrations for 5 tissue types

    FOR tissue=0,4 DO BEGIN
        CASE model OF
          '2CFM': CONC_2CFM, [t,ca], PARS(tissue), conc
          '2CXM': CONC_2CXM, [t,ca], PARS(tissue), conc
        ENDCASE
        ct[tissue,*] = conc
    ENDFOR

    RETURN, ct
END