;Description
;-----------
;
;Simulates measurement with a given
;uniform sampling interval TR (sec)
;and Contrast-to-Noise Ratio (CNR).
;
;For more details see *Anonymous* et al (submitted)
;
;
;Syntax
;------
;
;msr_ct = MSR_CONC(t, ct, ca, TR, CNR, MSR_T=msr_t, MSR_CA=msr_ca)
;
;
;Returns
;-------
;
;msr_ct: the measured tissue concentrations
;
;
;Arguments
;---------
;
;t: pseudo-continuous times
;ct: pseudo-continuous tissue concentrations
;ca: pseudo-continuous arterial concentrations
;TR: Temporal resolution of the measurement (same units as t)
;CNR: Contrast-to-Noise ratio of the measurement
;    If CNR=0 then no noise is added
;
;
;Kewords
;-------
;
;MSR_T: named variable which upon return contains the time points of the measurement
;MSR_CA: named variable which upon return contains the measured arterial concentrations
;
;
;Example
;-------
;
;Measure 2CFM concentrations for patient 3 with TR=2s and CNR=50
;
;IDL> ct = EXACT_CONC(3, '2CFM', TIME=t, AIF=ca)
;IDL> msr_ct = MSR_CONC(t, ct, ca, 2.0, 50, MSR_T=msr_t, MSR_CA=msr_ca)
;IDL> plot, t, ct
;IDL> oplot, msr_t, msr_ct, psym=4
;IDL> plot, t, ca
;IDL> oplot, msr_t, msr_ca, psym=4
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


FUNCTION MSR_CONC, NPIX=npix, TIME=t, CONC=ct, AIF=ca, TR=TR, CNR=CNR, MSR_T=msr_t, MSR_CA=msr_ca

    nt = floor(max(t)/TR)
    t0 = TR*randomu(seed1)
    msr_t = t0 + TR*findgen(nt)

    msr_ca = INTERPOL(ca, t, msr_t)
    msr_ct_pix = INTERPOL(ct, t, msr_t)

    msr_ct = fltarr(npix, nt)
    for i=0L, npix-1 do msr_ct[i,*] = msr_ct_pix

    IF CNR EQ 0 THEN RETURN, msr_ct

    SD = max(ca)/CNR

    msr_ca = msr_ca + SD*RANDOMN(seed2, nt)
    msr_ct = msr_ct + SD*RANDOMN(seed3, nt*npix)

    RETURN, msr_ct

END