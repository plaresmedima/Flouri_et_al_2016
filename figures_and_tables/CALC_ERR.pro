;Description
;-----------
;
;Calculates the error for each reconstraction of a parameter
;and the goodness-of-fit as a percentage of the exact value,
;following the formulas in Eq.[22] and [23] respectively
;from Flouri et al [submitted paper]
;
;For more details see Figure 3 from [submitted paper]
;
;Syntax
;------
;
;CALC_ERR, nSim=nSim
;
;
;Arguments
;---------
;
;None
;
;
;Keywords
;-------
;
;nSim: the number of simulations to perform
;
;
;Example
;-------
;
;IDL> CALC_ERR, nSim=10^2D
;FP: LLS Median Error (%)    -0.497695
;FP: LLS 90% Confidence Interval (%)     -4.62087      3.79374
;FP: NLLS Median Error (%)      1.03936
;FP: NLLS 90% Confidence Interval (%)     -1.52252      3.91198


;-----------------------------------------------------------------------------
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
;-----------------------------------------------------------------------------



PRO CALC_ERR, nSim=nSim

    DEVICE, Decomposed = 0
    DEVICE, Bypass_Translation=0

    if n_elements(nSim) eq 0 then nSim=1000.

	ct = EXACT_CONC(Model='2CFM', Tacq=300.0, TIME=t, AIF=ca) ;Tissue concentration (mM)


    ;GENERATE RESULTS
    Error = fltarr(5,3,nSim) ;(5 parameters) x (3 methods)
    Tissue = fix(5*randomu(seed, nSim))
    FOR i=0L, nSim-1 DO BEGIN
        ;Measure
        p_exact = PARS(tissue[i])
        msr_ct = MSR_CONC(NPIX=1, TIME=t, AIF=ca, CONC=ct[tissue[i],*], TR=1.25, CNR=1000.0, MSR_T=msr_t, MSR_CA=msr_ca)  ;Data (mM)
       ;Fit
        Error[0:3,0,i] = 100*(LLS_2CFM(msr_t-msr_t[0], msr_ct, msr_ca, FIT=fit_LLS)-p_exact)/p_exact
        Error[0:3,1,i] = 100*(LLS_2CFM(msr_t-msr_t[0], msr_ct, msr_ca, FIT=fit_WLLS, WEIGHTS=msr_ca)-p_exact)/p_exact
        Error[0:3,2,i] = 100*(NLLS_FIT(msr_t-msr_t[0], msr_ct, msr_ca, FIT=fit_NLLS, PARS(0)/2, Model='2CFM')-p_exact)/p_exact
        Error[4,0,i] = 100*norm(fit_LLS-msr_ct)/norm(msr_ct)
        Error[4,1,i] = 100*norm(fit_WLLS-msr_ct)/norm(msr_ct)
        Error[4,2,i] = 100*norm(fit_NLLS-msr_ct)/norm(msr_ct)
    ENDFOR


    Percent = fltarr(3,5,3) ;(Percentiles) x (parameters) x (methods)
    FOR method=0,2 DO BEGIN
       FOR par=0,4 DO BEGIN
          Percent[*,par,method] = PERC(Error[par,method,*], [5., 50., 95.])
       ENDFOR
    ENDFOR

    print, 'FP: LLS Median Error (%)', Percent[1,0,0]
    print, 'FP: LLS 90% Confidence Interval (%)', Percent[0,0,0], Percent[2,0,0]

    print, 'FP: NLLS Median Error (%)', Percent[1,0,2]
    print, 'FP: NLLS 90% Confidence Interval (%)', Percent[0,0,2], Percent[2,0,2]

END