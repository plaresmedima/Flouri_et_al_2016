;Description
;-----------
;
;Creates the plots of the error distribution at fixed CNR but
;variable TR as shown in Figure 4 of Flouri et al [submitted paper]
;
;
;Syntax
;------
;
;FIG_TR, nSim=nSim
;
;
;Arguments
;---------
;
;None
;
;
;Keywords
;--------
;
;nSim: the number of simulations to perform
;
;
;Example
;-------
;Plot the figure for fixed CNR but variable TR at 10^2D runs
;
;IDL> FIG_TR, nSim=10^2D



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
;------------------------------------------------------------------------------



PRO FIG_TR, nSim=nSim

    Device, Decomposed = 0
    DEVICE, Bypass_Translation=0

    if n_elements(nSim) eq 0 then nSim=10.

    TRmin = 5.0   ;sec
    TRmax = 10.0
    nTR = 10.

    dTR = (TRmax-TRmin)/(nTR-1.)
    TR = TRmin + dTR*findgen(nTR)

    ct = EXACT_CONC(Model='2CFM', Tacq=300.0, TIME=t, AIF=ca)

    ;GENERATE RESULTS
    Percent = fltarr(3,nTR,5,3) ;(Percentiles) x (TR's) x (parameters) x (methods)
    Tissue = fix(5*randomu(seed, nSim))
    FOR k=0L, nTR-1 DO BEGIN

        Error = fltarr(5,3,nSim) ;(5 parameters) x (3 methods)
        FOR i=0L, nSim-1 DO BEGIN
            p_exact = PARS(tissue[i])
            msr_ct = MSR_CONC(NPIX=1, TIME=t, AIF=ca, CONC=ct[tissue[i],*], TR=TR[k], CNR=1000.0, MSR_T=msr_t, MSR_CA=msr_ca)
            Error[0:3,0,i] = 100*(LLS_2CFM(msr_t-msr_t[0], msr_ct, msr_ca, FIT=fit_LLS)-p_exact)/p_exact
            Error[0:3,1,i] = 100*(LLS_2CFM(msr_t-msr_t[0], msr_ct, msr_ca, FIT=fit_WLLS, WEIGHTS=msr_ca)-p_exact)/p_exact
            Error[0:3,2,i] = 100*(NLLS_FIT(msr_t-msr_t[0], msr_ct, msr_ca, FIT=fit_NLLS, PARS(0)/2, Model='2CFM')-p_exact)/p_exact
            Error[4,0,i] = 100*norm(fit_LLS-msr_ct)/norm(msr_ct)
            Error[4,1,i] = 100*norm(fit_WLLS-msr_ct)/norm(msr_ct)
            Error[4,2,i] = 100*norm(fit_NLLS-msr_ct)/norm(msr_ct)
        ENDFOR

        FOR method=0,2 DO BEGIN
            FOR par=0,4 DO BEGIN
                Percent[*,k,par,method] = PERC(Error[par,method,*], [5., 50., 95.])
            ENDFOR
        ENDFOR
    ENDFOR


    ;DISPLAY RESULTS
    Parameters = ['FP', 'TP', 'PS', 'TE', 'Fit']
    Methods = ['LLS', 'WLLS', 'NLLS']

    window, 1, xsize=300, ysize=200, xpos=0, ypos=0
    window, 2, xsize=1500, ysize=600, xpos=0, ypos=0
    pos = 0
    FOR method=0,2 DO BEGIN
        FOR par=0,4 DO BEGIN

            WSET, 1
            PLOT, [min(TR)-0.5,max(TR)+0.5], [min(Percent[0,*,par,*]),max(Percent[2,*,par,*])], $
                /nodata, color=0, background=255, /YSTYLE, /XSTYLE, $
                XTITLE = 'TR (sec)', YTITLE = 'Error (%)', $
                TITLE=Parameters[par] + ' (' + Methods[method] + ')'
            OPLOT, TR, Percent[1,*,par,method], color=0, psym=4
            OPLOT, TR, TR*0, color=0, linestyle=1
            ERRPLOT, TR, Percent[0,*,par,method], Percent[2,*,par,method], color=0
            img = tvrd()

            WSET, 2
            TVSCL, img, pos
            pos = pos+1
        ENDFOR
    ENDFOR
    WDELETE, 1
    write_tiff, SIM_PATH() + 'TR_Error.tif', reverse(tvrd(),2)

END