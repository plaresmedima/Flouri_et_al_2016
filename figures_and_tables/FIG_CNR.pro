;Description
;-----------
;
;Creates the plots of the error distribution at fixed TR but
;variable CNR as shown in Figure 3 of Flouri et al [submitted paper]
;
;
;Syntax
;------
;
;FIG_CNR, nSim=nSim
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
;
;Plot the figure for fixed TR but variable CNR at 10^2D runs
;
;IDL> FIG_CNR, nSim=10^2D
;


;---------------------------------------------------------------------------------
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
;---------------------------------------------------------------------------------


PRO FIG_CNR, nSim=nSim

    Device, Decomposed = 0
    DEVICE, Bypass_Translation=0

    if n_elements(nSim) eq 0 then nSim=1000.

    NCRmin = 1./1000
    NCRmax = 1./50   ;High-noise limit
    nNCR = 10.0

    dNCR = (NCRmax-NCRmin)/(nNCR-1.)
    NCR = NCRmin + dNCR*findgen(nNCR)

    ct = EXACT_CONC(Model='2CFM', Tacq=300.0, TIME=t, AIF=ca)

    ;GENERATE RESULTS
    Percent = fltarr(3,nNCR,5,3) ;(Percentiles) x (CNR's) x (parameters) x (methods)
    Tissue = fix(5*randomu(seed, nSim))
    FOR k=0L, nNCR-1 DO BEGIN

        Error = fltarr(5,3,nSim) ;(5 parameters) x (3 methods)
        FOR i=0L, nSim-1 DO BEGIN
            p_exact = PARS(tissue[i])
            msr_ct = MSR_CONC(NPIX=1, TIME=t, AIF=ca, CONC=ct[tissue[i],*], TR=1.25, CNR=1/NCR[k], MSR_T=msr_t, MSR_CA=msr_ca)
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
            PLOT, [0,max(100*NCR)+0.1], [min(Percent[0,*,par,*]),max(Percent[2,*,par,*])], $
                /nodata, color=0, background=255, /YSTYLE, /XSTYLE, $
                XTITLE = '100/CNR', YTITLE = 'Error (%)', $
                TITLE=Parameters[par] + ' ('+Methods[method]+')'
            OPLOT, 100*NCR, Percent[1,*,par,method], color=0, psym=4
            OPLOT, 100*NCR, NCR*0, color=0, linestyle=1
            ERRPLOT, 100*NCR, Percent[0,*,par,method], Percent[2,*,par,method], color=0
            img = tvrd()

            WSET, 2
            TVSCL, img, pos
            pos = pos+1

        ENDFOR
    ENDFOR
    WDELETE, 1
    write_tiff, SIM_PATH() + 'CNR_Error.tif', reverse(tvrd(),2)

END