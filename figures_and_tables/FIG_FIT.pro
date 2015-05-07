;Description
;-----------
;
;Creates the plots of the exact concentrations and model
;fits for protocol 1 at TR=2.0s and CNR=50 as shown in
;Figure 2 of Flouri et al [submitted paper]
;
;
;Syntax
;------
;
;FIG_FIT
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
;None
;
;
;Example
;-------
;
;FIG_FIT

;------------------------------------------------------------------------------
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
;-------------------------------------------------------------------------------



PRO FIG_FIT

    DEVICE, Decomposed = 0
    DEVICE, Bypass_Translation=0

	Tissue=0

    ;Generate data
    ct = EXACT_CONC(Model='2CFM', Tacq=300.0, TIME=t, AIF=ca) ;Tissue concentration (mM)
    msr_ct = MSR_CONC(NPIX=1, TIME=t, AIF=ca, CONC=ct[Tissue,*], TR=2.0, CNR=50.0, MSR_T=msr_t, MSR_CA=msr_ca)  ;Data (mM)

    ;Fit
    pars_LLS = LLS_2CFM(msr_t-msr_t[0], msr_ct, msr_ca, FIT=fit_LLS)
    pars_NLLS = NLLS_FIT(msr_t-msr_t[0], msr_ct, msr_ca, PARS(0)/2, FIT=fit_NLLS, Model='2CFM')

    ;Errors
    pars_exact = PARS(Tissue)
    FitErr_LLS = 100*norm(fit_LLS-msr_ct)/norm(msr_ct)
    Err_LLS = 100*(pars_LLS-pars_exact)/pars_exact
    FitErr_NLLS = 100*norm(fit_NLLS-msr_ct)/norm(msr_ct)
    Err_NLLS = 100*(pars_NLLS-pars_exact)/pars_exact
    FoM = abs([Err_NLLS,FitErr_NLLS]) - abs([Err_LLS,FitErr_LLS])

    WINDOW, 1, xsize=400, ysize=300, xpos=0, ypos=0
    WINDOW, 2, xsize=1200, ysize=300, xpos=0, ypos=0

    WSET, 1
    PLOT, [0,max(t)], [-0.5,0.5+max(ca)], $
        /nodata, color=0, background=255, /YSTYLE, /XSTYLE, $
        XTITLE = 'Time (sec)', YTITLE = 'Concentration (mM)', $
        TITLE = 'Arterial plasma concentration'
    OPLOT, t, t*0, color=0
    OPLOT, t, ca, color=0, linestyle=2
    OPLOT, msr_t, msr_ca, color=0, psym=4
    img = tvrd()
    WSET, 2
    TVSCL, img, 0

    WSET, 1
    PLOT, [0,max(t)], [-0.2,0.2+max(ct)], $
        /nodata, color=0, background=255, /YSTYLE, /XSTYLE, $
        XTITLE = 'Time (sec)', YTITLE = 'Concentration (mM)', $
        TITLE = 'Tissue concentration + LLS fit'
    OPLOT, t, t*0, color=0
    OPLOT, t, ct[Tissue,*], color=0, linestyle=2
    OPLOT, msr_t, msr_ct, color=0, psym=4
    OPLOT, msr_t, fit_LLS, color=0, thick=2
    img = tvrd()
    WSET, 2
    TVSCL, img, 1

    WSET, 1
    PLOT, [0,max(t)], [-0.2,0.2+max(ct)], $
        /nodata, color=0, background=255, /YSTYLE, /XSTYLE, $
        XTITLE = 'Time (sec)', YTITLE = 'Concentration (mM)', $
        TITLE = 'Tissue concentration + NLLS fit'
    OPLOT, t, t*0, color=0
    OPLOT, t, ct[Tissue,*], color=0, linestyle=2
    OPLOT, msr_t, msr_ct, color=0, psym=4
    OPLOT, msr_t, fit_NLLS, color=0, thick=2
    img = tvrd()
    WSET, 2
    TVSCL, img, 2

    WDELETE, 1

    WSET, 2
    x = 0.15
    y = 0.8
    dy = 0.05
    xyouts, x, y, 'Figures of Merit (Accuracy)', /normal, color=0 & y=y-dy
    xyouts, x, y, '-----------------', /normal, color=0 & y=y-dy
    str = strsplit(strcompress(FoM[0]),'.',/extract) & str=str[0]+'.'+strmid(str[1],0,2)
    xyouts, x, y, 'FP: ' + str + ' %', /normal, color=0  & y=y-dy
    str = strsplit(strcompress(FoM[1]),'.',/extract) & str=str[0]+'.'+strmid(str[1],0,2)
    xyouts, x, y, 'TP: ' + str + ' %', /normal, color=0  & y=y-dy
    str = strsplit(strcompress(FoM[2]),'.',/extract) & str=str[0]+'.'+strmid(str[1],0,2)
    xyouts, x, y, 'PS: ' + str + ' %', /normal, color=0  & y=y-dy
    str = strsplit(strcompress(FoM[3]),'.',/extract) & str=str[0]+'.'+strmid(str[1],0,2)
    xyouts, x, y, 'TE: ' + str + ' %', /normal, color=0  & y=y-dy
    str = strsplit(strcompress(FoM[4]),'.',/extract) & str=str[0]+'.'+strmid(str[1],0,2)
    xyouts, x, y, 'Fit: '+ str + ' %', /normal, color=0 & y=y-dy

    write_tiff, SIM_PATH() + 'Plots.tif', reverse(tvrd(),2)

END