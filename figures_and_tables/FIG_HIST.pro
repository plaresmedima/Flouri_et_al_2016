PRO FIG_HIST, nSim=nSim

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
        msr_ct = MSR_CONC(NPIX=1, TIME=t, AIF=ca, CONC=ct[tissue[i],*], TR=2.0, CNR=500.0, MSR_T=msr_t, MSR_CA=msr_ca)  ;Data (mM)
       ;Fit
        Error[0:3,0,i] = 100*(LLS_2CFM(msr_t-msr_t[0], msr_ct, msr_ca, FIT=fit_LLS)-p_exact)/p_exact
        Error[0:3,1,i] = 100*(LLS_2CFM(msr_t-msr_t[0], msr_ct, msr_ca, FIT=fit_WLLS, WEIGHTS=msr_ca)-p_exact)/p_exact
        Error[0:3,2,i] = 100*(NLLS_FIT(msr_t-msr_t[0], msr_ct, msr_ca, FIT=fit_NLLS, PARS(0)/2, Model='2CFM')-p_exact)/p_exact
        Error[4,0,i] = 100*norm(fit_LLS-msr_ct)/norm(msr_ct)
        Error[4,1,i] = 100*norm(fit_WLLS-msr_ct)/norm(msr_ct)
        Error[4,2,i] = 100*norm(fit_NLLS-msr_ct)/norm(msr_ct)
    ENDFOR



    ;GENERATE HISTOGRAMS
    nbins = 5 + nSim/20.
    xhist = fltarr(5,3,nbins)
    yhist = fltarr(5,3,nbins)
    FOR par=0,4 DO BEGIN
  	  minb = PERC(Error[par,*,*],0.5)   ;Exclude outliers
	  maxb = PERC(Error[par,*,*],99.5)
	  binsize = (maxb-minb)/nbins
      FOR method=0,2 DO BEGIN
	     hist = histogram(Error[par,method,*],min=minb,max=maxb,nbins=nbins)
	     xhist[par,method,*] = binsize*findgen(nbins) + minb + binsize/2
	     yhist[par,method,*] = 100.0*hist/total(hist)
      ENDFOR
    ENDFOR



    ;DISPLAY HISTOGRAMS

    Parameters = ['FP', 'TP', 'PS', 'TE', 'Fit']
    Methods = ['LLS', 'WLLS', 'NLLS']

    window, 1, xsize=300, ysize=200, xpos=0, ypos=0
    window, 2, xsize=1500, ysize=600, xpos=0, ypos=0
    pos = 0
    FOR method=0,2 DO BEGIN
        FOR par=0,4 DO BEGIN

            WSET, 1
 	        plot $
 	        ,   [min(xhist[par,*,*]),max(xhist[par,*,*])] $
 	        ,   [0,max(yhist[par,*,*])]  $
 	        ,	Title = Parameters[par] + ' (' + Methods[method] +')' $
	        , 	/nodata, /xstyle, /ystyle, background=255, color=0 $
	        , 	xtitle = 'Error (%)', ytitle='Frequency (%)' $
	        , 	charsize=1.0, charthick=1.0, xthick=1.0, ythick=1.0, thick=1.0
	        oplot, xhist[par,method,*], yhist[par,method,*], psym=10, color=0
	        oplot, [0,0], [0,max(yhist[par,*,*])], color=0, thick=1, linestyle=2
            img = tvrd()

            WSET, 2
            TVSCL, img, pos
            pos = pos+1

        ENDFOR
    ENDFOR
    WDELETE, 1
    write_tiff, SIM_PATH() + 'Histograms.tif', reverse(tvrd(),2)

END