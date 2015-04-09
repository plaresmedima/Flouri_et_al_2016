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
        msr_ct = MSR_CONC(NPIX=1, TIME=t, AIF=ca, CONC=ct[tissue[i],*], TR=2.0, CNR=500.0, MSR_T=msr_t, MSR_CA=msr_ca)  ;Data (mM)
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