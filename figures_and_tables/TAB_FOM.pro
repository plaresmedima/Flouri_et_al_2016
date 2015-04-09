PRO TAB_FOM, nSim=nSim, Protocol=Protocol

    if n_elements(nSim) eq 0 then nSim = 1000.  ;Set to 100000 for final results
    if n_elements(Protocol) eq 0 then Protocol = 1

    ct = EXACT_CONC(Model='2CFM', Tacq=300.0, TIME=t, AIF=ca)

    CASE Protocol OF
      1: P = {TR: 2.0,  CNR: 50.}
      2: P = {TR: 10.0, CNR: 500.}
      3: P = {TR: 2.0, CNR: 500.}
    ENDCASE


    Error = fltarr(4,3,nSim) ;(4 parameters) x (3 methods)
    Tissue = fix(5*randomu(seed, nSim))
    FOR i=0L, nSim-1 DO BEGIN
        p_exact = PARS(tissue[i])
        msr_ct = MSR_CONC(NPIX=1, TIME=t, AIF=ca, CONC=ct[tissue[i],*], TR=P.TR, CNR=P.CNR, MSR_T=msr_t, MSR_CA=msr_ca)
        Error[0:3,0,i] = 100*(LLS_2CFM(msr_t-msr_t[0], msr_ct, msr_ca, FIT=fit_LLS)-p_exact)/p_exact
        Error[0:3,1,i] = 100*(LLS_2CFM(msr_t-msr_t[0], msr_ct, msr_ca, FIT=fit_WLLS, WEIGHTS=msr_ca)-p_exact)/p_exact
        Error[0:3,2,i] = 100*(NLLS_FIT(msr_t-msr_t[0], msr_ct, msr_ca, FIT=fit_NLLS, PARS(0)/2, Model='2CFM')-p_exact)/p_exact
    ENDFOR



    Percent = fltarr(3,4,3) ;(Percentiles) x (4 parameters) x (3 methods)
    FOR method=0,2 DO BEGIN
        FOR par=0,3 DO BEGIN
            Percent[*,par,method] = PERC(Error[par,method,*], [5., 50., 95.])
        ENDFOR
    ENDFOR



    FOM = fltarr(2,4,2) ; (SE or RE) x (4 parameters) x (LLS or WLLS)
    FOR method=0,1 DO BEGIN
       FOM[0,*,method] = abs(Percent[1,*,2]) - abs(Percent[1,*,method])
       FOM[1,*,method] = (Percent[2,*,2]-Percent[0,*,2]) - (Percent[2,*,method]-Percent[0,*,method])
    ENDFOR


    ;EXPORT RESULTS
    FOM_string = strarr(5,6)
    FOM_string[1:*,0] = ['LLS',      'LLS',      'WLLS',     'WLLS']
    FOM_string[1:*,1] = ['Accuracy', 'Precision','Accuracy', 'Precision']
    FOM_string[0,2:*] = ['FP', 'TP', 'PS', 'TE']
    FOM_string[1:2,2:*] = strcompress(FOM[*,*,0],/remove_all)
    FOM_string[3:4,2:*] = strcompress(FOM[*,*,1],/remove_all)

    print, FOM_string
    write_csv, SIM_PATH() + 'FoM_'+ strcompress(regime,/remove_all) +'.csv', FoM_string

END