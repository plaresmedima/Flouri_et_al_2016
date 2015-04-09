PRO TAB_CALC, nSim=nSim

    if n_elements(nSim) eq 0 then nSim=1000.

    ;CALCULATE DATA

    ct = EXACT_CONC(Model='2CFM', Tacq=300.0, TIME=t, AIF=ca)
    msr_ct = MSR_CONC(NPIX=nSim, TIME=t, AIF=ca, CONC=ct[0,*], TR=2.0, CNR=50., MSR_T=msr_t, MSR_CA=msr_ca)

    ;DETERMINE CALC TIME

    current_time = systime(1)
    FOR i=0L, nSim-1 DO fitpars = LLS_2CFM(msr_t, msr_ct[i,*], msr_ca)
    calc_time_LLS = (systime(1)-current_time)/nSim

    current_time = systime(1)
    FOR i=0L, nSim-1 DO fitpars = NLLS_FIT(msr_t, msr_ct[i,*], msr_ca, PARS(0)/2, Model='2CFM')
    calc_time_NLLS = (systime(1)-current_time)/nSim

    print, 'Calculation time for LLS (sec) ', calc_time_LLS*256.^2
    print, 'Calculation time for NLLS (min) ', calc_time_NLLS*256.^2/60.
    print, 'FoM', calc_time_NLLS/calc_time_LLS

END