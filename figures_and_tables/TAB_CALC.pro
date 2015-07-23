;Description
;-----------
;
;Determines the length of time needed to calculate the parameters
;of the 2CFM or 2CXM via LLS or NLLS method
;
;
;Syntax
;------
;
;TAB_CALC, nSim=nSim
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
;IDL> TAB_CALC, nSim=10^2D
;Calculation time for LLS (sec)        10.485781
;Calculation time for NLLS (min)        23.276203
;FoM       133.18723


;--------------------------------------------------------------------------------
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
;--------------------------------------------------------------------------------



PRO TAB_CALC, nSim=nSim

    if n_elements(nSim) eq 0 then nSim=1000.

    ;CALCULATE DATA

    ct = EXACT_CONC(Model='2CFM', Tacq=300.0, TIME=t, AIF=ca)
    msr_ct = MSR_CONC(NPIX=nSim, TIME=t, AIF=ca, CONC=ct[0,*], TR=1.25, CNR=50., MSR_T=msr_t, MSR_CA=msr_ca)

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