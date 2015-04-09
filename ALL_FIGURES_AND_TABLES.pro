;Create all figures and tables in the paper

PRO ALL_FIGURES_AND_TABLES

    nSim = 10^3D

    FIG_FIT

    FIG_HIST, nSim=nSim
    FIG_TR, nSim=nSim
    FIG_CNR, nSim=nSim

    TAB_FOM, nSim=nSim, Protocol=1
    TAB_FOM, nSim=nSim, Protocol=2
    TAB_FOM, nSim=nSim, Protocol=3
    TAB_CALC, nSim=10^2D

END