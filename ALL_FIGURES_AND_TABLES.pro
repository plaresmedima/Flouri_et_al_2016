;Description
;-----------
;
;Create all figures and tables in the paper
;
;
;Syntax
;------
;
;ALL_FIGURES_AND_TABLES
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
;ALL_FIGURES_AND_TABLES


;-----------------------------------------------------------------------------------
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
;-----------------------------------------------------------------------------------



PRO ALL_FIGURES_AND_TABLES

    nSim = 10^4D

    FIG_FIT

    FIG_HIST, nSim=nSim
    FIG_TR, nSim=nSim
    FIG_CNR, nSim=nSim

    TAB_FOM, nSim=nSim, Protocol=1
    TAB_FOM, nSim=nSim, Protocol=2
    TAB_FOM, nSim=nSim, Protocol=3
    TAB_CALC, nSim=10^2D

END