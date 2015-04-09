;Description
;-----------
;
;Defines five whole-kidney tissues:
;+ One representing normal kidneys using parameter values
;  measured in healthy volunteers [Sourbron et al Invest Radiol 2008].
;+ Four pathological kidneys taken from a
;  recent patient study on renal artery stenosis [Lim et al AJP 2013].
;  Patient cases were selected by identifying the kidneys
;  corresponding to the 10th and 90th percentiles in TE and VP
;All values are displayed in table 1 [submitted paper]
;
;
;Syntax
;------
;
;parameters = PARS(tissue)
;
;
;Arguments
;---------
;
;tissue: number from 0 to 4 where
;   0 = Healthy volunteer
;   1 = Patient 1
;   2 = Patient 2
;   3 = Patient 3
;   4 = Patient 4
;
;Returns
;-------
;
;parameters
;    A 4-element array which contains
;    the values [FP, V, E, xE] for the required tissue
;
;
;Example
;-------
;
;Display the parameters of Patient 1:
;
;IDL> print, PARS(1)
;    0.017894737      0.41000000      0.13148789      0.41463415


;---------------------------------------------------------------------------
;    Copyright (C) 2014 *Anonimised*
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
;---------------------------------------------------------------------------


FUNCTION PARS, tissue

    CASE tissue OF

    0: p = [6.50D, 125.D, 0.24D, 0.62D] ;Volunteer
    1: p = [9.50D, 102.D, 0.17D, 0.24D] ;Patient 1
    2: p = [13.9D, 153.D, 0.31D, 0.24D] ;Patient 2
    3: p = [7.27D, 117.D, 0.19D, 0.26D] ;Patient 3
    4: p = [10.3D, 214.D, 0.29D, 0.18D] ;Patient 4
    5: p = [100., 100., 0.29D, 0.18D] ;Testcase

    ENDCASE

    TP = p[0]
    TE = p[1]
    VP = p[2]
    VE = p[3]

    FP = VP/TP
    PS = VE/TE

    RETURN, [FP, TP, PS, TE]
END