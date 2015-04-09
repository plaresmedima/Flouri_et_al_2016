;Description
;-----------
;
;Calculates the parameters of the 2CFM via a
;non-linear least squares method.
;
;For more details see *Anonymous* et al (submitted)
;
;
;Syntax
;------
;
;pars = NLLS_2CFM(t, ct, ca, init, FIT=fit)
;
;
;Arguments
;---------
;
;t: measured time points
;ct: tissue concentrations at the measured time points
;ca: arterial concentrations at the measured time points
;init: 4-element array with initial values for [TP, TE, VP, VE]
;
;
;Keywords
;--------
;
;FIT: optional, named variable which returns the best fit to the measured concentrations
;
;
;Returns
;-------
;
;pars: 4-element floating point array with the values [TP, TE, VP, VE]
;
;
;Example
;-------
;
;Reconstruct 2CFM parameters for patient 3,
;with initial values derived from the healthy volunteer
;
;IDL> ct = EXACT_CONC(3, '2CFM', TIME=t, AIF=ca)
;IDL> print, 'Exact parameters: ', PARS(3)
;Exact parameters:        7.2700000       117.00000      0.19000000      0.26000000
;IDL> print, 'Reconstruction: ', NLLS_2CFM(t, ct, ca, PARS(0)/2)
;Reconstruction:        7.2700000       117.00000      0.19000000      0.26000000
;

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


FUNCTION NLLS_FIT, t, ct, ca, init, FIT=fit, MODEL=model

    params = init
	weights = 1 + 0*t ; no weighting
	fit = MPCURVEFIT([t, ca], ct, weights, params, FUNCTION_NAME='CONC_'+model, /NODERIVATIVE, /QUIET)
    return, params
END