;Description
;-----------
;
;Integrates a vector using the trapezoidal rule.
;;
;
;Syntax
;------
;
;int = INT_TRAP(t, c)
;
;
;Arguments
;---------
;
;t: X-coordinates of vector values
;c: vector to be integrated (same nr of elements as t)
;

;Returns
;-------
;
;pars: floating point array with the same nr of elements as t
;
;
;Example
;-------
;
;Integrate concentrations for patient 3
;
;IDL> ct = EXACT_CONC(3, '2CFM', TIME=t, AIF=ca)
;IDL> ct_int = INT_TRAP(t, ct)
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

function INT_TRAP, t, c

    n = n_elements(t)
    dt = t[1:n-1]-t[0:n-2]

    return, [0, TOTAL( dt*(c[0:n-2]+c[1:n-1])/2.0 , /cumulative)]
end