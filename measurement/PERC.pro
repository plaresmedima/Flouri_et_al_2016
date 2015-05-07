;Description
;-----------
;
;Calculates the 5th, 50th and 95th percentile
;
;Syntax
;------
;
;Percent = PERC(X, P)
;
;Arguments
;---------
;
;X: n-element array of data values
;P: percentages for which to compute percentiles
;
;
;Keywords
;--------
;
;None
;
;
;Returns
;-------
;
;Percent: 3-element floating point array with the values [5., 50., 95.]
;
;
;Example
;-------
;
;Define the n-element array of sample data
;X = [65, 63, 67, 64, 68, 62, 70, 59, 68, 67, 69, 71, 66, 60]
;Calculate the 5th, 50th and 95th percentile
;IDL> print, 'result:',PERC(X, [0.5,50.,95.])
;result:      59      66      70


;------------------------------------------------------------------------------
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
;-------------------------------------------------------------------------------



FUNCTION PERC, X, P

	X = X[sort(X)]
	n = n_elements(X)
	return, X[floor((P/100D)*(n-1D))]
end

