FUNCTION PERC, X, P

	X = X[sort(X)]
	n = n_elements(X)
	return, X[floor((P/100D)*(n-1D))]
end