R_M = [74.6, 87.9, 64.2, 105.4, 78.1, 119.2, 151.8, 201.3, 110.6, 164.5, 158.5, 120.8]
A_N_min = [9.61, 3.22, 12.98, 3.74, 33.96, 29.96, 69.35, 16.45, 50.38, 60.05, 29.23, 14.88, 0.71]
;
n_win = 0
start_plot, n_win, 'test'
   plot, A_N_min[0:10], R_M[1:11], psym = 1
end_plot
end
