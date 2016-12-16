function s = compute_slope(up,um,c)
s = up.*upwind_flux_1d(1,1,1,circshift(c,0))./upwind_flux_1d(1,1,1,circshift(c,-1));
s = s - um.*upwind_flux_1d(1,1,1,circshift(c,-2))./upwind_flux_1d(1,1,1,circshift(c,-1));