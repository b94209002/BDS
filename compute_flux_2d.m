function ux = compute_flux_2d(m,dx,u,c)
IR = [m 1:m-1];
IL = 1:m;
ux = -1/dx*u.*(c(IR,:)-c(IL,:));