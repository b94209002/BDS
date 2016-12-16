clear

exact_sol =1;
%exact_sol = 0;mx = fliplr([100 200 400 800 1600]/2); hx= 1./mx;
for m = [64 128 256 512]
m
%m = 100; % number of grid points
h = 1/m; % Grid size
x = linspace(h/2,1-h/2,m)';y = linspace(h/2,1-h/2,m)';
[xx yy] = meshgrid(x,y);
x = linspace(0,1-h,m)';y = linspace(0,1-h,m)';
[xf yf] = meshgrid(x,y);

% x-coordinate
%xi = linspace(h/2,1-h/2,m)';% x-coordinate on the interface
a = 1; % advection coiefficient
T = 1;
mu = 0.05; % diffusion coefficient
L = 1;
nu = 0.5; 
dt = nu*h;

% u(1,1) is the velocity at (0,0)
[u v] = fixed_velocity(m,m,a,0);
%[u v] = fixed_velocity2(xx,yy,a,a);
[ufun vfun pxfun pyfun] = taylor_vortex_function(0,0,1);


%[c c2]=taylor_vortex(ufun,vfun,xx,yy,0);
%c = (sin(pi*xx/L).*sin(pi*yy/L)).^100; %tracer initial condition.
c = zeros(m);
c((xx' < .6 & xx' >.4)) = 1;

c0 =c;
%c0 =  (sin(pi*(xx-u')/L).*sin(pi*(yy-v')/L)).^10;
for t =dt:dt:T

%[u v] = taylor_vortex(ufun,vfun,xf,yf,t);
%[uc vc] = taylor_vortex(ufun,vfun,xx,yy,t);
%[uc1 vc1] = taylor_vortex(ufun,vfun,xx,yy,t-dt);
% advection tracer 
%c1 = Lax_wendroff_update_2d(dt,h,h,circshift(u,[-1 0]),circshift(v,[0 -1]),c1);
%c1 = Fromm_update_2d(dt,h,h,u,v,c1);
px = -pxfun(xx,yy,t);px=0;
c = BDS_update_2d(dt,h,h,u,v,px,c);
%c2 = BDS_update_2d(dt,h,h,u,v,c2);
%c0 = (sin(pi*(xx-t)/L).*sin(pi*(yy-t)/L)).^100;

%pcolor(xx,yy,c');shading flat;colorbar;%hold on;quiver(xx,yy,u',v','w');hold off;
%title([' t = ' num2str(t) ' max = ' num2str(max(max(c))) ', min = ' num2str(min(min(c))) ] );drawnow
end

%eval(['w' num2str(log2(m/64)+1) '= c;'])
eval(['w_err(' num2str(log2(m/64)+1) ')= compute_error(h^2,c,c0);']);

end


%compute_norm(wnorm1,wnorm2,wnorm3,wnorm4,wnorm5)

return
