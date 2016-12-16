clear

exact_sol =1;
%exact_sol = 0;mx = fliplr([100 200 400 800 1600]/2); hx= 1./mx;
for m = [64 128 256 512 1024]
%m = 100; % number of grid points
h = 1/m; % Grid size
x = linspace(h/2,1-h/2,m)';
xx = x + h/2;
% x-coordinate
%xi = linspace(h/2,1-h/2,m)';% x-coordinate on the interface
a = 1; % advection coiefficient
T = 1;
L = 1;
nu = 0.5; 
vavg = 1;
dt = nu*h;

u = sin(2*pi*xx);

c = sin(2*pi*x); %tracer initial condition.
c0 =c ;

for t = dt:dt:T

% advection tracer 
%c = Lax_wendroff_update_1d(dt,h,u,c);
c = Fromm_update_1d(dt,h,u,c);
%c = central_update_1d(dt,h,u,c);
%c = update_tracer_1d(dt,h,u,c);

%c0 = sin(2*pi*(x-t));
%
plot(x,c);title([num2str(t) ', max = ' num2str(max(c)) ]);drawnow
end

%eval(['w' num2str(log2(m/64)+1) '= c;'])
eval(['w_err(' num2str(log2(m/64)+1) ')= compute_error(h,c,c0);']);

end


%compute_norm(wnorm1,wnorm2,wnorm3,wnorm4,wnorm5)

return
