clear

exact_sol =1;
%exact_sol = 0;mx = fliplr([100 200 400 800 1600]/2); hx= 1./mx;
for m = [8 16 32 64 128 256 512]
%m = 128;
%m = 100; % number of grid points
h = 1/m; % Grid size
x = linspace(0,1-h,m)';
y = linspace(0,1-h,m)';
[xx yy] = meshgrid(x,y);
% x-coordinate

k = 2*pi*1i*[0:(m/2-1) (-m/2):(-1)]; % Vector of wavenumbers
[KX KY]  = meshgrid(k,k); % Matrix of (x,y) wavenumbers corresponding
                          % to Fourier mode (m,n)
%KX(1,1)=1;KY(1,1)=1;
L = 1;
c = (sin(2*pi*xx/L).*sin(2*pi*yy/L)).^100; %tracer initial condition.
%c = sin(2*pi*xx/L).*sin(2*pi*yy/L);

u = ones(m);v=u;

%uh = fft(u);vh = fft(v);
delsq = (KX.^2 + KY.^2);
poisson = delsq; poisson(1,1) = 4*pi^2;

dealias=KX/i<2/3*m&KY/i<2/3*m; % Cutting of frequencies using the 2/3 rule
%dealias = ones(m);
ch = fft2(c);

%test diffusion
dc = only_diffusion(delsq,ch);
c0 = 4*(-200*pi^2*c +9900*pi^2* (sin(2*pi*xx).*sin(2*pi*yy)).^98.*(((sin(2*pi*xx/L).*cos(2*pi*yy/L)).^2+(cos(2*pi*xx/L).*sin(2*pi*yy/L)).^2)));
%c0 = -8*pi^2*c;
c1 = real(ifft2(dc));

ach = compute_adv(KX,KY,dealias,u,v,ch);
ac = real(ifft2(ach));
a0 = (200*pi*(sin(2*pi*xx).*sin(2*pi*yy)).^99.*(sin(2*pi*yy).*cos(2*pi*xx)+cos(2*pi*yy).*sin(2*pi*xx)));

%eval(['w' num2str(log2(m/50)+1) '= w;'])
eval(['w_err(' num2str(log2(m/8)+1) ')= compute_error(h^2,c1,c0);']);
eval(['v_err(' num2str(log2(m/8)+1) ')= compute_error(h^2,ac,a0);']);

ach = compute_adv(KX,KY,ones(m),u,v,ch);
ac = real(ifft2(ach));

eval(['u_err(' num2str(log2(m/8)+1) ')= compute_error(h^2,ac,a0);']);

end

%compute_norm(wnorm1,wnorm2,wnorm3,wnorm4,wnorm5)

mx=[8 16 32 64 128 256 512];

paper_x = 21; paper_y = 21;
fig_x = 21*.8; fig_y = 21*.8; fig_lf = (paper_x-fig_x)/2; fig_bt = (paper_y-fig_y)/2;
axis_x=0.75; axis_y=0.8; axis_lf=0.15; axis_bt=(1-axis_y)/2;
resol = get(0,'screensize'); rat=ceil(0.8*resol(4)/fig_y);
scr_x = fig_x*rat; scr_y= fig_y*rat; scr_lf=(resol(3)-scr_x)/2; scr_bt=(resol(4)-scr_y)/2;

f=figure('visible','off');
set(gcf,'PaperUnits','centimeters');
set(gcf,'Position',[scr_lf scr_bt scr_x scr_y],'PaperPosition',[fig_lf fig_bt fig_x fig_y],'resize','on')
axes('position',[axis_lf axis_bt axis_x axis_y]);
fn='Times new roman'; xl='\itx \rm(\itkm\rm)'; yl='\ity \rm(\itkm\rm)';

subplot('position',[.1 .1 .8 .8]) 

loglog(mx,w_err,'-o',mx,u_err,'--o',mx,v_err,'o-.','linewidth',2)
legend('Diffusion','Advection','Advection w/ Filter','location','southwest')
axis([6 1024 1e-12 1e6]);grid on
ylabel('L_\infty error')
xlabel('number of grid')
title('spatial operator verification')
print(f,'-r200','-depsc','../fig/Verify_spatial')

close(f)

if 0
f=figure('visible','off');
set(gcf,'PaperUnits','centimeters');
set(gcf,'Position',[scr_lf scr_bt scr_x scr_y],'PaperPosition',[fig_lf fig_bt fig_x fig_y],'resize','on')
axes('position',[axis_lf axis_bt axis_x axis_y]);
fn='Times new roman'; xl='\itx \rm(\itkm\rm)'; yl='\ity \rm(\itkm\rm)';

subplot('position',[.1 .1 .8 .8]) 
[c h] = contour(xx,yy,u-u0);hold on
clabel(c,h);
fill([0.025 0.175 0.175 0.025 0.025],[.85 .85 .95 .95 .85],'w')
text(.05,.9,'\it{u_n-u_e}','fontsize',16)
set(gca,'xtick',[0:.1:1],'ytick',[0:.1:1])
title('u_n-u_e, m =64, max(u) = 0.75')
print(f,'-r200','-depsc',['../fig/UnsubUe' num2str(m)])

close(f)

f=figure('visible','off');
set(gcf,'PaperUnits','centimeters');
set(gcf,'Position',[scr_lf scr_bt scr_x scr_y],'PaperPosition',[fig_lf fig_bt fig_x fig_y],'resize','on')
axes('position',[axis_lf axis_bt axis_x axis_y]);
fn='Times new roman'; xl='\itx \rm(\itkm\rm)'; yl='\ity \rm(\itkm\rm)';

subplot('position',[.1 .1 .8 .8]) 
[c h] = contour(xx,yy,v-v0);hold on
clabel(c,h);
fill([0.025 0.175 0.175 0.025 0.025],[.85 .85 .95 .95 .85],'w')
text(.05,.9,'\it{v_n-v_e}','fontsize',16)
set(gca,'xtick',[0:.1:1],'ytick',[0:.1:1])
title('v_n-v_e, m =64, max(v) = 0.75')
print(f,'-r200','-depsc',['../fig/VnsubVe' num2str(m)])

close(f)

end
return


subplot('position',[.1 .51 .39 .39]) 
w = real(ifft2(wh));
[c h] = contour(xx,yy,w,[-8:2:8]);hold on
clabel(c,h);
fill([0.025 0.125 0.125 0.025 0.025],[.85 .85 .95 .95 .85],'w');
text(.05,.9,'\omega','fontsize',16)
set(gca,'xtick',[0:.1:1],'ytick',[0:.1:1],'xticklabel',[])

subplot('position',[.51 .51 .39 .39]) 
[c h] = contour(xx,yy,u,[0:.2:2]);hold on
clabel(c,h);
fill([0.025 0.125 0.125 0.025 0.025],[.85 .85 .95 .95 .85],'w')
text(.05,.9,'\it{u}','fontsize',16)
set(gca,'xtick',[0:.1:1],'ytick',[0:.1:1],'xticklabel',[],'yticklabel',[])


subplot('position',[.1 .1 .39 .39]) 
[c h] = contour(xx,yy,v,[0:.2:2]);hold on
clabel(c,h);
fill([0.025 0.125 0.125 0.025 0.025],[.85 .85 .95 .95 .85],'w')
text(.05,.9,'\it{v}','fontsize',16)
set(gca,'xtick',[0:.1:1],'ytick',[0:.1:1])

subplot('position',[.51 .1 .39 .39]) 
c = real(ifft2(ch));
[c h] = contour(xx,yy,c,[0.1:.1:1]);hold on
clabel(c,h);
fill([0.025 0.125 0.125 0.025 0.025],[.85 .85 .95 .95 .85],'w')
text(.05,.9,'c','fontsize',16)
set(gca,'xtick',[0:.1:1],'ytick',[0:.1:1],'yticklabel',[])

print(f,'-r200','-depsc',['../fig/Peudo_Spectrum_' num2str(m)])

close(f)
