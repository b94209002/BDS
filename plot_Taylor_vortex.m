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

cinter = [-2:0.1:2];
ax = [0 1 0 1 ]; tick = [0:.2:1];cax= [-2 2];


subplot('position',[.1 .54 .4 .4]) 
contourf(x4,x4',w4,cinter);caxis(cax);axis(ax) 
set(gca,'xtick',tick,'xticklabel',[])
set(gca,'ytick',tick)
text(0.31,0.575,'m=512','fontsize',16)

subplot('position',[.55 .54 .4 .4])
contourf(x3,x3',w3,cinter);caxis(cax);axis(ax) 
set(gca,'xtick',tick,'xticklabel',[])
set(gca,'ytick',tick,'yticklabel',[])
text(0.31,0.575,'m=256','fontsize',16)

subplot('position',[.1 .12 .4 .4]) 
contourf(x2,x2',w2,cinter);caxis(cax);axis(ax) 
set(gca,'xtick',tick,'ytick',tick)
text(0.31,0.575,'m=128','fontsize',16)

subplot('position',[.55 .12 .4 .4]) 
contourf(x1,x1',w1,cinter);caxis(cax);axis(ax) 
set(gca,'xtick',tick)
set(gca,'ytick',tick,'yticklabel',[])
text(0.31,0.575,'m=64','fontsize',16)

h = colorbar('h');caxis([0 1])
set(h,'position',[.1 .06 .85 .02]);
print(f,'-r200','-depsc',['../fig/Steady_Taylor_vortex'])

close(f)

return
f=figure('visible','off');
set(gcf,'PaperUnits','centimeters');
set(gcf,'Position',[scr_lf scr_bt scr_x scr_y],'PaperPosition',[fig_lf fig_bt fig_x fig_y],'resize','on')
axes('position',[axis_lf axis_bt axis_x axis_y]);
fn='Times new roman'; xl='\itx \rm(\itkm\rm)'; yl='\ity \rm(\itkm\rm)';

cinter = [0.0:0.05:1];
ax = [0 1 0 1]; tick = [.1:.1:.9];
m=32;h = 1/m;
x = linspace(h/2,1-h/2,m)';y = linspace(h/2,1-h/2,m)';
[xx yy] = meshgrid(x,y);
[uc vc] = fixed_velocity2(xx,yy,xx,yy,a,a);
c = (sin(pi*x4/L).*sin(pi*x4'/L)).^100;
subplot('position',[.1 .10 .85 .85]) 
contourf(x4,x4',c,cinter);caxis([0 1]);hold on
quiver(xx,yy,uc,vc,'w');
axis(ax) 

set(gca,'xtick',tick)
set(gca,'ytick',tick)
%text(0.31,0.575,'m=512','fontsize',16)
title('Initial condition for varied velocity experiment')
h = colorbar('h');
set(h,'position',[.1 .03 .85 .02]);caxis([0,1])
print(f,'-r200','-depsc',['../fig/Varied_velocity_ini'])

close(f)
