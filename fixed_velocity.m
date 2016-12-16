function [u v] = fixed_velocity(mx,my,u0,v0)

u = u0*ones(mx,my)';
v = v0*ones(mx,my)';