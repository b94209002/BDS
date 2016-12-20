function u = Crank_Nicolson(backward,forward,u,b)

% call compute diffusion
%diff = forward*reshape(u,m2,1);
diff = forward.*u;
%diff = Diffusion(length(b),4e-4,0.15,0.02,u);
%RHS = u +.5 * D + A
RHS =diff +b;
%RHS = diff + reshape(b,m2,1);
%u = u+ b;
%u = reshape(backward\RHS,m,m);
u = RHS./backward ;

