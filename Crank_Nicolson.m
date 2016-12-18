function u = Crank_Nicolson(backward,forward,u,b)

% call compute diffusion
diff = forward.*u;
%diff = Diffusion(length(b),4e-4,0.15,0.02,u);
% RHS = u +.5 * D + A
RHS = diff + b;
%u = u+ b;
u = backward.*RHS;

