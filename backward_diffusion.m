function wh = backward_diffusion(m,dt,mu,delsq,ddw,N)

wh = (ddw+N)./(ones(m)/dt -.5*mu*delsq);