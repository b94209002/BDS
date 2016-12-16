function ddw = forward_diffusion(m,dt,mu,delsq,wh)

ddw =( ones(m)/dt +.5*mu*delsq ).* wh;