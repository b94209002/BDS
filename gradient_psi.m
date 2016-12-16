function [uh vh] = gradient_psi(KX,KY,psi)
 
uh = KY.*psi;
vh = -KX.*psi;