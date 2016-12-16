function [wxh wyh] = gradient_omega(KX,KY,wh)

wxh = KX.*wh;
wyh = KY.*wh;

