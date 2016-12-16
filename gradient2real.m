function [wx wy] = gradient2real(wxh,wyh)

wx = real(ifft2(wxh));
wy = real(ifft2(wyh));