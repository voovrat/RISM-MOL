function r_cs = rHNC_closure(r,rGamma_in, EXPONENTA, u_short)
% !! size(r) = NxNumAtom*2

r_cs =  r.*exp(rGamma_in./r).*EXPONENTA - r - rGamma_in;