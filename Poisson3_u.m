function u = Poisson3_u(xx,yy)
    omegax = 0.24; omegay = 0.49; 
    omegapx = 0.331;omegapy = 0.667;
    u = sin(omegax*pi*xx+omegay*pi*yy) + cos(omegapx*pi*xx + omegapy*pi*yy);
return