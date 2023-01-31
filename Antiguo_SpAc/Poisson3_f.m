function rhs = Poisson3_f(xx,yy)
    omegax = 0.23; omegay = 0.49; omegapx = 0.331;omegapy = 0.667;
    rhs =  -pi*pi*(omegax*omegax + omegay*omegay)*sin(omegax*pi*xx + omegay*pi*yy)...
    -pi*pi*(omegapx*omegapx + omegapy*omegapy)*cos(omegapx*pi*xx + omegapy*pi*yy);
return