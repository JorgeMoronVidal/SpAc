function BC = GBC(N_fou,theta_l,x,y)
    n = (N_fou-1)/2;
    tt = atan2(y,x);
    if(tt <= 0)
        tt = 2*pi + tt;
    end
    BC = zeros(size(x));
    for k = -n:n
       BC = BC + cos(k*(theta_l-tt));
    end
    BC = -BC/N_fou;
return 