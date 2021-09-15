function theta = Intersection_Square(Cs)
    theta = NaN(16,1);
    theta_h = acos(1/(sqrt(2)*Cs));
    theta_v = asin(1/(sqrt(2)*Cs));
    theta_t =  acos(1/Cs);
    theta(1) = theta_h;
    theta(2) = pi*0.25 - theta_t;
    theta(3) = pi*0.25 + theta_t;
    theta(4) = theta_v;
    theta(5) = pi - theta_v;
    theta(6) = pi*0.75 - theta_t;
    theta(7) = pi*0.75 + theta_t;
    theta(8) = pi -theta_h;
    theta(9) = pi + theta_h;
    theta(10) = pi*1.25 - theta_t;
    theta(11) = pi*1.25 + theta_t;
    theta(12) = pi + theta_v;
    theta(13) = 2*pi - theta_v;
    theta(14) = pi*1.75 - theta_t;
    theta(15) = pi*1.75 + theta_t;
    theta(16) = 2*pi -theta_h;
end