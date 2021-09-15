function [N,dtheta,dist] = Compute_Optimal_N(Cs,Nmin,Nmax)
    thetas = Intersection_Square(Cs);
    thetas = sort(thetas);
    figure 
    hold on 
    x = cos(2*pi*0.01*(1:101));
    y = sin(2*pi*0.01*(1:101));
    plot(x,y);
    scatter(cos(thetas),sin(thetas));
    N = Nmin:Nmax;
    dtheta = 2*pi*1./(N);
    dist = NaN(length(N),length(thetas));
    for i = 1:length(N)
        theta_test = thetas(1) + dtheta(i) *(0:N(i)-1);
        for j = 1:length(thetas)
            dist(i,j) = min(abs(theta_test-thetas(j)));
        end
    end
    figure 
    plot(N,sum(dist,2)'./dtheta)
    index = find(~((sum(dist,2)'./dtheta) - min((sum(dist,2)'./dtheta))));
    aux_theta = thetas(1) - 0.5*dtheta(index) + dtheta(index) *(0:N(index)-1);
    theta_0 = min(aux_theta(aux_theta>=0));
    figure 
    hold on
    title("N = "+string(N(index)) +" d\theta = "+string(dtheta(index)) + " \theta_0 = " + string(theta_0))
    plot(x,y,'blue');
    plot(x + sqrt(2)/Cs,y,'green');
    plot(x + sqrt(2)/Cs,y + sqrt(2)/Cs,'green');
    plot(x,y + sqrt(2)/Cs,'green');
    plot(x - sqrt(2)/Cs,y,'green');
    plot(x - sqrt(2)/Cs,y + sqrt(2)/Cs,'green');
    plot(x,y - sqrt(2)/Cs,'green');
    plot(x - sqrt(2)/Cs,y - sqrt(2)/Cs,'green');
    plot(x + sqrt(2)/Cs,y - sqrt(2)/Cs,'green');
    theta_test = theta_0 + dtheta(index) *(0:N(index)-1);
    scatter(cos(theta_test),sin(theta_test),'m','filled');
    scatter(cos(theta_test) + sqrt(2)/Cs,sin(theta_test),'b','filled');
    scatter(cos(theta_test) + sqrt(2)/Cs,sin(theta_test) + sqrt(2)/Cs,'b','filled');
    scatter(cos(theta_test),sin(theta_test) + sqrt(2)/Cs,'b','filled');
    scatter(cos(theta_test) - sqrt(2)/Cs,sin(theta_test) + sqrt(2)/Cs,'b','filled');
    scatter(cos(theta_test),sin(theta_test) - sqrt(2)/Cs,'b','filled');
    scatter(cos(theta_test) -sqrt(2)/Cs,sin(theta_test),'b','filled');
    scatter(cos(theta_test) - sqrt(2)/Cs,sin(theta_test) - sqrt(2)/Cs,'b','filled');
    scatter(cos(theta_test) + sqrt(2)/Cs,sin(theta_test),'b','filled');
    scatter(cos(theta_test) + sqrt(2)/Cs,sin(theta_test) - sqrt(2)/Cs,'b','filled');
    %scatter(cos(thetas),sin(thetas),'black','filled');
    axis equal
end