function ui = Interpolant_circle_u(r,theta,u,r_i,theta_i)
    N = length(r(:,1));
    %n = N/2;
    %A_k = NaN(size(u));
    %k_vec = NaN(N,1);
    %for column = 1:length(r(1,:))
        %for row_A = 1:N
            %A_k(row_A,column) = 0;
            %k = -n + row_A -1;
            %k_vec(row_A) = k;
            %for row_u = 1:N
                %A_k(row_A,column) = A_k(row_A,column) + cos(-k*theta(row_u,column))*u(row_u,column);
            %end
            %A_k(row_A,column) = A_k(row_A,column) * (1.0/N);
        %end
    %end
    Ni = length(theta_i);
    ui = zeros(Ni,1);
    r_vec = [r(1,:), -flip(r(1,:))]';
    %r_vec = [r(1,:)]';
    len_r_vec = length(r_vec);
    for i = 1:Ni
        %disp(i)
        aux = NaN(len_r_vec,1);
        %We compute the value of the function at every r for a given
        %theta_i
        for col_A = 1:len_r_vec/2
        %for col_A = 1:len_r_vec
            aux(col_A) = 0;
            aux(col_A) = real(PDDTrigInterpolator(theta_i(i),theta(:,col_A),u(:,col_A)));
        end
        for col_A = 1:len_r_vec/2
            %aux(col_A + len_r_vec/2) = 0;
            %for row_A = 1:N
                %k = -n + row_A -1;
                aux(col_A + len_r_vec/2) = real(PDDTrigInterpolator(theta_i(i)+pi,theta(:,len_r_vec/2 -col_A +1),u(:,len_r_vec/2 -col_A +1)));
        end
        %end
        %figure
        %hold on 
        %scatter(r_vec,aux);
        %ref = [Poisson3_u(r_vec(1:len_r_vec/2)*cos(theta_i(i)),r_vec(1:len_r_vec/2)*sin(theta_i(i)))',...
               %Poisson3_u(-1*r_vec(len_r_vec/2 +1:len_r_vec)*cos(theta_i(i)+pi),...
               %-1*r_vec(len_r_vec/2 +1:len_r_vec)*sin(theta_i(i)+pi))']';
        %ref = [Poisson3_u(r_vec*cos(theta_i(i)),r_vec*sin(theta_i(i)))]';
        %plot(r_vec,ref-aux);
        ui(i) = interp1(r_vec,aux,r_i(i),'spline');
        %ui(i) = lagrange(r_i(i),r_vec,aux);
    end
return