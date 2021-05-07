%p29.m (Modified by Jorge Moron) - Poisson eq. on unit circle with homogeneous or inhomogeneous BC's
% M is the discretization in the angle t (Has to be even)
M = 120;
%N is the discretization in the radious r (Has to be odd)
N = 61;
% Laplacian in polar coordinates
[D,r] = cheb(N); N2 = (N-1)/2; D2 = D^2;
D1 = D2(1:N2+1,1:N2+1); D2 = D2(1:N2+1,N+1:-1:N2+2);
E1 = D(1:N2+1,1:N2+1); E2 = D(1:N2+1,N+1:-1:N2+2);
dt = 2*pi/M; t = dt*(1:M)'; M2 = M/2;
D2t = toeplitz([-pi^2/(3*dt^2)-1/6 ...
.5*(-1).^(2:M)./sin(dt*(1:M-1)/2).^2]);
R = diag(1./r(1:N2+1)); Z = zeros(M2); I = eye(M2);
L = kron(D1+R*E1,eye(M))+kron(D2+R*E2,[Z I;I Z])+kron(R^2,D2t);

% Coordinates are set as vectors
[rr,tt] = meshgrid(r(1:N2+1),t([1:M]));
[xx,yy] = pol2cart(tt,rr);
xx = xx(:); yy = yy(:);
rr = rr(:); tt = tt(:);

%b contains the index of the knots which support the stencil 
b = find( rr == max(rr));
L(b,:) = zeros(size(L(b,:))); L(b,b) = eye(size(L(b,b)));

%Right-hand side of Lu - c(x)*u = f(x) 
%rhs = 4*ones(size(tt));
%Poisson2_f(xx,yy);
rhs = Poisson3_f(xx,yy);


%We impose BC's
%rhs(b) = rr(b).^2;
%rhs(b) = sin(kx*pi*xx(b)).*sin(ky*pi*yy(b))+C;
%rhs(b) = Poisson3_u(xx(b),yy(b));
rhs(b) = zeros(size(b));
%for i = 1:length(b)
       %rhs(b(i)) = GBC(xx(b(i)),yy(b(i)));
%end
u = L\rhs; 
uu = reshape(u,M,N2 + 1);
[rr,tt] = meshgrid(r(1:N2+1),t([1:M]));
[xx,yy] = pol2cart(tt,rr);
[xxx,yyy] = meshgrid(-1:.04:1,-1:.04:1);
%uuu = interp2(xx,yy,uu,xxx,yyy,'spline');
clf, subplot(1,3,1)
hold on
mesh(xx,yy,uu), colormap([0 0 0])
%scatter3(p_matlab(1,:),p_matlab(2,:),u_matlab)
%scatter3(xx,yy,u);
%axis([-1 1 -1 1 0.9 3]), view(-20,45),
title('Pseudopspectrarl solution')
%u_a = rr.^2;
%u_a = sin(kx*pi*xx).*sin(ky*pi*yy)+C;
u_a = Poisson3_u(xx,yy);
subplot(1,3,2),
mesh(xx,yy,u_a), colormap([0 0 0])
%scatter3(xx,yy,u_a);
%axis([-1 1 -1 1 0.9 3]),
title('Analytical solution')
error = abs(u_a - uu);
subplot(1,3,3),
mesh(xx,yy,error), colormap([0 0 0])
title('Absolute error')
