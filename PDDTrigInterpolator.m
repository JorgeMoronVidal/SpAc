function P = PDDTrigInterpolator(xi,x,y)
% TRIGINTERP Trigonometric interpolation.
% Input:
%   xi  evaluation points for the interpolant (vector)
%   x   equispaced interpolation nodes (vector, length N)
%   y   interpolation values (vector, length N)
% Output:
%   P   values of the trigonometric interpolant (vector)
N = length(x);
n = N/2;
%l = linspace(0,N-1,N);
%theta_l = l*2*pi/N;
P = zeros(length(xi),1);
for i = 1:length(xi)
    for l = 1:N
        for k = -n:n-1
        P(i) = P(i) + 1/N * y(l) * cos(k*(xi(i)-x(l)));
        end
    end
end
return 