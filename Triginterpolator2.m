function P = Triginterpolator2(xi,x,y)
% TRIGINTERP Trigonometric interpolation.
% Input:
%   xi  evaluation points for the interpolant (vector)
%   x   equispaced interpolation nodes (vector, length N)
%   y   interpolation values (vector, length N)
% Output:
%   P   values of the trigonometric interpolant (vector)
N = length(x);
%n = (N-1)/2;
n = N/2;
%l = linspace(0,2*n,2*n+1);
l = linspace(0,2*n-1,2*n);
%tl = l*2*pi/(2*n+1);
tl = l*2*pi/(2*n);
clear l;
c = zeros(N,1);
for k_i = 1:N
    k = -n + k_i -1;
    %k = -n + k_i;
    for l = 1:N
        c(k_i) = c(k_i) + exp(-1i*k*x(l))*y(l);
    end
    %c(k_i) = c(k_i) * (1/(2*n+1));
    c(k_i) = c(k_i) * (1/(2*n));
end
P = zeros(length(xi),1);
xfft= fftshift(fft(y))'/(N);
for i = 1:length(xi)
    for l = 1:N
        k = -n + l -1;
        %ll = k_vec(l);
        P(i) = P(i) + c(l)*exp(1i*xi(i)*k);
    end
end

%figure 
%hold on 
%plot(real(c))
%plot(real(xfft))
return 