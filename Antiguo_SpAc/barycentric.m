% Example
f = @(x) exp(x); % function
m = 30; % number of interpolation points
ts = -1+1/m:2/m:1-1/m; xs = cos(pi*(ts+1)/2);
ys = f(xs);
x = linspace(-1,1,5000);
% Barycentric interpolation
c = (-1).^(0:m-1).*sin(pi*(ts+1)/2);
numer = zeros(size(x));
denom = zeros(size(x));
exact = zeros(size(x));
for j = 1:m
xdiff = x-xs(j);
temp = c(j)./xdiff;
numer = numer + temp*ys(j);
denom = denom + temp;
exact(xdiff==0) = j;
end
y = numer./denom; jj = find(exact); y(jj) = ys(exact(jj));
% Check error
max(abs(y - f(x)))