%An interpolator invented by me to interpolate functions defined in a grid 
% in a circle 
radius = 0.7071;
M = 96;
N = 61;
M_test = 131;
N_test = 55;
N2 = (N-1)/2;
[~,r] = cheb(N,radius);
dt = 2*pi/M; t = dt*(1:M)';
[rr,tt] = meshgrid(r(1:N2+1),t([1:M]));
dt_test = 2*pi/M_test; t = dt_test * (1:M_test)';
dr_test = radius/N_test; r = dr_test * (1:N_test)';
%r = r(1:2);
[rrr,ttt] = meshgrid(r,t);
rrr_interp = rrr(:);
ttt_interp = ttt(:);
uu = Poisson3_u(rr.*cos(tt),rr.*sin(tt));
uuu = Poisson3_u(rrr.*cos(ttt),rrr.*sin(ttt));
figure
hold on
surf(rrr.*cos(ttt),rrr.*sin(ttt),uuu)
uuu_i = Interpolant_circle_u(rr,tt,uu,rrr_interp,ttt_interp);
scatter3(rrr_interp.*cos(ttt_interp),rrr_interp.*sin(ttt_interp),uuu_i);
comp = Poisson3_u(rrr_interp.*cos(ttt_interp),rrr_interp.*sin(ttt_interp))
figure
hold on
scatter3(rrr_interp.*cos(ttt_interp),rrr_interp.*sin(ttt_interp),uuu_i-comp)