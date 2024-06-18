function fhatPSF_0 = fn0_H0_Analytical(Rg,Zc,kmax,kmin,gamma,N)
%-- This code computes H_0 from Eq (21)

c  = physconst('Lightspeed');
os_factor = 1;
N2 = N*os_factor+1; % Should be large enough to satisfy Nyquist (i.e. N > 2*2*kmax/(2*pi))

kxy = linspace(-os_factor*2*kmax,os_factor*2*kmax,N2);

dxy = 2*pi/(2*max(kxy));
xy = dxy*(-(N2-1)/2:(N2-1)/2);

[x2D,y2D] = ndgrid(xy,xy);
[kx2D,ky2D] = ndgrid(kxy,kxy);

Ro = sqrt(Rg.^2+Zc.^2);
r_0 = sqrt((x2D).^2+(y2D).^2+Zc.^2);

ft_term = ftx(fty(exp(-2*1i*kmax*Ro).*exp(2*1i*kmax*r_0)./(r_0-Ro) - exp(-2*1i*kmin*Ro).*(exp(2*1i*kmin*r_0)./(r_0-Ro))));

fhatPSF_0 = (c/(16*pi*1i*gamma)).*(1/(Ro^2)).*ft_term.*besselj(0,sqrt(kx2D.^2+ky2D.^2).*Rg);

%%
end
