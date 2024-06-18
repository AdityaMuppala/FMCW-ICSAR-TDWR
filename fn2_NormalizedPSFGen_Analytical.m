function [fhatPSF,x2,y2] = fn2_NormalizedPSFGen_Analytical(fhatPSF_0,xPeakPos,yPeakPos,kmax,N,M1)

%-- This code generates the general PSF function H from H_0.

os_factor = 1;
N2 = N*os_factor+1; % Should be large enough to satisfy Nyquist (i.e. N > 2*2*kmax/(2*pi))

kxy = linspace(-os_factor*2*kmax,os_factor*2*kmax,N2);

dxy = 2*pi/(2*max(kxy));
xy = dxy*(-(N2-1)/2:(N2-1)/2);

[kx2D,ky2D] = ndgrid(kxy,kxy);

fhat_BP_FT4 = fhatPSF_0.*exp(-1i*kx2D.*xPeakPos*(M1/N2)-1i*ky2D.*yPeakPos*(M1/N2));

fhat_BP_FT5 = fftInterpolate(fhat_BP_FT4,[M1 M1]);
x2 = linspace(min(xy),max(xy),M1);
y2 = linspace(min(xy),max(xy),M1);

fhatPSF = iftx(ifty(fhat_BP_FT5));


% figure(2); clf
% surf(x2D,y2D,abs(fhatPSF));
% shading interp;
% xlabel('X (m)');
% ylabel('Y (m)');
% zlabel('f(x,y) - Target Scene');
% hold on;
% axis equal; axis tight
% colormap gray;
% view([0 90])
% % xlim([-1 1])
% % ylim([-1 1])
% % set(gca, 'clim', [0 1]);
% colorbar
% w = waitforbuttonpress;

%%
end
