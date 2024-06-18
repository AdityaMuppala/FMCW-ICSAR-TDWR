%-- Code written by Aditya Varma Muppala for the paper titled: FMCW Inverse Circular Synthetic Aperture Radar Using a Fast Time-Domain Reconstruction
%-- Last edited on 06/18/2024. 
%-- The data and NUFFT codes are in the zipped file.
%-- Extract the zip file and place files in the same path as the MATLAB codes.
 

clc
clear
clear path

path = 'NUFFT_code';
addpath(path)

Rg = 0.92; % Array radius
Zc = 0.23; % Distance between array plane and object plane

%-- Radar System Operational Parameters
fBW = 8e9;                 % bandwidth
wBW = 2*pi*fBW;

fc = 80e9;                 % carrier frequency
wc = 2*pi*fc;

c  = 3e8;                  % RF propagation speed
lambda_min = c/(fc+fBW/2); % Wavelength at highest frequency
lambda_max = c/(fc-fBW/2); % Wavelength at lowest frequency
kmin = 2*pi/lambda_max;    % Wavenumber at lowest frequency
kmax = 2*pi/lambda_min;    % Wavenumber at highest frequency

%-- Fast-Time domain parameters and arrays
phip_deg = 52-(0:0.1:359.9); phip = deg2rad(phip_deg);

fs = 50*1e6;               %Sampling rate;
T = 48*1e-6;               % Chirp duration in microseconds
N = fs * T;                % number of fast-time samples
gamma = -wBW/T;            % Chirp rate (System uses a down-chirp so gamma is negative)
M = length(phip);          % Number of slow time measurements

t = linspace(-T/2,T/2,N);  % Fast time
R0 = 0.5;                  % Beam radius on ground

krho_min = 0;
krho_max = 2*kmax;
k_rho = linspace(krho_min,krho_max,N);

%-- Data acquisition
fileID = fopen('Measurement_1_pistol.ch0','rb');
ch1 = fread(fileID,'double');
fclose(fileID);

fileID2 = fopen('Bkg_Measurement_1_pistol.ch0','rb');
ch1_bkg = fread(fileID2,'double');
fclose(fileID2);

valid_samples_per_trigger = fs*T;
ch1b = squeeze(reshape(ch1,valid_samples_per_trigger,[],M));
sdata_time = ch1b - ch1_bkg;

%-- Phase center calibration. See Section III.B in reference [13]
del_lt = 21e-2;
exp_shift = transpose(exp(-1i*(wc*2*del_lt/c-2*gamma*del_lt*t/c)));
s_b = sdata_time.*exp_shift; % Calibrated beat signals ready for TDWR


%% Step 0: Pre-computing G_0 and Gamma_0

aux1 = ones(1,N);
Kt2D = (2/c)*(gamma*t.'+wc)*aux1;
k_rho2D = aux1.'*k_rho;

aux2 = ones(1,M);
phip2D = aux1.'*phip;
k_rho2D_b = k_rho.'*aux2;
Gamma_0 = k_rho2D_b.*exp(-1i*k_rho2D_b.*Rg.*cos(phip2D));
Gamma_0_FT = fty(Gamma_0);

kx = k_rho2D_b.*cos(phip2D);
ky = k_rho2D_b.*sin(phip2D);

downSampleFactor = 1;
dx = 2*pi/(2*max(k_rho));
M1 = 2*R0/dx; M1 = round(M1/32)*32;
M1 = M1/downSampleFactor; M1 = M1+1;

knots=[kx(:),ky(:)]; Nx=M1; Ny=M1; Desired_accuracy = 6;

dx = 2*pi/(2*max(k_rho));
dy = dx;
x = dx*(-M1/2:M1/2-1);
y = dy*(-M1/2:M1/2-1);

[x2D,y2D] = ndgrid(x,y);

G_0 = exp(-1i*(Zc*sqrt(Kt2D.^2-k_rho2D.^2)));
G_0(Kt2D<=k_rho2D) = 0;

G_0_inv = G_0'; % Could also use T-SVD

%% TDWR

tic
Gamma_FT = fty(G_0_inv*s_b);
Fp = ifty(Gamma_FT.*conj(Gamma_0_FT));

% NUFFT
f1=FGG_2d_type1(Fp(:),knots,[Nx,Ny],Desired_accuracy);
fhat_original = flip(flip(f1),2); % Restoring flips introduced between Eq. (10) and (11)

fhat = fhat_original;
time1 = toc

%% CLEAN

downPad = 2; UpPad = 100; % Upsampling and downsampling to accurately find the peaks
fhatPSF_0 = fn0_H0_Analytical(Rg,Zc,kmax,kmin,gamma,N); % Precomputing H_0

nComp = 25; % Number of CLEAN iterations to run (User defined)
Clean_Component = zeros(3,nComp);

tic
for ind = 1:nComp
    ind
    [fhatPeak,xPeakPos,yPeakPos] = fn1_peakSearch(fhat,x,y,downPad,UpPad); % Peak search
    Clean_Component(:,ind) = [fhatPeak,xPeakPos,yPeakPos]; % Populating CLEAN components

    [fhatPSF,x,y] = fn2_NormalizedPSFGen_Analytical(fhatPSF_0,xPeakPos,yPeakPos,kmax,N,M1); % Find PSF at peak location
    [fhatPeak2,~,~] = fn1_peakSearch(fhatPSF,x,y,downPad,UpPad); % Find amplitude of this PSF
    
    fhatPSFNorm = fhatPSF*(fhatPeak/fhatPeak2); % Normalize PSF to complex peak magnitude
    
    fhat = fhat-fhatPSFNorm; % Subtract PSF

end

fhat_residual = fhat; % Residual image after CLEAN

fhat_Clean = zeros(size(fhat));

%-- Restoring clean components into image

for ind = 1:nComp
    xp = Clean_Component(2,ind); yp = Clean_Component(3,ind); fp = Clean_Component(1,ind);
    [~,xpind] = min(abs(x-xp)); [~,ypind] = min(abs(y-yp));
    fhat_Clean(xpind,ypind) = fp;
end

time2 = toc

%% Plotting

fhat_final = fhat_Clean+fhat_residual;
fIm1 = abs(fhat_original)./max(abs(fhat_original(:)));
fIm2 = abs(fhat_final)./max(abs(fhat_final(:)));

% Plotting
figure(1); clf
surf(x2D*100,y2D*100,fIm1);
shading interp;
xlabel('X (cm)');
ylabel('Y (cm)');
% title('f(x,y) - Target Scene');
hold on;
axis equal; axis tight
colormap gray;
view([0 90])
xlim([-20 20])
ylim([-20 20])
set(gca, 'clim', [0.2 0.5]);
set(gcf,'position',[550,50,800,800])
set(gca,'FontSize',20)
set(gca,'GridAlpha',1)
set(gca,'GridLineStyle','--')
set(gca,'fontname','ariel')
set(gcf,'color','w');

figure(2); clf
surf(x2D*100,y2D*100,fIm2);
shading interp;
xlabel('X (cm)');
ylabel('Y (cm)');
hold on;
axis equal; axis tight
colormap gray;
view([0 90])
xlim([-20 20])
ylim([-20 20])
set(gca, 'clim', [0.2 0.5]);
set(gcf,'position',[550,50,800,800])
set(gca,'FontSize',20)
set(gca,'GridAlpha',1)
set(gca,'GridLineStyle','--')
set(gca,'fontname','ariel')
set(gcf,'color','w');

clear path
