%% Initialization and Computation Parameters

% Clearing the workspace
clear;
clc;
% close all hidden;

%%% Domain (integers with sum up to 10)
scope = 4;
oversampling = 4;

% Defining the interval of the domain (treated as square)
x = -(2^(oversampling)-1/(2^scope)):(1/(2^scope)):(2^(oversampling)); %-1/(2^scope)

% Print intermediate Steps
verbose = false;

%% Main Parameters

%%% Wavefront definition

% Phase: OSA/ANSI Zernike Polynomials (line1: index, line2: coefficient)
wf_phs = [0     1     2     3     4     5;
          1  -0.5  -0.5     0     0     0];

% Amplitude: polar coordinates (with apodization)
fwhm = 0.4; %Full width at half maximum
wf_amp = @(r, theta) exp(-0.5*(r/(fwhm/(sqrt(2*log(2))))).^2); %gassian beam

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Difraction Grating and Spectral Response Definition

% Dispersion factor
disp = 1; % lines per mm

% Spectrum of the source
src = [0 0 0 0 0 1 0 0 0 0 0]; %from 300 to 1000 nm

% Sectral Response of System
sys = [0 0 0 0 0 1 0 0 0 0 0]; %from 300 to 1000 nm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Scaling Factors

% Plate scale
f = 1000; %focal lenght in mm
plt_scl = 206265/f; %arcseconds per mm

% Intensity of source
src_int = 1; %Watts per mm^2

%% Generating Wavefront from Polynomials and Apodization

wf = wavefront(wf_phs, wf_amp, x, verbose);
wf = wf/sum(sum(abs(wf)));


%% Taking the Fourier Transform of the Wavefront

figure(2);
F = fft2(wf);
F3 = fftshift(F);
PSF = F3.*conj(F3);

surf(PSF); shading interp;
colorbar;
view(2);
axis tight; axis square;

%% Convoluting the Transform with the Spectrum of the Source times Spectral Response of the System

%% Generating Sensor Simulation

%number of pixels in one direction
pixn = 2^(scope);
pix = (2*2^oversampling*2^scope)/(pixn);

step = (length(PSF))/pixn;
for ii = 1:pixn
    for jj = 1:pixn
        CCD(ii,jj) = sum(sum(PSF((2 + (ii-1)*step):((ii)*step),(2 + (jj-1)*step):((jj)*step))));
        jj = jj + 1;
    end
    ii = ii + 1;
end

wsx = 0; wsy = 0;
for ii = 1:length(CCD(1,:))
    wsx = ii*sum(CCD(ii,:)) + wsx;
    wsy = ii*sum(CCD(:,ii)) + wsy;
end
centx = wsx/sum(sum(CCD));
centy = wsy/sum(sum(CCD));

figure(3)
% imshow(flip(CCD./max(max(CCD))),'InitialMagnification','fit');
% imshow(flip(CCD/282),'InitialMagnification','fit');
image = bar3(flip(CCD/282), 1); view(2);
for k = 1:length(image)
    zdata = image(k).ZData;
    image(k).CData = zdata;
end
set(image, 'EdgeColor', 'none');
colormap gray; colorbar;
axis tight; axis square;
hold on;
plot(centy,pixn-centx+1,'o')