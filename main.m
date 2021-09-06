%% Initialization and Computation Parameters

% Clearing the workspace
clear;
clc;
% close all hidden;

%%% Domain (integers with sum up to 10)
scope = 6;
oversampling = 10-scope;

% Defining the interval of the domain (treated as square)
x = -(2^(oversampling)-1/(2^scope)):(1/(2^scope)):(2^(oversampling)); %-1/(2^scope)

% Print intermediate Steps
verbose = true;
verbose_level_down = false;

%% Main Parameters

%%% Wavefront definition

% Phase: OSA/ANSI Zernike Polynomials
wf_phs = [0     1     2     3     4     7; %index
          1     0     0     0     0     0]; %coefficient

% Amplitude: polar coordinates (with apodization)
% fwhm = 1000; %0.3; %Full width at half maximum
% wf_amp = @(r, theta) exp(-0.5*(r/(fwhm/(sqrt(2*log(2))))).^2); %gassian beam

I_e2 = 0.5; %0.3; %Full width at 1/e^2 of maximum
wf_amp = @(r, theta) exp(-2*(r/I_e2).^2); %gassian beam

% wf_amp = @(r, theta) 1; %uniform beam

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Difraction Grating and Spectral Response Definition

% Dispersion factor
disp = 1; % lines per mm

% Spectrum of the source
src_fwhm = 1*10^(-9); %full width at half maximum spread of light source in [m]
src_spt = @(f) exp(-0.5*(f/(src_fwhm/(sqrt(2*log(2))))).^2); %from 300 to 1000 [nm]
src_wvl = 1000*10^(-9); % central wavelenght in [m]
src_frc = 299792458/(src_wvl);%299792458/(300*10^(-9)); %frequency of source in [Hz]

% Sectral Response of System
sys_spt = [0 0 0 0 0 1 0 0 0 0 0]; %from 300 to 1000 [nm]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Scaling Factors

% Plate scale
fcl_len = 1.08; %594*10^(-3); %focal lenght in [m]
ext_ppl = 0.18; %262*10^(-3); %exit pupil diameter in [m]

% Domain physical size
extent = (2^scope)*src_wvl*fcl_len/ext_ppl;
sample_size = extent/(2^(oversampling+scope));
positions = -extent+sample_size:sample_size:extent;

% Exposure
exp_tme = 0.000000001; %Seconds

% Power output of source
src_pow = 0.000001; %Watts
num_pht = (src_pow*exp_tme)/(6.62607015*10^(-34)*src_frc); %Numer of photons to estimate

% CCD background noise
pht_prm = 10000000; %Spurious photons per second at operation temperature

%% Generating Wavefront from Polynomials and Apodization

wf = wavefront(wf_phs, wf_amp, x, verbose_level_down);
wf = wf/sum(sum(abs(wf)));

%% Taking the Fourier Transform of the Wavefront

F = fft2(wf);
F3 = flip(fftshift(F));
PSF = F3.*conj(F3);
PSF = [PSF(:,2:length(PSF)), PSF(:,1)];

if verbose
    figure(2);
    surf(positions, positions, PSF, PSF); shading interp;
    colorbar;
    view(2);
    axis tight; axis square;
    xlabel("distance from center [m]", 'Interpreter', 'latex'); ylabel("distance from center [m]", 'Interpreter', 'latex');
    title("Point Spread Function mid-Spectrum", 'Interpreter', 'latex');
end

%% Convoluting the Transform with the Spectrum of the Source times Spectral Response of the System
% H = fspecial('motion',256,0);
H = [2*ones(1,512) 10*ones(1,512) ones(1,512)]; H = H/sum(H);
% H = [1];
Dispersed_PSF = imfilter(PSF,H,'conv');

if verbose
    figure(3);
    surf(positions, positions, Dispersed_PSF, Dispersed_PSF);shading interp;
    colorbar;
    view(2);
    axis tight; axis square;
    xlabel("distance from center [m]", 'Interpreter', 'latex'); ylabel("distance from center [m]", 'Interpreter', 'latex'); 
    title("Projection after Diffraction Grating", 'Interpreter', 'latex');
end

%% Generating Sensor Simulation

%number of pixels in one direction
pixn = 2^(scope);
pix = (2*2^oversampling*2^scope)/(pixn);

step = (length(Dispersed_PSF))/pixn;
for ii = 1:pixn
    for jj = 1:pixn
        CCD(ii,jj) = sum(sum(Dispersed_PSF((1 + (ii-1)*step):((ii)*step),(1 + (jj-1)*step):((jj)*step))));
        jj = jj + 1;
    end
    ii = ii + 1;
end

% for exposure = 1:10000
    
%simulating single photon hits
CCD = photon_integration(CCD,round(num_pht),verbose_level_down);
CCD = CCD + photon_integration(ones(pixn),round(pht_prm*exp_tme),verbose_level_down);

% imshow(flip(integration))
% axis tight; axis square;

wsx = 0; wsy = 0;
for ii = 1:length(CCD(1,:))
    wsx = ii*sum(CCD(ii,:)) + wsx;
    wsy = ii*sum(CCD(:,ii)) + wsy;
end
centx = wsx/sum(sum(CCD));
centy = wsy/sum(sum(CCD));

if verbose
    figure(5)
    % imshow(flip(CCD./max(max(CCD))),'InitialMagnification','fit');
    % imshow(flip(CCD/282),'InitialMagnification','fit');
    % CCD = CCD/sum(sum(PSF));
    image = bar3(flip(CCD), 1); view(2);
    for k = 1:length(image)
        zdata = image(k).ZData;
        image(k).CData = zdata;
    end
    set(image, 'EdgeColor', 'none');
    colormap gray; colorbar;
    axis tight; axis square;
    hold on;
    plot3(centy,pixn-centx+1,max(max(CCD)),'o'); hold off
    title("CCD Photon Count", 'Interpreter', 'latex');
end

[centy,pixn-centx+1]

% writematrix([centy,pixn-centx+1],'Baricenter_TEST2.txt','WriteMode','append');
% end