function wf = wavefront(phase, amplitude, domain, verbose)

[X,Y] = meshgrid(domain,domain);
[theta,r] = cart2pol(X,Y);
idx = r<=1;
phs = nan(size(X));
amp = nan(size(X));

y = zernfun2(phase(1, 1),r(idx),theta(idx)).*phase(2, 1);
for ii = 2:length(phase(1, :))
    y = zernfun2(phase(1, ii),r(idx),theta(idx)).*phase(2, ii) + y;
end

phs(idx) = y(:,1);
g = amplitude(r, theta);
amp(idx) = g(idx);


wf = amplitude(r, theta).*exp((2*pi.*phs).*1i);

wf(isnan(wf)) = 0;

if verbose == true
    figure
    subplot(1,2,1);
    surf(domain,domain,phs); 
    shading interp;
    colormap cool;
    colorbar;
    view(2)
    axis square;
    
    subplot(1,2,2);
    surf(domain,domain,amp); 
    shading interp;
    colormap cool;
    colorbar;
    view(2)
    axis square;
end

end