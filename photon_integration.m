function integration = photon_integration(Probability_density_matrix,Number_of_photons,verbose)

integration = zeros(size(Probability_density_matrix));

Uniform_probability = unifrnd(0,sum(sum(Probability_density_matrix)), 1, Number_of_photons);

total = 0;
Cumulative_probability_matrix =  zeros(size(Probability_density_matrix));
for jj = 1: size(Probability_density_matrix, 2)
    for ii = 1:size(Probability_density_matrix, 1)
        Cumulative_probability_matrix(ii, jj) = Probability_density_matrix(ii, jj) + total;
        total = Cumulative_probability_matrix(ii, jj);
    end
end

for count = 1:Number_of_photons

last = 0;
photon = [0, 0];
for jj = 1:size(Cumulative_probability_matrix, 2)
    for ii = 1:size(Cumulative_probability_matrix, 1)
        if (Uniform_probability(count) <= Cumulative_probability_matrix(ii, jj) && Uniform_probability(count) > last)
            photon = [ii, jj];
        end
        last = Cumulative_probability_matrix(ii, jj);
    end
end

integration(photon(1), photon(2)) = integration(photon(1), photon(2)) + 1;

end

if verbose
    figure;
    image = bar3(integration, 1);
    for k = 1:length(image)
        zdata = image(k).ZData;
        image(k).CData = zdata;
    end

    colormap gray; colorbar;
    axis tight; axis square;
end