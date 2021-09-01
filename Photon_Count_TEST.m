clear; clc; close all hidden;

verbose = true;

N = 1000000;

A = [0.1, 0.0, 0.2;
     0.3, 0.05, 0.1;
     0.05, 0.15, 0.05];
CCD = zeros(size(A));

r = unifrnd(0,sum(sum(A)), 1, N);

total = 0;
B =  zeros(size(A));
for jj = 1: size(A, 2)
    for ii = 1:size(A, 1)
        B(ii, jj) = A(ii, jj) + total;
        total = B(ii, jj);
    end
end

for count = 1:N

last = 0;
photon = [0, 0];
for jj = 1:size(B, 2)
    for ii = 1:size(B, 1)
        if (r(count) <= B(ii, jj) && r(count) > last)
            photon = [ii, jj];
        end
        last = B(ii, jj);
    end
end

CCD(photon(1), photon(2)) = CCD(photon(1), photon(2)) + 1;

end

if verbose
    figure
    image = bar3(CCD/(sum(sum(CCD))), 1); view(2);
    for k = 1:length(image)
        zdata = image(k).ZData;
        image(k).CData = zdata;
    end
    % set(image, 'EdgeColor', 'none');
    colormap gray; colorbar;
    axis tight; axis square;
end