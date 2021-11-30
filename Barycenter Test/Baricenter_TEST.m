clear; clc; %close all hidden;

A1 = readmatrix('Barycenter_sat1_noise0.txt');
A2 = readmatrix('Barycenter_sat10_noise0.txt');
A3 = readmatrix('Barycenter_sat100_noise0.txt');
A4 = readmatrix('Barycenter_sat1000_noise0.txt');
A5 = readmatrix('Barycenter_sat10000_noise0.txt');

figure(1);
scatter(A1(1:1000,1)-mean(A1(1:1000,1)), A1(1:1000,2)-mean(A1(1:1000,2)), '.'); hold on;
scatter(A2(1:1000,1)-mean(A2(1:1000,1)), A2(1:1000,2)-mean(A2(1:1000,2)), '.'); hold on;
scatter(A3(1:1000,1)-mean(A3(1:1000,1)), A3(1:1000,2)-mean(A3(1:1000,2)), '.'); hold on;
scatter(A4(1:1000,1)-mean(A4(1:1000,1)), A4(1:1000,2)-mean(A4(1:1000,2)), '.'); hold on;
scatter(A5(:,1)-mean(A5(:,1)), A5(:,2)-mean(A5(:,2)), '.'); hold on;
axis equal; %axis([28, 31, 28, 31]);

figure(2);
lfvar1 = fitlm(log10([1000 10000 100000 1000000 10000000]), log10([var(A1(:,1)) var(A2(:,1))  var(A3(:,1))  var(A4(:,1))  var(A5(:,1))]))
plotAdded(lfvar1)
xlabel("log10(Number of photons)");
ylabel("log10(Variance of the barycenter)");
title("");
% loglog([1000 10000 100000 1000000], [var(A1(:,1)) var(A2(:,1))  var(A3(:,1))  var(A4(:,1))], 'o')

figure(3);
lfvar2 = fitlm(log10([1000 10000 100000 1000000 10000000]), log10([var(A1(:,2)) var(A2(:,2))  var(A3(:,2))  var(A4(:,2))  var(A5(:,2))]))
plotAdded(lfvar2)
xlabel("log10(Number of photons)");
ylabel("log10(Variance of the barycenter)");
title("");
% loglog([1000 10000 100000 1000000], [var(A1(:,2)) var(A2(:,2))  var(A3(:,2))  var(A4(:,2))], 'o')

%%%%%%%%%%%%

% A = A4;
% 
% figure;
% histogram2(A(:,1), A(:,2))
% axis tight; axis square;
% 
% mean(A(:,1))
% mean(A(:,2))
% 
% var(A(:,1))
% var(A(:,2))
% 
% figure
% fit = @(x) exp(-0.5*(x.^2)/var(A(:,2)));
% x_ = -100*var(A(:,2)):0.01*var(A(:,2)):100*var(A(:,2));
% plot(x_, fit(x_))
% axis tight; axis square;
