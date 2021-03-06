clear; clc; %close all hidden;

A1 = readmatrix('Barycenter_sat10_noise10.txt');
A2 = readmatrix('Barycenter_sat10_noise1.txt');
A3 = readmatrix('Barycenter_sat10_noise0.txt');

figure;
scatter(A1(:,1), A1(:,2), '.'); hold on;
scatter(A2(:,1), A2(:,2), '.'); hold on;
scatter(A3(:,1), A3(:,2), '.'); hold off;
axis square; %axis([28, 31, 28, 31]);

A = A3;

figure;
histogram2(A(:,1), A(:,2))
axis tight; axis square;

mean(A(:,1))
mean(A(:,2))

var(A(:,1))
var(A(:,2))

figure
fit = @(x) exp(-0.5*(x.^2)/var(A(:,2)));
x_ = -100*var(A(:,2)):0.01*var(A(:,2)):100*var(A(:,2));
plot(x_, fit(x_))
axis tight; axis square;
