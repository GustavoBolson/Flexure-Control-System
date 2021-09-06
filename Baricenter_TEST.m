clear; clc; close all hidden;

A = readmatrix('Baricenter_TEST1.csv');

figure;
scatter(A(:,1), A(:,2), '.');
axis tight; axis square;

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
