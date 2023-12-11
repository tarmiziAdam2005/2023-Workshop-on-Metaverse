

clc; clearvars;
close all;
% Created on 10/12/2023 by Tarmizi bin Adam
% This file demonstrate the time-series smoothing (trend filtering)
% using Tsengs's Alternating Minimization Algorithm (AMA)

%Load Malaysian Covid-19 daily death time-series
% Data from 17/2/2020 - 2/12/2023 (1356 Days)

load testCov19Case.mat;

y = testCov19Case';

timeSeries_y = y;  

lam = 70;
beta = 100;
tol = 1e-5;
Nit = 500;

out = amtv(timeSeries_y,lam,beta,tol, Nit); % The AMA solver for TV denoising



%% Some result plotting %%
%{
figure;
subplot(2,1,1)
plot(timeSeries_y,'Linewidth',2.5,'Color','black');

xlabel('Day','FontSize',25,'interpreter','latex');
ylabel('Deaths','FontSize',25,'interpreter','latex');
axis tight;
grid on;
set(gca, 'FontSize',20)

subplot(2,1,2)
plot(out.sol,'Linewidth',2.5,'Color','red');

xlabel('Day','FontSize',25,'interpreter','latex');
ylabel('Deaths','FontSize',25,'interpreter','latex');
axis tight;
grid on;
set(gca, 'FontSize',20)
%}

figure;
plot(timeSeries_y,'Linewidth',2,'Color','black');
hold on;
plot(out.sol,'Linewidth',2.5,'Color','red');

xlabel('Day','FontSize',25,'interpreter','latex');
ylabel('Deaths','FontSize',25,'interpreter','latex');
axis tight;
grid on;
set(gca, 'FontSize',20)


