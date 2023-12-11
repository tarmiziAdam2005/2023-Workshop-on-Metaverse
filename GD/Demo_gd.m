clc;
clear all;
close all;

N = 256; n = N^2;
x = double(imread('foot.png'));
%x = rgb2gray(x);
x = imresize(x,[N N]);
sigma = 0.09;  % This is the additive gaussian noise level. Th larger, the noisier the image
y = x(:) + sigma*max(x(:))*randn(n,1);




%{

figure; 
imshow(y,[]); %Show Noisy Image
title(sprintf('Noisy Image (PSNR = %3.3f dB)',psnr_fun(y,x)));
%}

mu = 0.05;
nit = 1000;
tol = 1e-10;
lamda = 1;
[v,Err] = GD(y,mu,lamda,nit,tol, n);

figure;
colormap gray;
subplot(221); imagesc(x); axis image; title(sprintf('Original')) ;

subplot(222); imagesc(reshape(y,N,N)); axis image; title(sprintf('Noisy(PSNR = %3.3f dB,SSIM = %3.3f)',psnr_fun(y,x),ssim_index(reshape(y,N,N),x)));

subplot(223); imagesc(reshape(v,N,N)); axis image; 
title(sprintf('Tikhonov(PSNR = %3.3f dB,SSIM = %3.3f)',psnr_fun(v,x),ssim_index(reshape(v,N,N),x)));

figure;
semilogy(Err,'Linewidth',2.5,'Color','black');
xlabel('Iterations (k)','FontSize',20,'Interpreter', 'latex');
ylabel('Relative Error','FontSize',20,'Interpreter', 'latex');
axis tight;
grid on;
l = legend('Gradient Descent', 'latex');
set(l,'interpreter','latex','FontSize', 25);
set(gca, 'FontSize',20)
