clc;
%clear all;
%close all;

N = 256; n = N^2;
x = double(imread('foot.png'));
x = imresize(x,[N N]);
sigma = 0.09;  % This is the additive gaussian noise level. Th larger, the noisier the image
y = x(:) + sigma*max(x(:))*randn(n,1);

mu = 0.005;
nit = 1000;
tol = 1e-10;
lamda = 2.2;
[z,Err3,tg1] = HeavyBall(y,mu,lamda,nit,tol, n);

figure;
colormap gray;
subplot(221); imagesc(x); axis image; title(sprintf('Original')) ;

subplot(222); imagesc(reshape(y,N,N)); axis image; title(sprintf('Noisy(PSNR = %3.3f dB,SSIM = %3.3f)',psnr_fun(y,x),ssim_index(reshape(y,N,N),x)));
subplot(223); imagesc(reshape(z,N,N)); axis image; 
title(sprintf('Tikhonov(PSNR = %3.3f dB,SSIM = %3.3f)',psnr_fun(z,x),ssim_index(reshape(z,N,N),x)));

figure;
semilogy(Err3,'Linewidth',2.5,'Color','green');
xlabel('Iterations (k)','FontSize',20,'Interpreter', 'latex');
ylabel('Relative Error','FontSize',20,'Interpreter', 'latex');
axis tight;
grid on;

l = legend('HB','Gradient Descent');
set(l,'interpreter','latex','FontSize', 25);
set(gca, 'FontSize',20)

