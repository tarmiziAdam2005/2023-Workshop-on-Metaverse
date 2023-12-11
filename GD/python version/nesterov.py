# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 18:43:22 2023

@author: Ong
"""

import numpy as np
import scipy.sparse as sp
from skimage import io, transform
from skimage.metrics import structural_similarity as ssim
import matplotlib.pyplot as plt

def psnr_fun(image1, image2):
    mse = np.mean((image1 - image2) ** 2)
    if mse == 0:
        return float('inf')
    max_pixel = np.max(image1)
    psnr = 20 * np.log10(max_pixel / np.sqrt(mse))
    return psnr

def ssim_fun(image1, image2):
    return ssim(image1, image2)

def DiffOper(N):
    data = [-np.ones(N), np.ones(N)]
    offsets = [0, 1]
    D = sp.diags(data, offsets, shape=(N, N), format="csc")
    D[0, 0] = 0
    return D

def tikhonov_agd(g, mu, lamda, nit, tol, n):
    x =g.copy()
    t=1
    y=x.copy()
    D = DiffOper(int(np.sqrt(n)))
    RelErr = np.zeros(nit)

    for k in range(nit):
        x_old = x.copy()
        x = y - mu * (y-g + 2 * lamda * D.T * D * y)
        t_old=t
        
        t=(k-1)/(k+2)
        y=x+t*(x-x_old)
        
        RelErr[k] = np.linalg.norm(x - x_old, "fro") / np.linalg.norm(x, "fro")

        if RelErr[k] < tol:
            RelErr = RelErr[: k + 1]
            break

    return x, RelErr

#main 
N = 256
n = N * N
x = io.imread('foot.png', as_gray=True)
x = transform.resize(x, (N, N), anti_aliasing=True)
sigma = 0.09
y = x + sigma * np.max(x) * np.random.randn(N, N)

mu = 0.0001
nit = 10000
tol = 1e-10
lamda = 1
v, Err = tikhonov_agd(y, mu, lamda, nit, tol, n) 

# Display the results and create figures
plt.figure(figsize=(12, 8))

# Original Image
plt.subplot(231)
plt.imshow(x.reshape(N, N), cmap='gray')
plt.axis('image')
plt.title('Original')

# Noisy Image
plt.subplot(232)
plt.imshow(y.reshape(N, N), cmap='gray')
plt.axis('image')
psnr = psnr_fun(y, x)
ssim_value = ssim_fun(y, x)
plt.title(f'Noisy (PSNR = {psnr:.3f} dB, SSIM = {ssim_value:.3f})')

plt.subplot(233)
plt.imshow(v.reshape(N, N), cmap='gray')
plt.axis('image')
psnr = psnr_fun(v, x)
ssim_value = ssim_fun(v, x)
plt.title(f'Tikhonov (PSNR = {psnr:.3f} dB, SSIM = {ssim_value:.3f})')

# Plot Relative Error
plt.subplot(234)
plt.semilogy(Err, linewidth=2.5, color='black')
plt.xlabel('Iterations (k)', fontsize=12)
plt.ylabel('Relative Error', fontsize=12)
plt.axis('tight')
plt.grid()
plt.legend(['Nesterov Gradient Descent'], loc='best', fontsize=12)

plt.tight_layout()
plt.show()
