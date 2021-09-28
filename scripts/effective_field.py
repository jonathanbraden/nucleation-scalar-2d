#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def load_log_data(fName,nSim):
    d = np.genfromtxt(fName)
    # To add: Automatically determine number of time steps
    nv = d.shape[-1]
    return np.genfromtxt(fName).reshape(nSim,-1,nv)

def add_vp_labels(a):
    a.set_xlabel(r'$\left\langle\phi\right\rangle_{\rm V}$')
    a.set_ylabel(r"$V'_{\rm eff}$")

def plot_veff(d,vp=lambda x: x,vp3=lambda x : 0, vp4=lambda x : 0, a=None, i_phi=1, i_vp=8, i_sig=9):
    if a==None:
        f,a = plt.subplots()

    sig0 = np.mean(d[:,0,i_sig])  # Use numerics to estimate sigma^2
    phi = d[:,:,i_phi].flatten()
    d_vp = d[:,:,i_vp].flatten()
    sig = d[:,:,i_sig].flatten()

    sns.kdeplot(phi, d_vp,                          shade=True, shade_lowest=False, ax=a, alpha=0.5,label=r'$\left\langle V^\prime\right\rangle$')
    sns.kdeplot(phi, d_vp-vp(phi),                  shade=True, shade_lowest=False, ax=a, alpha=0.5,label=r'$\left\langle V^\prime\right\rangle - V_{\rm T}^\prime$')
    sns.kdeplot(phi, d_vp-vp(phi)-0.5*sig0*vp3(phi), shade=True, shade_lowest=False, ax=a, alpha=0.5, label=r'$\left\langle V^\prime\right\rangle - V_{\rm T}^\prime - \frac{\sigma^2_0}{2}V^{\prime\prime\prime}$')
    #sns.kdeplot(phi, d_vp-vp(phi)-0.5*sig*vp3(phi), shade=False, shade_lowest=False, ax=a, alpha=0.5, label=r'$\left\langle V^\prime\right\rangle - V_{\rm T}^\prime - \frac{\sigma^2}{2}V^{\prime\prime\prime}$')
    # Subtract Gaussian piece
    # Add scatter plot to show KDE estimate is ok
    # Add curve for V'_tree
    
    add_vp_labels(a)
    return f,a

def plot_veff_gauss(d, vp, vp_g, vp4=lambda x : 0, a=None, i_phi=1, i_vp=8, i_sig = 9):
    if a==None:
        f,a = plt.subplots()

    i_phi = 1; i_vp = 8; i_sig = 9
    sig0 = np.mean(d[:,0,i_sig]) 
    phi = d[:,:,i_phi].flatten(); d_vp = d[:,:,i_vp].flatten(); sig = d[:,:,i_sig].flatten()
    sns.kdeplot(phi, d_vp,               shade=True, shade_lowest=False, ax=a, alpha=0.5, label=r'$\left\langle V^\prime\right\rangle$')
    sns.kdeplot(phi, d_vp-vp(phi),       shade=True, shade_lowest=False, ax=a, alpha=0.5, label=r'$\left\langle V^\prime\right\rangle - V_{\rm T}^\prime$')
    sns.kdeplot(phi, d_vp-vp_g(phi,sig0), shade=True, shade_lowest=False, ax=a, alpha=0.5, label=r'$\left\langle V^\prime\right\rangle - V_{\rm T}^\prime - \Delta V^\prime_{\rm G}(\sigma_0)$')
    sns.kdeplot(phi, d_vp-vp_g(phi,sig), shade=True, shade_lowest=False, ax=a, alpha=0.5, label=r'$\left\langle V^\prime\right\rangle - V_{\rm T}^\prime - \Delta V^\prime_{\rm G}(\sigma)$')
    # Add the first nonGaussian correction
    # Add scatter plot to show KDE has converged
    
    add_vp_labels(a)
    return f,a

def veff_plot_dw(fName='log-dw-512.out',eps=0.):
    f,a = plt.subplots()
    d = load_log_data(fName,nSim=10)
    f,a = plot_veff(d,lambda x : x*(x**2-1.)+eps*(x**2-1.), lambda x : 6.*x + 2.*eps, lambda x : 6.)
    return f,a

def veff_plot_quartic(fName='log-quartic-l1.out',lam=0.):
    f,a = plt.subplots()
    d = load_log_data(fName,nSim=1)
    f,a = plot_veff(d, lambda x : x+lam*x**3, lambda x : 6*lam*x, lambda x : 6*lam)
    return f,a

def veff_plot_sg(fName='log-sg-n512-kc64-phi1.out'):
    f,a = plt.subplots()
    d = load_log_data(fName,nSim=1)
    f,a = plot_veff(d, lambda x : 0.5*np.sin(2.*x), lambda x : -2.*np.sin(2.*x), lambda x : -4.*np.cos(2.*x))
    return f,a

if __name__=="__main__":
    f,a = veff_plot_quartic()
    #f.legend()
    pass
