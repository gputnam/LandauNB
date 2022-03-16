import numpy as np
from scipy.special import erf

# constants -- units MeV, cm, s
# MeV * cm
# hbarc = 1.973269804 * 1e-11
# cm / s
# c = 2.998 * 1e10

# CHANGE me
# muonE = 1000.

Mmuon = 105.6
Melec = 0.5110
Relec = 2.817940 * 1e-13

# motion -- beta
beta = np.sqrt(1 - (Mmuon * Mmuon) / (muonE * muonE))
# motion -- gamma
gamma = muonE / Mmuon
# motion -- (gamma * beta)^2
gammaB2 = gamma*gamma*beta*beta

# electron ionization energy [MeV]
I0 = 188.0 * 1e-6

# maximum energy transfer [Mev]
# maxE = 2 * Melec * gammaB2
maxE = 2 * Melec * gammaB2 / \
    ( 1 + 2 * gamma * Melec / Mmuon + (Melec/Mmuon)**2)

# get the density of argon 
# density -- g / mL
LAr_density_gmL = 1.396
# molar mass -- g / mol
Ar_molar_mass = 39.9623
# avogadro number
mole = 6.0221409*1e23
# charge number / atomic mass (Z / A)
# 99.6% of LAr is Ar40
Ar_ZA = 18. / Ar_molar_mass

# electron number density (N / cm^3)
density = Ar_ZA * mole * LAr_density_gmL

# wirep (cm)
wirep = 0.3

# MeV cm^2 / mol
K = 4*np.pi*mole*Relec**2*Melec # 0.307075

# outputdir = "/home/grayputnam/Work/Summer2020/July13BetheBlochBarkas/"
zeta = (K/2.)*Ar_ZA*(1./beta**2) * LAr_density_gmL / maxE

def f_kappa(thickness):
    return zeta * thickness

mean_dEdx = LAr_density_gmL * K * Ar_ZA * (1 / (beta * beta)) * (\
        (1./2)*np.log( maxE*2 * gammaB2*Melec / (I0**2)) \
        - beta * beta)

def f_mpv(thickness):
    kappa = f_kappa(thickness)*maxE
    j=0.200
    return (kappa/thickness)*\
    (np.log(2 * Melec * gammaB2 / I0) + np.log(kappa/I0) + j - beta*beta)

def smeared_dep(x0, w, sigma):
    if sigma < 1e-4:
        return 1*((x0 > -w/2) & (x0 <= w/2))
    return (1./2)*(erf((w/2+x0)/(np.sqrt(2)*sigma)) +\
                   erf((w/2-x0)/(np.sqrt(2)*sigma)))

def f_thickness(sigma, a=wirep):
    return a*np.exp(quad(lambda x: -smeared_dep(x, a, sigma) * np.log(smeared_dep(x, a, sigma))/a, 
                       -(a/2) - 5*sigma, (a/2) + 5*sigma)[0])

Dtransverse = 5.85e-3 # cm^2/ms
def smearing(driftT):
    return np.sqrt(2*Dtransverse*driftT)

def xsec(e):
    # front constant of xsec
    const = (2 * np.pi  * Relec * Relec * Melec) / (beta * beta)
    A = 1
    C = -beta * beta / maxE
    ret = const * (A / (e * e) + C / e)
    ret[e > maxE] = 0
    return ret

def xsec_integral(e): 
    if e > maxE:
        return xsec_integral(maxE, sign, correction)
    # front constant of xsec
    const = (2 * np.pi  * Relec * Relec * Melec) / (beta * beta)
    A = 1
    C = -beta * beta / maxE
    ret = const * (A * (1./I0 - 1./e) + \
            C * np.log(e/I0))
    return ret
total_xsec = xsec_integral(maxE)

def selfconvolve(a):
    a_fft = np.fft.rfft(a)
    a_fft_2 = a_fft * a_fft
    del a_fft
    return np.fft.irfft(a_fft_2)

def doconvolve(a, b):
    a_fft = np.fft.rfft(a)
    b_fft = np.fft.rfft(b)
    ab_fft = a_fft*b_fft
    del a_fft
    del b_fft
    return np.fft.irfft(ab_fft)
