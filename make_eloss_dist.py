from tqdm.auto import tqdm
import numpy as np
import sys
import gc

downsample = 1
pitch = 1e-2
Nconvolve = 1000
wirep = 0.3

n_dist_sample = 50_000_000
# n_dist_sample = 100000000
this_maxE = 10
# this_maxE = 20
    
dE = this_maxE / n_dist_sample
dEw = dE*downsample

save_downsample = 20

def make_f0():
    energies = np.linspace(0, this_maxE, n_dist_sample, endpoint=False)
    
    dx0 = pitch
    i = 0
    while dx0 > 1e-8:
        i += 1
        dx0 = dx0 / 2
    
    f0 = np.zeros((n_dist_sample,))
    f0[energies >= I0] = dx0 * xsec(energies[energies >= I0]) * density 
    f0[0] += (1 - density * total_xsec * dx0) / dE

    dx = dx0
    f_fft = np.fft.rfft(f0)
    for i_conv in tqdm(range(i), total=i):
        f_fft = f_fft * f_fft * dE
        dx = dx * 2
        f_fft = f_fft / np.real(np.sum(f_fft*dE) / len(f_fft))
    
        gc.collect()
    
    return np.fft.irfft(f_fft)

def weighted_f(fdx, xs, wirep, smearing):
    fdx_fft = np.fft.rfft(fdx)
    fw = np.zeros(len(fdx))
    fw[0] = 1. / dEw
    fw_fft = np.fft.rfft(fw)
    
    ws = smeared_dep(xs, wirep, smearing)
    for w in tqdm(ws[ws >= 1. / (n_dist_sample / downsample)], total=len(ws[ws >= 1. / (n_dist_sample / downsample)])):
        f_x_fft = np.zeros(len(fdx_fft), dtype=np.complex_)
        if w < 1. / (n_dist_sample / downsample):
            continue
        else:
            toset = np.around(np.arange(n_dist_sample//downsample//2+1) * w).astype(int)    
            f_x_fft[:] = fdx_fft[toset]
            f_x_fft = f_x_fft / np.real(np.sum(f_x_fft*dEw) / len(f_x_fft))

            fw_fft = fw_fft * f_x_fft * dEw
            fw_fft = fw_fft / np.real(np.sum(fw_fft*dEw) / len(fw_fft))
            
    return np.fft.irfft(fw_fft)

def main(smearing):
    print("Making f0 distribution")
    f = make_f0()

    fdx = np.mean(f.reshape(-1, downsample), axis=1)
    
    fdx = fdx / np.sum(fdx*dEw)
    xs = np.linspace(-pitch*(Nconvolve/2), pitch*(Nconvolve/2), Nconvolve+1)
    xcenter = (xs[1:] + xs[:-1]) / 2.

    print("Making weighted distribution")
    fw = weighted_f(fdx, xcenter, wirep, smearing)

    outf = ("f0_muonE%.0f_smeared%.3fcm" % (muonE, smearing)).replace(".", "_") + ".txt"
    print("Saving to %s" % outf)
    with open(outf, "w") as f:
        for w in np.mean(fw.reshape(-1, save_downsample), axis=1):
            f.write(str(w) + "\n")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("USAGE: python make_eloss_dist.py <Muon Energy [MeV] <Transverse Smearing [cm]>")
    else:
        import builtins 
        builtins.muonE = float(sys.argv[1])
        from constants import *
        main(float(sys.argv[2]))
