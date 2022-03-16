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

save_downsample = 100

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
    f = f0
    i_conv = 0
    fs = []
    pitches = []
    for i_conv in tqdm(range(i), total=i):
        f = selfconvolve(f) * dE
        dx = dx * 2
        f = f / np.sum(f*dE)
   
    gc.collect()
   
    return f

def weighted_f(fdx, xs, wirep, smearing):
    fw = np.zeros(len(fdx))
    fw[0] = 1. / dEw
    ws = smeared_dep(xs, wirep, smearing)
    for w in tqdm(ws[ws >= 1. / (n_dist_sample / downsample)], total=len(ws[ws >= 1. / (n_dist_sample / downsample)])):
        f_x = np.zeros(len(fdx))
        if w < 1. / (n_dist_sample / downsample):
            continue
        else:
            nskip = int(1. / w)
            toset = np.unique(np.floor(np.arange((n_dist_sample // downsample))/w).astype(int)) // nskip
            toset = toset[toset < ((n_dist_sample // downsample) // nskip)]
            nset = len(toset)
            f_x[:nset] = np.sum(fdx[:len(fdx)-(len(fdx) % nskip)].reshape(-1, nskip), axis=1)[toset]
            f_x = f_x / np.sum(f_x*dEw)

            fw = doconvolve(fw, f_x)[:len(fdx)] * dEw
            fw = fw / np.sum(fw*dEw)
    return fw

def main(smearing):
    print("Making f0 distribution")
    f = make_f0()

    # fdx = f[::downsample]
    fdx = np.mean(f.reshape(-1, downsample), axis=1)
    
    
    fdx = fdx / np.sum(fdx*dEw)
    xs = np.linspace(-pitch*(Nconvolve/2), pitch*(Nconvolve/2), Nconvolve+1)
    xcenter = (xs[1:] + xs[:-1]) / 2.

    print("Making weighted distribution")
    fw = weighted_f(fdx, xcenter, wirep, smearing)

    outf = ("f0_muonE%.0f_smeared%.3fcm" % (muonE, smearing)).replace(".", "_") + ".txt"
    print("Saving to %s" % outf)
    with open(outf, "w") as f:
        for w in fw[::save_downsample]:
            f.write(str(w))

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("USAGE: python make_eloss_dist.py <Muon Energy [MeV] <Transverse Smearing [cm]>")
    else:
        import builtins 
        builtins.muonE = float(sys.argv[1])
        from constants import *
        main(float(sys.argv[2]))
