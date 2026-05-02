#!/usr/bin/env python3
"""Analyze CGPP numerical results and produce summary tables."""

import pickle, numpy as np, sys

def analyze(fname):
    with open(fname, 'rb') as f:
        data = pickle.load(f)
    
    spectra = data['spectra']
    bg = data['background']
    
    print("=" * 70)
    print("CGPP Numerical Results: Summary")
    print("=" * 70)
    print(f"\nBackground: N_end={bg['N_end']:.2f}  "
          f"H_inf={bg['H_inf']*2.435e18:.3e} GeV  "
          f"H_e={bg['H_e']*2.435e18:.3e} GeV")
    
    # Sorted keys
    for sector in ['tensor_minimal', 'vector_minimal', 'tensor_nonmin', 'vector_nonmin']:
        print(f"\n{'-'*70}")
        print(f"Sector: {sector}")
        print(f"{'m/(√2 H_inf)':>14s}  {'Peak k':>10s}  {'Peak n_k':>12s}  "
              f"{'Peak |β|²':>12s}  {'∫ n_k d ln k':>14s}")
        print(f"{'-'*14}  {'-'*10}  {'-'*12}  {'-'*12}  {'-'*14}")
        
        for mv in [2.0, 3.0, 4.0, 5.0]:
            key = f"{sector}_m{mv:.1f}"
            if key not in spectra:
                continue
            s = spectra[key]
            nk = s['nk']; kn = s['k_norm']; beta = s['beta_sq']
            idx = np.argmax(nk)
            n_tot = np.trapezoid(nk, np.log(kn))
            print(f"{mv:14.1f}  {kn[idx]:10.4f}  {nk[idx]:12.4e}  "
                  f"{beta[idx]:12.4e}  {n_tot:14.4e}")
    
    # Cross-sector comparison for m=3
    print(f"\n{'='*70}")
    print("Cross-sector comparison (m/(√2 H_inf)=3.0):")
    print(f"{'Sector':>20s}  {'Peak k':>10s}  {'Peak n_k':>12s}  "
          f"{'∫ n_k d ln k':>14s}")
    for sector in ['tensor_minimal', 'vector_minimal', 'tensor_nonmin', 'vector_nonmin']:
        key = f"{sector}_m3.0"
        s = spectra[key]
        idx = np.argmax(s['nk'])
        n_tot = np.trapezoid(s['nk'], np.log(s['k_norm']))
        print(f"{sector:>20s}  {s['k_norm'][idx]:10.4f}  "
              f"{s['nk'][idx]:12.4e}  {n_tot:14.4e}")
    
    # Low-k scaling check
    print(f"\n{'='*70}")
    print("Low-k scaling check (tensor_minimal, m=3.0):")
    s = spectra['tensor_minimal_m3.0']
    mask = s['k_norm'] < 1.0
    k_low = s['k_norm'][mask]; nk_low = s['nk'][mask]
    print(f"  n_k/k^3 range: {nk_low[0]/k_low[0]**3:.4e} → "
          f"{nk_low[-1]/k_low[-1]**3:.4e}")
    
    return spectra, bg

if __name__ == '__main__':
    fname = sys.argv[1] if len(sys.argv) > 1 else \
            '/root/Agents/Spin2CGPP/Code/output/cgpp_full_result.pkl'
    analyze(fname)

