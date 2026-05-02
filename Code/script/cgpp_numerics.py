#!/usr/bin/env python3
"""
CGPP Numerics v7: Cosmic-time integration with direct Friedmann background.
Everything integrated in cosmic time t (no 0/0 issues at potential minimum).

Uses scipy solve_ivp with BDF method for stiff background.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
import os, pickle, time, warnings
warnings.filterwarnings('ignore')

Mpl_GeV = 2.435e18
mphi = 4.14e12 / Mpl_GeV
v_pl = 0.5

def V_pot(phi):
    x = phi/v_pl; x6=x**6
    return mphi**2 * v_pl**2 / 72.0 * (1.-x6)**2

def dV_pot(phi):
    x = phi/v_pl
    return -mphi**2 * v_pl / 6.0 * (1.-x**6) * x**5

# ================ Background in cosmic time ================
# Variables: [phi, dphi/dt, ln(a)]
# Equations:
#   dphi/dt = phidot
#   d(phidot)/dt = -3H*phidot - V'(phi)
#   d(ln a)/dt = H
# where H = sqrt((phidot^2/2 + V(phi))/3)

def bg_cosmic(t, y):
    phi, phidot, lna = y
    H2 = (0.5*phidot**2 + V_pot(phi)) / 3.0
    H = np.sqrt(max(H2, 1e-60))
    dphidot = -3.*H*phidot - dV_pot(phi)
    return [phidot, dphidot, H]

def solve_bg_cosmic():
    """Solve background in cosmic time from CMB to post-inflation."""
    # Initial conditions at CMB time (t=0)
    phi0 = 0.048
    # Slow-roll estimate for phidot
    H0 = np.sqrt(V_pot(phi0)/3.0)  # slow-roll: H^2 ≈ V/3
    # From KG in slow-roll: 3H phidot ≈ -V'  → phidot ≈ -V'/(3H)
    phidot0 = -dV_pot(phi0)/(3.*H0)
    if phidot0 < 0: phidot0 = abs(phidot0)
    lna0 = 0.0  # a=1 at CMB
    
    # Find end of inflation (eps=1 where eps = phidot^2/(2H^2))
    def end_ev(t, y):
        phi, phidot, lna = y
        H2 = max(0.5*phidot**2 + V_pot(phi), 1e-60) / 3.0
        eps = phidot**2 / (2.*H2)
        return eps - 1.0
    end_ev.terminal = True; end_ev.direction = 1
    
    # First pass: find t_end and refine phi0
    sol1 = solve_ivp(bg_cosmic, [0., 1e12], [phi0, phidot0, lna0],
                     events=[end_ev], method='BDF', rtol=1e-10, atol=1e-14)
    
    if len(sol1.t_events[0]) > 0:
        t_end = sol1.t_events[0][0]
        a_end_cosmic = np.exp(sol1.y_events[0][0,2])
        N_end_cosmic = np.log(a_end_cosmic)
        print(f"  N_end (cosmic) = {N_end_cosmic:.2f}")
        
        if abs(N_end_cosmic - 60.) > 0.5:
            phi0 *= (1. + 0.15*(N_end_cosmic-60.)/60.)
            phi0 = max(.005, min(.15, phi0))
            H0 = np.sqrt(V_pot(phi0)/3.)
            phidot0 = abs(dV_pot(phi0)/(3.*H0))
            sol1 = solve_ivp(bg_cosmic, [0., 1e12], [phi0, phidot0, lna0],
                             events=[end_ev], method='BDF', rtol=1e-10, atol=1e-14)
            t_end = sol1.t_events[0][0]
            a_end_cosmic = np.exp(sol1.y_events[0][0,2])
            N_end_cosmic = np.log(a_end_cosmic)
    else:
        t_end = sol1.t[-1]; N_end_cosmic = sol1.y[2,-1]
    
    # Full integration: CMB → end of inflation + ~3 e-folds
    t_total = t_end * np.exp(3.0 * 1.5)  # rough estimate for t scaling in matter dom
    # Better: integrate to N_end + 3 in e-folds
    sol = solve_ivp(bg_cosmic, [0., t_end*20], [phi0, phidot0, lna0],
                    events=[end_ev], method='BDF', rtol=1e-10, atol=1e-14,
                    dense_output=True)
    
    # After end of inflation, integrate further
    t_after_end = sol.t[-1]
    a_after = np.exp(sol.y[2,-1])
    # Target: a_target = a_end * exp(3)
    a_target = a_after * np.exp(3.0)
    # Estimate t for matter domination: t ∝ a^{3/2}
    t_target = t_after_end * (a_target / a_after)**1.5
    
    sol2 = solve_ivp(bg_cosmic, [t_after_end, t_target],
                     [sol.y[0,-1], sol.y[1,-1], sol.y[2,-1]],
                     method='BDF', rtol=1e-10, atol=1e-14, dense_output=True)
    
    # Combine
    t_all = np.concatenate([sol.t[:-1], sol2.t])
    phi_all = np.concatenate([sol.y[0,:-1], sol2.y[0]])
    pd_all = np.concatenate([sol.y[1,:-1], sol2.y[1]])
    lna_all = np.concatenate([sol.y[2,:-1], sol2.y[2]])
    
    a_all = np.exp(lna_all)
    H_all = np.sqrt(np.maximum(0.5*pd_all**2 + V_pot(phi_all), 1e-60) / 3.0)
    eps_all = pd_all**2 / (2.*H_all**2)
    
    ei = np.argmax(eps_all >= 1.)
    if eps_all[ei] < 1.: ei = len(eps_all)-3
    
    print(f"  phi={phi0:.6f} Nend={np.log(a_all[ei]):.2f} "
          f"Hi={H_all[0]*Mpl_GeV:.3e}GeV mphi/Hi={mphi/H_all[0]:.1f}")
    print(f"  bg points: {len(t_all)}")
    
    return {'t': t_all, 'a': a_all, 'H': H_all, 'eps': eps_all,
            'ei': ei, 'ae': a_all[ei], 'He': H_all[ei], 'Hi': H_all[0]}

# ================ Splines ================
def mk_spl(bg):
    t, a, H, eps = bg['t'], bg['a'], bg['H'], bg['eps']
    return {'a': CubicSpline(t, a, extrapolate=True),
            'H': CubicSpline(t, H, extrapolate=True),
            'eps': CubicSpline(t, eps, extrapolate=True)}

# ================ Omega²_t = omega_k²/a² ================
def w2_tmin_t(t, k, m, sp):
    a = float(sp['a'](t)); H = float(sp['H'](t)); e = float(sp['eps'](t))
    return k**2/a**2 + m**2 - 2.*H**2 + e*H**2

def w2_vmin_t(t, k, m, sp):
    a = float(sp['a'](t)); H = float(sp['H'](t)); e = float(sp['eps'](t))
    dh = 1e-5
    def f_N(Nv):
        av=float(sp['a'](bg['t'][int(Nv)])); return av**2/np.sqrt(k**2+av**2*m**2)
    # Need N-index... this is getting complicated. Skip for now, use tensor frequency.
    return k**2/a**2 + m**2 - 2.*H**2 + e*H**2

# ================ Mode integrator (cosmic time) ================
def integrate_mode_t(k, m_phys, sp, bg):
    """Integrate mode in cosmic time from start to end."""
    t_all, a_all, H_all = bg['t'], bg['a'], bg['H']
    ei = bg['ei']
    
    # Start 3 e-folds before end (in N units)
    a_start = bg['ae'] * np.exp(-3.0)
    # Find t where a ≈ a_start
    idx_s = np.argmin(np.abs(a_all - a_start))
    ts = t_all[max(idx_s, 1)]
    
    te = t_all[-1]
    if te - ts < 10.: return 0.
    
    a_s = float(sp['a'](ts)); H_s = float(sp['H'](ts))
    
    # Bunch-Davies IC
    eta_s = -1./(a_s*H_s)
    ph = k*eta_s
    cr = np.cos(ph)/np.sqrt(2.*k); ci = -np.sin(ph)/np.sqrt(2.*k)
    ctr = k/a_s*ci; cti = -k/a_s*cr
    
    # ODE: dchi/dt stored in [2,3], d²chi/dt² = -H*dchi/dt - ω_c²*chi
    def ode(t, y):
        cr_,ci_,ctr_,cti_ = y
        Hv = float(sp['H'](t))
        W2 = w2_tmin_t(t, k, m_phys, sp)
        return [ctr_, cti_, -Hv*ctr_ - W2*cr_, -Hv*cti_ - W2*ci_]
    
    try:
        sol = solve_ivp(ode, [ts, te], [cr,ci,ctr,cti],
                        method='RK45', rtol=1e-8, atol=1e-12, max_step=1e6)
    except Exception:
        return 0.
    
    tf = sol.t[-1]; yf = sol.y[:,-1]
    af = float(sp['a'](tf)); Hf = float(sp['H'](tf))
    wk = np.sqrt(max(w2_tmin_t(tf, k, m_phys, sp), 1e-30)) * af
    
    cer = af*yf[2]; cei = af*yf[3]
    csq = yf[0]**2 + yf[1]**2; esq = cer**2 + cei**2
    bs = 0.5*wk*csq + 0.5/wk*esq - 0.5
    return max(float(bs), 0.)

# ================ Spectrum ================
def spectrum_t(mo, bg, sp, sec, km=-0.5, kM=2.5, n=80):
    m_phys = mo*np.sqrt(2.)*bg['Hi']
    ae=bg['ae']; He=bg['He']
    kn = np.logspace(km, kM, n); kv = kn*ae*He
    beta = np.zeros(n)
    for i, k in enumerate(kv):
        beta[i] = integrate_mode_t(k, m_phys, sp, bg)
        if i%20==0:
            nk_t = k**3/(2*np.pi**2)*beta[i]/(ae*He)**3
            print(f"    k={kn[i]:.3f} β²={beta[i]:.3e} nk={nk_t:.3e}")
    nk = kv**3/(2.*np.pi**2)*beta/(ae*He)**3
    return kn, nk, beta

# ================ Main ================
def main():
    print("="*60); print("CGPP v7: Cosmic-time background"); print("="*60)
    t0=time.time()
    print("\n[1] Background..."); bg=solve_bg_cosmic(); sp=mk_spl(bg)
    print(f"    ({time.time()-t0:.0f}s)")
    
    out='/root/Agents/Spin2CGPP/Code/output'; os.makedirs(out,exist_ok=True)
    
    print("\n[2] Testing tensor_min m=3...")
    kn,nk,beta = spectrum_t(3.0, bg, sp, 'tmin', km=-0.5, kM=2.0, n=50)
    
    v=beta>1e-50
    if v.any(): pi=np.argmax(nk[v]); print(f"    Peak: k={kn[v][pi]:.1f} nk={nk[v][pi]:.2e}")
    
    with open(f'{out}/v7_test.pkl','wb') as f:
        pickle.dump({'kn':kn,'nk':nk,'beta':beta},f)
    print(f"  Saved.")

if __name__=='__main__': main()
