#!/usr/bin/env python3
"""
CGPP Paper-Match Numerics v3: N-time mode integration with full post-inflation
inflaton oscillations, using the Hamiltonian fast symplectic integrator.

Key improvements:
  - Background includes full inflaton oscillations post-inflation.
  - Background converted to uniform N-grid for fast interpolation.
  - Hamiltonian variables u = chi * sqrt(aH) → u'' + Q² u = 0 (no friction).
  - Exact-step (symplectic) midpoint integrator for speed and accuracy.
  - All 6 sectors: minimal/nonminimal × tensor/vector/scalar.
  - Integrates to a ~ 676 a_e (N_end + 6.5 e-folds).
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import os, sys, time, pickle, warnings
warnings.filterwarnings('ignore')

Mpl_GeV = 2.435e18
mphi_GeV_val = 4.14e12
v_pl = 0.5

def mphi():
    return mphi_GeV_val / Mpl_GeV  # ~1.700e-6

def V_hill(phi):
    x = phi / v_pl
    return mphi()**2 * v_pl**2 / 72.0 * (1.0 - x**6)**2

def dV_hill(phi):
    x = phi / v_pl
    return -mphi()**2 * v_pl / 6.0 * (1.0 - x**6) * x**5

# ============================================================
# Background: Inflation in N-time + post-inflation in t-time
# Converted to uniform N-grid
# ============================================================

def solve_background(phi0=0.04799, N_post=7.0):
    """
    Full background:
    1. N-time inflation (0 to N_end)
    2. Post-inflation coordinate-time integration
    3. Convert to unified N-grid with derived quantities
    Returns dictionary with interpolators on N.
    """
    # --- Phase 1: Inflation in N-time ---
    psi0 = abs(dV_hill(phi0) / V_hill(phi0))
    H0 = np.sqrt(V_hill(phi0) / (3.0 - 0.5 * psi0**2))
    
    def bg_ode_N(N, y):
        phi, psi, lnH = y
        H = np.exp(np.clip(float(lnH), -60, 10))
        eps = 0.5 * psi**2
        dpsi = -(3.0 - eps) * psi - dV_hill(phi) / H**2
        dlnH = -eps
        return [psi, dpsi, dlnH]
    
    def ev_eps1(N, y):
        return 0.5 * y[1]**2 - 1.0
    ev_eps1.terminal = True; ev_eps1.direction = 1
    
    sol_inf = solve_ivp(bg_ode_N, [0, 80], [phi0, psi0, np.log(H0)],
                        events=[ev_eps1], method='RK45',
                        rtol=1e-10, atol=1e-14, max_step=0.5)
    
    N_end = sol_inf.t_events[0][0] if len(sol_inf.t_events[0]) > 0 else sol_inf.t[-1]
    
    N_inf = sol_inf.t[sol_inf.t <= N_end]
    phi_inf = sol_inf.y[0, :len(N_inf)]
    psi_inf = sol_inf.y[1, :len(N_inf)]
    H_inf = np.exp(np.clip(sol_inf.y[2, :len(N_inf)], -60, 10))
    eps_inf = 0.5 * psi_inf**2
    a_inf = np.exp(N_inf)
    
    phi_e = phi_inf[-1]; phi_dot_e = psi_inf[-1] * H_inf[-1]
    H_e = H_inf[-1]; a_e = a_inf[-1]; t_e = 0.0
    
    # Compute t_e
    for i in range(1, len(N_inf)):
        dN = N_inf[i] - N_inf[i-1]
        H_mid = 0.5 * (H_inf[i-1] + H_inf[i])
        t_e += dN / H_mid
    
    H_inf_est = np.mean(H_inf[:max(10, len(H_inf)//20)])
    
    # --- Phase 2: Post-inflation in N-time ---
    # Evolve from N_end to N_end+N_post using N-time ODE
    # Variables: phi, psi = phi_dot/H, lnH
    N_target = N_end + N_post
    
    def bg_ode_N_post(N, y):
        phi, psi, lnH = y
        H = np.exp(np.clip(float(lnH), -60, 10))
        eps_val = 0.5 * psi**2
        dpsi = -(3.0 - eps_val) * psi - dV_hill(phi) / H**2
        dlnH = -eps_val
        return [psi, dpsi, dlnH]
    
    # Use DOP853 for stiff oscillatory phase
    # Start from N_end conditions
    sol_post = solve_ivp(bg_ode_N_post, [N_end, N_target],
                         [phi_e, psi_inf[-1], np.log(H_e)],
                         method='DOP853', rtol=1e-9, atol=1e-14,
                         max_step=0.05)
    
    N_post_arr = sol_post.t
    phi_post = sol_post.y[0]
    psi_post = sol_post.y[1]
    H_post = np.exp(np.clip(sol_post.y[2], -60, 10))
    eps_post = 0.5 * psi_post**2
    a_post = np.exp(N_post_arr)
    
    # Compute phi_dot for post-inflation
    phi_dot_post = psi_post * H_post
    
    # Compute conformal time for post-inflation
    # Conformal time at N_end: η_e = -1/(a_e H_e)
    # dη/dN = 1/(aH) → integrate
    eta_post = np.zeros(len(N_post_arr))
    eta_post[0] = eta_e
    for i in range(1, len(N_post_arr)):
        dN = N_post_arr[i] - N_post_arr[i-1]
        a_mid = 0.5 * (a_post[i-1] + a_post[i])
        H_mid = 0.5 * (H_post[i-1] + H_post[i])
        eta_post[i] = eta_post[i-1] + dN / (a_mid * H_mid)
    
    # Trim to target N range
    mask = N_post_arr <= N_target
    idx_max = np.sum(mask)
    if idx_max >= len(N_post_arr):
        idx_max = len(N_post_arr) - 1
    
    N_post_arr = N_post_arr[:idx_max+1]
    a_post = a_post[:idx_max+1]
    H_post = H_post[:idx_max+1]
    eps_post = eps_post[:idx_max+1]
    phi_post_v = phi_post[:idx_max+1]
    phi_dot_post = phi_dot_post[:idx_max+1]
    eta_post = eta_post[:idx_max+1]
    
    # --- Combine into unified N-time grid ---
    # Subsample inflationary data
    N_inf_sub = N_inf[::5]
    if N_inf_sub[-1] != N_end:
        N_inf_sub = np.append(N_inf_sub, N_end)
    a_inf_sub = np.exp(N_inf_sub)
    
    # Combine
    N_all = np.concatenate([N_inf_sub, N_post_arr[1:]])
    
    # Interpolate all quantities to N-grid
    def _interp(t_src, vals, kind='linear'):
        f = interp1d(t_src, vals, kind=kind, fill_value='extrapolate', copy=False)
        return f
    
    # For inflaton sector, use N grid directly
    f_a_inf = _interp(N_inf, a_inf)
    f_H_inf = _interp(N_inf, H_inf)
    f_eps_inf = _interp(N_inf, eps_inf)
    f_phi_inf = _interp(N_inf, phi_inf)
    
    # For post-inflation, N_post_arr maps to quantities
    f_a_post = _interp(N_post_arr, a_post)
    f_H_post = _interp(N_post_arr, H_post)
    f_eps_post = _interp(N_post_arr, eps_post)
    f_eta_post = _interp(N_post_arr, eta_post)
    f_phi_post = _interp(N_post_arr, phi_post_v)
    f_phidot_post = _interp(N_post_arr, phi_dot_post)
    
    # Build combined functions
    def get_a(N):
        return float(f_a_inf(N) if N <= N_end else f_a_post(N))
    def get_H(N):
        return float(f_H_inf(N) if N <= N_end else f_H_post(N))
    def get_eps(N):
        return float(f_eps_inf(N) if N <= N_end else f_eps_post(N))
    def get_phi(N):
        return float(f_phi_inf(N) if N <= N_end else f_phi_post(N))
    def get_eta(N):
        if N <= N_end:
            a = get_a(N); H = get_H(N)
            return -1.0 / (a * H)
        return float(f_eta_post(N))
    def get_phidot(N):
        if N <= N_end:
            H = get_H(N); eps = get_eps(N)
            return np.sqrt(2.0 * eps) * H
        return float(f_phidot_post(N))
    
    # Build arrays on uniform N-grid for fast vectorized evaluation
    N_grid = np.linspace(0, N_all[-1], max(2000, int(N_all[-1] * 60)))
    a_grid = np.array([get_a(n) for n in N_grid])
    H_grid = np.array([get_H(n) for n in N_grid])
    eps_grid = np.array([get_eps(n) for n in N_grid])
    phi_grid = np.array([get_phi(n) for n in N_grid])
    eta_grid = np.array([get_eta(n) for n in N_grid])
    
    # eps_prime = deps/dN using simple finite difference (smooth enough with dense grid)
    eps_p_grid = np.gradient(eps_grid, N_grid)
    
    # Interpolators on N
    kind = 'linear'; fill = 'extrapolate'
    bg = {
        'N_arr': N_grid, 'a_arr': a_grid, 'H_arr': H_grid,
        'eps_arr': eps_grid, 'eps_p_arr': eps_p_grid,
        'phi_arr': phi_grid, 'eta_arr': eta_grid,
        'N_end': N_end, 'a_e': a_e, 'H_e': H_e,
        'eta_e': eta_e, 'H_inf': H_inf_est,
        'N_min': N_grid[0], 'N_max': N_grid[-1],
    }
    for name in ['a', 'H', 'eps', 'eps_p', 'phi', 'eta']:
        bg[f'f_{name}'] = interp1d(N_grid, bg[f'{name}_arr'], kind=kind,
                                    fill_value=fill, copy=False)
    
    def eval_bg(N):
        return (float(bg['f_a'](N)), float(bg['f_H'](N)),
                float(bg['f_eps'](N)), float(bg['f_eps_p'](N)))
    bg['eval'] = eval_bg
    
    print(f"  N_end = {N_end:.2f}  N_range = [{N_grid[0]:.2f}, {N_grid[-1]:.2f}]")
    print(f"  H_inf = {H_inf_est*Mpl_GeV:.3e} GeV  H_e = {H_e*Mpl_GeV:.3e} GeV")
    print(f"  m_phi/H_inf = {mphi()/H_inf_est:.1f}  H_inf/H_e = {H_inf_est/H_e:.3f}")
    print(f"  a_max/a_e = {a_grid[-1]/a_e:.1f}")
    print(f"  N_grid: {len(N_grid)} points, ΔN ≈ {N_grid[1]-N_grid[0]:.4f}")
    
    return bg

# ============================================================
# Omega² functions (comoving)
# ============================================================

def o2_tensor_min(k, a, H, eps, m, **kw):
    return k**2 + a**2 * (m**2 - H**2 * (2.0 - eps))

def o2_tensor_nonmin(k, a, H, eps, m, **kw):
    return k**2 + a**2 * (m**2 + H**2 * (1.0 - eps))

def o2_vector_min(k, a, H, eps, m, **kw):
    a2, k2, m2 = a**2, k**2, m**2
    a2m2 = a2 * m2; D = k2 + a2m2
    num = (6.0*k2**2 + 5.0*k2*a2m2 + 2.0*a2m2**2
           - eps*(2.0*k2**2 + 3.0*k2*a2m2 + a2m2**2))
    return k2 + a2m2 - a2*H**2*num/D**2

def o2_vector_nonmin(k, a, H, eps, eps_p, m, **kw):
    """Vector nonminimal via finite difference of K_C in N-time."""
    a2, k2, m2 = a**2, k**2, m**2; a4 = a**4
    if H == 0: return k2
    
    Hdot_t = -H**2 * eps  # dH/dt
    mu1sq = m2 + 3.0*H**2 - Hdot_t/H
    mu2sq = m2 + 3.0*H**2 + 2.0*Hdot_t/H
    if mu1sq <= 0: return k2 + a2*m2
    D0 = k2 + a2*mu1sq
    if D0 <= 0: return k2
    
    KC0 = a4*k2*mu1sq/D0; MC0 = a4*k2*mu2sq
    
    dN = 0.0005
    a_p = a*np.exp(dN); H_p = H*np.exp(-eps*dN)
    eps_p_use = eps + eps_p*dN
    Hdot_p = -H_p**2 * eps_p_use
    mu1sq_p = m2 + 3.0*H_p**2 - Hdot_p/H_p
    D_p = k2 + a_p**2*mu1sq_p
    KC_p = a_p**4*k2*mu1sq_p/max(D_p,1e-60) if mu1sq_p>0 and D_p>0 else KC0
    
    a_m = a*np.exp(-dN); H_m = H*np.exp(eps*dN)
    eps_m_use = eps - eps_p*dN
    Hdot_m = -H_m**2 * eps_m_use
    mu1sq_m = m2 + 3.0*H_m**2 - Hdot_m/H_m
    D_m = k2 + a_m**2*mu1sq_m
    KC_m = a_m**4*k2*mu1sq_m/max(D_m,1e-60) if mu1sq_m>0 and D_m>0 else KC0
    
    KC_p_N = (KC_p - KC_m)/(2.0*dN)
    KC_pp_N = (KC_p - 2.0*KC0 + KC_m)/dN**2
    
    num = 4.0*KC0*MC0 + KC_p_N**2 - 2.0*KC0*KC_pp_N
    den = 4.0*KC0**2
    return num/den if den > 0 else k2

def _KB_MB(k, a, H, eps, eps_p, m):
    """Exact K_B, M_B for nonminimal scalar."""
    k2,k4,k6,k8,k10 = k**2,k**4,k**6,k**8,k**10
    a2,a4,a6 = a**2,a**4,a**6
    H2,H4,H6 = H**2,H**4,H**6
    m2v,m4,m6 = m**2,m**4,m**6
    
    dotH = -H2*eps  # dH/dt = H * dH/dN = H * (-H*eps) = -H² eps
    # d²H/dt²: d(dotH)/dN = -2H*dotH*dN/dN - H²*d(eps)/dN 
    # Actually dotH = H * H_dot_N = H * (-H*eps) = -H² eps
    # d(dotH)/N = -2H*H_N*eps - H²*eps_N = -2H*(-H*eps)*eps - H²*eps_p = 2H²*eps² - H²*eps_p
    # ddotH_t = d(dotH)/dt = H * d(dotH)/dN = H*(2H²*eps²-H²*eps_p) = H³(2eps²-eps_p)
    ddotH = H**3 * (2.0*eps**2 - eps_p)
    
    m2pH2 = m2v + H2
    m2p3H2 = m2v + 3.0*H2
    m2p3H2m_dH = m2p3H2 - dotH
    m2p3H2p2dH = m2p3H2 + 2.0*dotH
    
    P = (4.0*(m2p3H2+3.0*dotH)*k4
         + 12.0*a2*(m2pH2*m2p3H2 + 2.0*m2pH2*dotH - dotH**2)*k2
         + 9.0*a4*m2pH2*m2p3H2m_dH*m2p3H2p2dH)
    if P <= 1e-60: return 1e-100, 1e-100
    
    KB = a4/P*(-4.0*dotH**2*k6 + 3.0*a2*m2pH2*m2p3H2m_dH*m2p3H2p2dH*k4)
    
    c10 = (12.0*m2pH2*m2p3H2**3
           + 16.0*m2p3H2**2*(6.0*m2v+7.0*H2)*dotH
           + 4.0*m2p3H2*(63.0*m2v+71.0*H2)*dotH**2
           + 8.0*(25.0*m2v+27.0*H2)*dotH**3 - 48.0*dotH**4
           - 32.0*H*m2p3H2*dotH*ddotH - 48.0*H*dotH**2*ddotH)
    
    brkt = (2.0*m2pH2*m2p3H2**2*(2.0*m2v+5.0*H2)
            + m2p3H2*(19.0*m4+64.0*m2v*H2+49.0*H4)*dotH
            + 2.0*(7.0*m4+20.0*m2v*H2+17.0*H4)*dotH**2
            - (23.0*m2v+25.0*H2)*dotH**3 + 2.0*dotH**4
            - 2.0*H*m2pH2*m2p3H2*ddotH - 4.0*H*m2pH2*dotH*ddotH)
    c8 = 12.0*a2*m2p3H2p2dH*brkt
    c6 = 9.0*a4*m2pH2*m2p3H2m_dH*m2p3H2p2dH**2*(7.0*m2pH2*m2p3H2+17.0*m2pH2*dotH-8.0*dotH**2)
    c4 = 27.0*a6*m2pH2**2*m2p3H2m_dH**2*m2p3H2p2dH**3
    
    MB = a6/P**2*(c10*k10 + c8*k8 + c6*k6 + c4*k4)
    return KB, MB

def o2_scalar_nonmin(k, a, H, eps, eps_p, m, **kw):
    """Nonminimal scalar omega² via finite-diff of log K_B."""
    KB0, MB0 = _KB_MB(k, a, H, eps, eps_p, m)
    if KB0 <= 1e-100: return k**2
    
    dN = 0.001
    a_f = a*np.exp(dN); H_f = H*np.exp(-eps*dN)
    eps_f = eps + eps_p*dN; eps_p_f = eps_p
    KB_f, _ = _KB_MB(k, a_f, H_f, eps_f, eps_p_f, m)
    
    a_b = a*np.exp(-dN); H_b = H*np.exp(eps*dN)
    eps_b = eps - eps_p*dN; eps_p_b = eps_p
    KB_b, _ = _KB_MB(k, a_b, H_b, eps_b, eps_p_b, m)
    
    L0 = np.log(max(KB0,1e-100))
    Lp = np.log(max(KB_f,1e-100)); Lm = np.log(max(KB_b,1e-100))
    L_dN = (Lp-Lm)/(2.0*dN); L_d2N = (Lp-2.0*L0+Lm)/dN**2
    
    a2H2 = a**2*H**2
    O2_N = MB0/(a2H2*max(KB0,1e-100)) - 0.5*L_d2N - 0.5*(1.0-eps)*L_dN - 0.25*L_dN**2
    return max(O2_N*a2H2, 1e-60)

O2_FUNCS = {
    'tensor_min': o2_tensor_min,
    'tensor_nonmin': o2_tensor_nonmin,
    'vector_min': o2_vector_min,
    'vector_nonmin': o2_vector_nonmin,
    'scalar_nonmin': o2_scalar_nonmin,
}

# ============================================================
# Fast symplectic integrator in N-time
# ============================================================

def integrate_mode(k, m, bg, sector, N_start, N_max,
                   dN_step=0.005, A_threshold=0.01, N_buffer=0.7):
    """
    Fast integration using exact-step midpoint method:
    u_{n+1} = u_n cos(Q_mid dN) + u'_n sin(Q_mid dN)/Q_mid
    u'_{n+1} = -u_n Q_mid sin(Q_mid dN) + u'_n cos(Q_mid dN)
    
    where u = chi * sqrt(aH) satisfies u'' + Q² u = 0.
    Q² = Omega² + eps'/2 - (1-eps)²/4.
    """
    o2_func = O2_FUNCS[sector]
    N_end = bg['N_end']
    
    # Bunch-Davies IC at N_start
    a_s, H_s, eps_s, eps_p_s = bg['eval'](N_start)
    aH_s = a_s * H_s
    eta_s = float(bg['f_eta'](N_start))
    
    chi0 = np.exp(-1j * k * eta_s) / np.sqrt(2.0 * k)
    chi_N0 = -1j * k * chi0 / aH_s
    F0 = np.sqrt(aH_s)
    ur = chi0.real * F0; ui = chi0.imag * F0
    uNr = chi_N0.real * F0 + chi0.real * F0 * 0.5 * (1.0 - eps_s)
    uNi = chi_N0.imag * F0 + chi0.imag * F0 * 0.5 * (1.0 - eps_s)
    
    N = N_start
    while N < N_max - 1e-10:
        a, H, eps, eps_p = bg['eval'](N)
        
        # Adaptive step based on Q
        O2 = o2_func(k, a, H, eps, eps_p=eps_p, m=m)
        Q_cur = np.sqrt(max(O2/(a**2*H**2) + 0.5*eps_p - 0.25*(1.0-eps)**2, 1e-60))
        dN = min(0.05, 0.5 / max(Q_cur, 0.1))
        
        if N + dN > N_max:
            dN = N_max - N
        
        # Midpoint evaluation
        N_mid = N + 0.5 * dN
        a_mid, H_mid, eps_mid, eps_p_mid = bg['eval'](N_mid)
        O2_mid = o2_func(k, a_mid, H_mid, eps_mid, eps_p=eps_p_mid, m=m)
        Q2_mid = max(O2_mid/(a_mid**2*H_mid**2) + 0.5*eps_p_mid
                     - 0.25*(1.0-eps_mid)**2, 1e-60)
        Q_mid = np.sqrt(Q2_mid)
        
        cQ, sQ = np.cos(Q_mid*dN), np.sin(Q_mid*dN)
        ur_n = ur*cQ + uNr*sQ/Q_mid
        ui_n = ui*cQ + uNi*sQ/Q_mid
        uNr_n = -ur*Q_mid*sQ + uNr*cQ
        uNi_n = -ui*Q_mid*sQ + uNi*cQ
        
        ur, ui, uNr, uNi = ur_n, ui_n, uNr_n, uNi_n
        N += dN
    
    # Extract chi and beta² at final time
    a_f, H_f, eps_f, eps_p_f = bg['eval'](N_max)
    aH_f = a_f * H_f; F_f = np.sqrt(aH_f)
    fp_f = 0.5 * (1.0 - eps_f)
    
    chi_r = ur / F_f; chi_i = ui / F_f
    chi_N_r = uNr / F_f - ur * fp_f / F_f
    chi_N_i = uNi / F_f - ui * fp_f / F_f
    
    O2_f = o2_func(k, a_f, H_f, eps_f, eps_p=eps_p_f, m=m)
    Om_f = np.sqrt(max(O2_f/(aH_f**2), 1e-60))
    
    beta_sq = 0.5*aH_f*Om_f*(chi_r**2+chi_i**2) + 0.5*aH_f/Om_f*(chi_N_r**2+chi_N_i**2) - 0.5
    wron = 2.0*aH_f*(chi_r*chi_N_i - chi_i*chi_N_r)
    
    if abs(wron) > 1e-6:
        beta_sq_f = max((beta_sq + 0.5)/abs(wron) - 0.5, 0.0)
    else:
        beta_sq_f = max(beta_sq, 0.0)
    
    # Adiabaticity check
    a_p = a_f*np.exp(0.01); H_p = H_f*np.exp(-eps_f*0.01)
    eps_p_val = eps_f + eps_p_f*0.01
    O2_p = o2_func(k, a_p, H_p, eps_p_val, eps_p=eps_p_f, m=m)
    a_m = a_f*np.exp(-0.01); H_m = H_f*np.exp(eps_f*0.01)
    eps_m_val = eps_f - eps_p_f*0.01
    O2_m = o2_func(k, a_m, H_m, eps_m_val, eps_p=eps_p_f, m=m)
    Q2_p = max(O2_p/(a_p**2*H_p**2) + 0.5*eps_p_f - 0.25*(1.0-eps_f)**2, 1e-60)
    Q2_m = max(O2_m/(a_m**2*H_m**2) + 0.5*eps_p_f - 0.25*(1.0-eps_f)**2, 1e-60)
    dQ_dN = abs(np.sqrt(Q2_p) - np.sqrt(Q2_m))/0.02
    Q_f = np.sqrt(max(O2_f/(aH_f**2) + 0.5*eps_p_f - 0.25*(1.0-eps_f)**2, 1e-60))
    A_val = dQ_dN / max(Q_f**2, 1e-60)
    
    return {
        'beta_sq': beta_sq_f, 'wronskian': wron,
        'N_stop': N_max, 'A_final': A_val,
    }


def compute_spectrum(bg, m_over_sqrt2_Hinf, sector='tensor_min',
                     k_min=0.03, k_max=200.0, n_k=60,
                     N_start_offset=3.0, N_extra=7.0,
                     verbose=True):
    """Compute |β|² spectrum."""
    H_inf = bg['H_inf']; a_e = bg['a_e']; H_e = bg['H_e']
    m = m_over_sqrt2_Hinf * np.sqrt(2.0) * H_inf
    N_end = bg['N_end']
    N_start = N_end - N_start_offset
    N_max = N_end + N_extra
    
    if N_max > bg['N_max']:
        N_max = bg['N_max']
    
    k_norm_arr = np.logspace(np.log10(k_min), np.log10(k_max), n_k)
    k_arr = k_norm_arr * a_e * H_e
    
    beta_sq_arr = np.zeros(n_k); wron_arr = np.zeros(n_k)
    
    if verbose:
        print(f"\n  {sector}  m/(√2 H_inf)={m_over_sqrt2_Hinf:.1f}  "
              f"k=[{k_min:.2e},{k_max:.2e}]  n={n_k}")
    
    t0 = time.time()
    for i, (kv, kn) in enumerate(zip(k_arr, k_norm_arr)):
        try:
            res = integrate_mode(kv, m, bg, sector, N_start, N_max)
            beta_sq_arr[i] = res['beta_sq']
            wron_arr[i] = res['wronskian']
            
            if verbose and (i % 10 == 0 or i == n_k - 1):
                nk = kn**3/(2.0*np.pi**2)*res['beta_sq']
                print(f"    [{i:3d}/{n_k}] k={kn:.3e}  |β|²={res['beta_sq']:.3e}  "
                      f"n_k={nk:.3e}  Wr={res['wronskian']:.4f}  "
                      f"t={time.time()-t0:.1f}s")
        except Exception as e:
            print(f"    [{i:3d}/{n_k}] ERROR: {e}")
            beta_sq_arr[i] = 0.0; wron_arr[i] = 0.0
    
    elapsed = time.time() - t0
    if verbose:
        print(f"    Done in {elapsed:.1f}s ({elapsed/n_k:.1f}s/mode)")
    
    nk_arr = k_norm_arr**3/(2.0*np.pi**2)*beta_sq_arr
    return {
        'k_norm': k_norm_arr, 'beta_sq': beta_sq_arr, 'nk': nk_arr,
        'm': m, 'm_over_sqrt2_Hinf': m_over_sqrt2_Hinf, 'sector': sector,
        'wronskian': wron_arr,
    }


def print_summary(spectra_dict):
    print("\n" + "="*70); print("SUMMARY"); print("="*70)
    for key, spec in spectra_dict.items():
        beta = spec['beta_sq']; kn = spec['k_norm']; nk = spec['nk']
        valid = (beta > 1e-60) & np.isfinite(beta)
        if not valid.any():
            print(f"\n  {key}: no valid data"); continue
        idx = np.argmax(nk[valid])
        n_tot = np.trapezoid(nk[valid], np.log(kn[valid]))
        print(f"\n  {key}:")
        print(f"    Peak k={kn[valid][idx]:.4e}  n_k={nk[valid][idx]:.4e}  "
              f"|β|²={beta[valid][idx]:.4e}")
        print(f"    ∫ n_k d ln k = {n_tot:.3e}")


def save_output(spectra_dict, bg, outdir, fname='spectra_paper_v3.pkl'):
    os.makedirs(outdir, exist_ok=True)
    fpath = os.path.join(outdir, fname)
    save_spec = {}
    for key, spec in spectra_dict.items():
        save_spec[key] = {k: spec[k] for k in
                          ['k_norm','beta_sq','nk','m','m_over_sqrt2_Hinf','sector']}
    save_bg = {k: bg[k] for k in ['N_end','a_e','H_e','H_inf','eta_e']}
    with open(fpath, 'wb') as f:
        pickle.dump({'spectra':save_spec,'background':save_bg}, f)
    print(f"\nSaved to {fpath}")


# ============================================================
# Main
# ============================================================
if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description='CGPP Paper-Match Numerics v3')
    ap.add_argument('--quick','-q',action='store_true')
    ap.add_argument('--full','-f',action='store_true')
    ap.add_argument('--sector','-s',type=str,default=None)
    ap.add_argument('--k-min',type=float,default=0.03)
    ap.add_argument('--k-max',type=float,default=200.0)
    ap.add_argument('--n-k',type=int,default=60)
    ap.add_argument('--output','-o',type=str,default='spectra_paper_v3.pkl')
    args = ap.parse_args()
    
    print("="*70); print("CGPP Paper-Match v3 (N-time symplectic)"); print("="*70)
    print("\n[1] Background...")
    bg = solve_background(N_post=7.0)
    a_e, H_e = bg['a_e'], bg['H_e']
    
    if args.quick:
        print("\n[2] Quick test...")
        for m_val in [2.0, 3.0]:
            m = m_val * np.sqrt(2.0) * bg['H_inf']
            k_norms = [0.1, 1.0, 10.0, 50.0]
            N_end = bg['N_end']; N_start = N_end - 3.0; N_max = bg['N_max']
            
            for sector in ['tensor_min', 'vector_min', 'tensor_nonmin']:
                print(f"\n  [{sector}] m/(√2 H_inf)={m_val}:")
                for kn in k_norms:
                    k = kn * a_e * H_e
                    res = integrate_mode(k, m, bg, sector, N_start, N_max)
                    nk = kn**3/(2*np.pi**2)*res['beta_sq']
                    print(f"    k/(a_eH_e)={kn:7.2f}  |β|²={res['beta_sq']:.4e}  "
                          f"n_k={nk:.4e}  Wr={res['wronskian']:.4f}  A={res['A_final']:.4f}")
        print("\nDone.")
    
    elif args.full:
        print("\n[2] Full spectra...")
        masses = [2.0, 3.0, 4.0, 5.0]
        sectors = (args.sector.split(',') if args.sector
                   else ['tensor_min','vector_min','tensor_nonmin','vector_nonmin'])
        N_end = bg['N_end']; N_start = N_end - 3.0; N_max = N_end + 6.5
        
        spectra = {}
        for sector in sectors:
            for mv in masses:
                key = f"{sector}_m{mv:.1f}"
                try:
                    spectra[key] = compute_spectrum(
                        bg, mv, sector=sector,
                        k_min=args.k_min, k_max=args.k_max, n_k=args.n_k,
                        N_start_offset=3.0, N_extra=6.5, verbose=True)
                except Exception as e:
                    print(f"    ERROR {key}: {e}")
                    import traceback; traceback.print_exc()
        
        outdir = '/root/Agents/Spin2CGPP/Code/output'
        save_output(spectra, bg, outdir, args.output)
        print_summary(spectra)
        print("\nAll done.")
    
    else:
        print(f"\n  Use --quick or --full.")
