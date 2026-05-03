#!/usr/bin/env python3
"""
CGPP Full Numerics v2: Fast mode integration using custom oscillatory integrator.
Uses the exact-step method (Q treated as constant over each step), which is
symplectic and conserves the Wronskian perfectly for slowly-varying Q.

Units: M_pl = 1.  Time: e-fold time N = ln(a).
"""

import numpy as np
import os, sys, time, pickle

np.set_printoptions(precision=6, linewidth=160)

# ============================================================
# Physical parameters
# ============================================================
Mpl_GeV = 2.435e18
mphi_GeV = 4.14e12
mphi = mphi_GeV / Mpl_GeV          # ~1.700e-6
v_pl = 0.5

def V_hill(phi):
    x = phi / v_pl
    return mphi**2 * v_pl**2 / 72.0 * (1.0 - x**6)**2

def dV_hill(phi):
    x = phi / v_pl
    return -mphi**2 * v_pl / 6.0 * (1.0 - x**6) * x**5

# ============================================================
# Background (same as before, but faster)
# ============================================================

def solve_background(N_extra=6.0, phi0=0.04799):
    psi0 = abs(dV_hill(phi0) / V_hill(phi0))
    H0 = np.sqrt(V_hill(phi0) / (3.0 - 0.5 * psi0**2))
    lnH0 = np.log(H0)

    from scipy.integrate import solve_ivp

    def bg_ode(N, y):
        phi, psi, lnH = y
        H = np.exp(np.clip(float(lnH), -60, 10))
        eps = 0.5 * psi**2
        dpsi = -(3.0 - eps) * psi - dV_hill(phi) / H**2
        dlnH = -eps
        return [psi, dpsi, dlnH]

    def ev_eps1(N, y):
        return 0.5 * y[1]**2 - 1.0
    ev_eps1.terminal = True
    ev_eps1.direction = 1

    sol = solve_ivp(bg_ode, [0, 80], [phi0, psi0, lnH0],
                    events=[ev_eps1], method='RK45',
                    rtol=1e-10, atol=1e-14, max_step=0.5)

    Ne = sol.t_events[0][0] if len(sol.t_events[0]) > 0 else sol.t[-1]

    N_inf = sol.t
    psi_inf = sol.y[1]; lnH_inf = sol.y[2]
    H_inf_arr = np.exp(np.clip(lnH_inf, -60, 10))
    eps_inf_arr = 0.5 * psi_inf**2
    a_inf_arr = np.exp(N_inf); phi_inf = sol.y[0]

    idx_end = np.argmax(eps_inf_arr >= 1.0)
    if eps_inf_arr[idx_end] < 1.0:
        idx_end = len(eps_inf_arr) - 5
    N_end = N_inf[idx_end]; a_e = a_inf_arr[idx_end]; H_e = H_inf_arr[idx_end]
    eta_inf_dS = -1.0 / (a_inf_arr * H_inf_arr)

    # Post-inflation analytic continuation
    eps_final, Delta = 1.5, 0.5
    N_post = np.linspace(N_end, N_end + N_extra, max(60, int(N_extra * 40)))
    dN_local = N_post[1] - N_post[0]
    a_post = np.exp(N_post)
    eps_post = 1.0 + (eps_final - 1.0)*(1.0 - np.exp(-(N_post - N_end)/Delta))
    H_post = np.zeros(len(N_post)); H_post[0] = H_e
    eta_post = np.zeros(len(N_post)); eta_post[0] = -1.0/(a_e*H_e)
    for i in range(1, len(N_post)):
        H_post[i] = H_post[i-1]*np.exp(-eps_post[i-1]*dN_local)
        eta_post[i] = eta_post[i-1] + dN_local/(a_post[i-1]*H_post[i-1])

    N_post, a_post, H_post, eps_post, eta_post = \
        N_post[1:], a_post[1:], H_post[1:], eps_post[1:], eta_post[1:]

    N_arr = np.concatenate([N_inf, N_post])
    a_arr = np.concatenate([a_inf_arr, a_post])
    H_arr = np.concatenate([H_inf_arr, H_post])
    eps_arr = np.concatenate([eps_inf_arr, eps_post])
    phi_arr = np.concatenate([phi_inf, np.full_like(N_post, phi_inf[idx_end])])
    eta_arr = np.concatenate([eta_inf_dS, eta_post])
    eps_prime = np.gradient(eps_arr, N_arr)

    H_inf_est = np.mean(H_inf_arr[:max(10, len(H_inf_arr)//20)])

    return {
        'N': N_arr, 'a': a_arr, 'H': H_arr,
        'eps': eps_arr, 'eps_prime': eps_prime,
        'phi': phi_arr, 'eta': eta_arr,
        'N_end': N_end, 'a_e': a_e, 'H_e': H_e,
        'eta_e': -1.0/(a_e*H_e), 'H_inf': H_inf_est,
    }


# ============================================================
# Fast background evaluator (interpolation class)
# ============================================================

class BgInterp:
    """Vectorized background interpolation - much faster than per-call."""
    def __init__(self, bg):
        # Create fine uniform grid for fast evaluation
        self.N_bg = bg['N']
        from scipy.interpolate import interp1d
        kind = 'linear'
        self._a = interp1d(bg['N'], bg['a'], kind=kind, copy=False,
                           bounds_error=False, fill_value=(bg['a'][0], bg['a'][-1]))
        self._H = interp1d(bg['N'], bg['H'], kind=kind, copy=False,
                           bounds_error=False, fill_value=(bg['H'][0], bg['H'][-1]))
        self._eps = interp1d(bg['N'], bg['eps'], kind=kind, copy=False,
                             bounds_error=False, fill_value=(bg['eps'][0], bg['eps'][-1]))
        self._eps_p = interp1d(bg['N'], bg['eps_prime'], kind=kind, copy=False,
                               bounds_error=False, fill_value=(bg['eps_prime'][0], bg['eps_prime'][-1]))
        self._eta = interp1d(bg['N'], bg['eta'], kind=kind, copy=False,
                             bounds_error=False, fill_value=(bg['eta'][0], bg['eta'][-1]))
        self.N_end = bg['N_end']

    def __call__(self, N):
        return (float(self._a(N)), float(self._H(N)),
                float(self._eps(N)), float(self._eps_p(N)))

    def eta_at(self, N):
        return float(self._eta(N))


# ============================================================
# Fast interpolated Q² function
# ============================================================

def make_Q2_interpolator(k, m, bg, sector, N_vals):
    """Precompute Q² = Omega² + eps'/2 - (1-eps)²/4 on a grid for fast lookups."""
    O2_func = SECTORS[sector]
    a_arr = bg['a']; H_arr = bg['H']; eps_arr = bg['eps']; eps_p_arr = bg['eps_prime']
    N_arr = bg['N']
    Q2_vals = np.empty(len(N_vals))
    for i, Nv in enumerate(N_vals):
        idx = np.searchsorted(N_arr, Nv)
        ilo = max(0, min(len(N_arr)-1, idx-1))
        ihi = max(0, min(len(N_arr)-1, idx))
        if ilo == ihi:
            a, H, eps, eps_p = a_arr[ilo], H_arr[ilo], eps_arr[ilo], eps_p_arr[ilo]
        else:
            t = (Nv - N_arr[ilo])/max(N_arr[ihi]-N_arr[ilo], 1e-30)
            a = a_arr[ilo] + t*(a_arr[ihi]-a_arr[ilo])
            H = H_arr[ilo] + t*(H_arr[ihi]-H_arr[ilo])
            eps = eps_arr[ilo] + t*(eps_arr[ihi]-eps_arr[ilo])
            eps_p = eps_p_arr[ilo] + t*(eps_p_arr[ihi]-eps_p_arr[ilo])
        O2 = O2_func(k, a, H, eps, eps_p, m)
        O2_d = max(O2/(a**2*H**2), 1e-60)
        Q2_vals[i] = O2_d + 0.5*eps_p - 0.25*(1.0-eps)**2
    
    from scipy.interpolate import interp1d
    return interp1d(N_vals, Q2_vals, kind='linear', copy=False,
                    bounds_error=False, fill_value=(Q2_vals[0], Q2_vals[-1]))


# ============================================================
# omega_k^2 functions
# ============================================================

def O2_tensor_minimal(k, a, H, eps, eps_prime, m, **kw):
    return k**2 + a**2*(m**2 - H**2*(2.0 - eps))

def O2_tensor_nonmin(k, a, H, eps, eps_prime, m, **kw):
    return k**2 + a**2*(m**2 + H**2*(1.0 - eps))

def O2_vector_minimal(k, a, H, eps, eps_prime, m, **kw):
    a2, k2, m2 = a**2, k**2, m**2
    D = k2 + a2*m2
    a2H2, a2m2 = a2*H**2, a2*m2
    num = (6.0*k2**2 + 5.0*k2*a2m2 + 2.0*a2m2**2
           - eps*(2.0*k2**2 + 3.0*k2*a2m2 + a2m2**2))
    return k2 + a2m2 - a2H2*num/D**2

def O2_vector_nonmin(k, a, H, eps, eps_prime, m, **kw):
    """Non-minimal vector omega_k^2. Paper Eq (3.66).
    mu1^2 = m^2 + 3H^2 + H*eps (since Hdot = -H^2*eps).
    mu2^2 = m^2 + 3H^2 - 2H*eps.
    K_C = a^4 k^2 mu1^2 / (k^2 + a^2 mu1^2), M_C = a^4 k^2 mu2^2.
    
    omega_k^2 = (4 K_C M_C + K_C'^2 - 2 K_C K_C'')/(4 K_C^2)
    where primes are d/dN (e-fold time), not d/dη!
    
    Actually, paper uses conformal time η for mode equation: χ'' + ω² χ = 0.
    The canonical form is χ = C / sqrt(2K_C), giving ω² = (4KC*MC + KC'² - 2KC*KC'')/(4KC²)
    where KC' = dKC/dη.
    
    We compute dKC/dη = aH * dKC/dN.
    """
    a2, k2, m2 = a**2, k**2, m**2; a4 = a**4; H2 = H**2
    Hdot = -H2*eps  # dH/dt
    
    # mu1^2 and mu2^2 in physical units
    # From paper: mu1² = m² + 3H² - \dot H = m² + 3H² + H² eps = m² + H²(3+eps)
    # But wait, \dot H = dH/dt. In paper: H' = a Hdot. With Hdot = -H² eps.
    # Paper eqs after (3.65): mu1² = m² - Λ + 3H² - a^{-1}H' 
    # H' = a * Hdot = -a H² eps. So a^{-1}H' = -H² eps.
    # mu1² = m² - Λ + 3H² - (-H² eps) = m² - Λ + 3H² + H² eps = m² - Λ + H²(3+eps)
    # mu2² = m² - Λ + 3H² + 2a^{-1}H' = m² - Λ + 3H² - 2H² eps = m² - Λ + H²(3-2eps)
    # With Λ=0: mu1² = m² + H²(3+eps), mu2² = m² + H²(3-2eps)
    
    mu1sq = m2 + H2*(3.0 + eps)
    mu2sq = m2 + H2*(3.0 - 2.0*eps)
    
    if mu1sq <= 0:
        return k2 + a2*m2
    
    D0 = k2 + a2*mu1sq
    if D0 <= 0:
        return k2
    
    KC0 = a4*k2*mu1sq/D0
    MC0 = a4*k2*mu2sq
    
    # N-time finite difference for K_C' = dK_C/dη = aH * dK_C/dN
    dN = 0.0005
    a_p = a*np.exp(dN); H_p = H*np.exp(-eps*dN)
    eps_p_val = eps + eps_prime*dN
    mu1sq_p = m2 + H_p**2*(3.0 + eps_p_val)
    D_p = k2 + a_p**2*mu1sq_p
    KC_p = a_p**4*k2*mu1sq_p/max(D_p,1e-100) if mu1sq_p>0 and D_p>0 else KC0
    
    a_m = a*np.exp(-dN); H_m = H*np.exp(eps*dN)
    eps_m_val = eps - eps_prime*dN
    mu1sq_m = m2 + H_m**2*(3.0 + eps_m_val)
    D_m = k2 + a_m**2*mu1sq_m
    KC_m = a_m**4*k2*mu1sq_m/max(D_m,1e-100) if mu1sq_m>0 and D_m>0 else KC0
    
    aH = a*H
    KCprime = aH*(KC_p - KC_m)/(2.0*dN)
    KCpp = aH**2*((KC_p - 2.0*KC0 + KC_m)/dN**2 + (1.0-eps)*(KC_p-KC_m)/(2.0*dN))
    
    numer = 4.0*KC0*MC0 + KCprime**2 - 2.0*KC0*KCpp
    denom = 4.0*KC0**2
    if denom <= 1e-100:
        return k2
    return numer/denom


def O2_scalar_nonmin(k, a, H, eps, eps_prime, m, **kw):
    """Non-minimal scalar omega_k^2 via finite difference."""
    dN_fd = 0.005
    KB0, MB0 = _KB_MB(k, a, H, eps, eps_prime, m)
    if KB0 <= 1e-60:
        return k**2
    KBp, _ = _KB_MB(k, a*np.exp(dN_fd), H*np.exp(-eps*dN_fd), eps+eps_prime*dN_fd, eps_prime, m)
    KBm, _ = _KB_MB(k, a*np.exp(-dN_fd), H*np.exp(eps*dN_fd), eps-eps_prime*dN_fd, eps_prime, m)
    
    L0 = np.log(max(KB0,1e-60)); Lp = np.log(max(KBp,1e-60)); Lm = np.log(max(KBm,1e-60))
    L_dN = (Lp - Lm)/(2.0*dN_fd); L_d2N = (Lp - 2.0*L0 + Lm)/dN_fd**2
    
    a2H2 = a**2*H**2
    Omega2 = MB0/(a2H2*max(KB0,1e-60)) - 0.5*L_d2N - 0.5*(1.0-eps)*L_dN - 0.25*L_dN**2
    return max(Omega2*a2H2, 1e-60*a2H2)


def _KB_MB(k, a, H, eps, eps_prime, m):
    """K_B and M_B for non-minimal scalar. Paper Eqs (3.69)-(3.70).
    dotH = dH/dt, ddotH = d²H/dt².
    dotH = -H²*eps, ddotH = H³(2*eps² - eps_prime).
    """
    dotH = -H**2*eps
    ddotH = H**3 * (2.0*eps**2 - eps_prime)
    
    a2, a4, a6 = a**2, a**4, a**6
    k2, k4, k6 = k**2, k**4, k**6
    k8, k10 = k**8, k**10
    H2, H4, H6 = H**2, H**4, H**6
    m2, m4, m6 = m**2, m**4, m**6
    
    m2pH2 = m2 + H2
    m2p3H2 = m2 + 3.0*H2
    m2p3H2m_dH = m2p3H2 - dotH
    m2p3H2p2dH = m2p3H2 + 2.0*dotH
    
    P = (4.0*(m2p3H2+3.0*dotH)*k4
         + 12.0*a2*(m2pH2*m2p3H2 + 2.0*m2pH2*dotH - dotH**2)*k2
         + 9.0*a4*m2pH2*m2p3H2m_dH*m2p3H2p2dH)
    if P <= 1e-60:
        return 1e-100, 1e-100
    
    # K_B
    KB = a4/P*(-4.0*dotH**2*k6 + 3.0*a2*m2pH2*m2p3H2m_dH*m2p3H2p2dH*k4)
    
    # M_B
    c10 = (12.0*m2pH2*m2p3H2**3
           + 16.0*m2p3H2**2*(6.0*m2+7.0*H2)*dotH
           + 4.0*m2p3H2*(63.0*m2+71.0*H2)*dotH**2
           + 8.0*(25.0*m2+27.0*H2)*dotH**3 - 48.0*dotH**4
           - 32.0*H*m2p3H2*dotH*ddotH - 48.0*H*dotH**2*ddotH)
    
    brkt = (2.0*m2pH2*m2p3H2**2*(2.0*m2+5.0*H2)
            + m2p3H2*(19.0*m4+64.0*m2*H2+49.0*H4)*dotH
            + 2.0*(7.0*m4+20.0*m2*H2+17.0*H4)*dotH**2
            - (23.0*m2+25.0*H2)*dotH**3 + 2.0*dotH**4
            - 2.0*H*m2pH2*m2p3H2*ddotH - 4.0*H*m2pH2*dotH*ddotH)
    c8 = 12.0*a2*m2p3H2p2dH*brkt
    
    c6 = (9.0*a4*m2pH2*m2p3H2m_dH*m2p3H2p2dH**2
          * (7.0*m2pH2*m2p3H2 + 17.0*m2pH2*dotH - 8.0*dotH**2))
    c4 = 27.0*a6*m2pH2**2*m2p3H2m_dH**2*m2p3H2p2dH**3
    
    MB = a6/P**2*(c10*k10 + c8*k8 + c6*k6 + c4*k4)
    return KB, MB


SECTORS = {
    'tensor_minimal':   O2_tensor_minimal,
    'vector_minimal':   O2_vector_minimal,
    'tensor_nonmin':    O2_tensor_nonmin,
    'vector_nonmin':    O2_vector_nonmin,
    'scalar_nonmin':    O2_scalar_nonmin,
}

# ============================================================
# Fast mode integration using exact-step (symplectic) integrator
# ============================================================

def integrate_mode_fast(k, m, bg_interp, sector, N_start, N_max, dN_step=0.002,
                        A_threshold=0.01, N_buffer=0.7):
    """
    Fast integration using exact-step (symplectic) method:
    u_{n+1} = u_n cos(Q_n dN) + u'_n sin(Q_n dN)/Q_n
    u'_{n+1} = -u_n Q_n sin(Q_n dN) + u'_n cos(Q_n dN)
    """
    a_s, H_s, eps_s, eps_p_s = bg_interp(N_start)
    eta_s = bg_interp.eta_at(N_start)
    aH_s = a_s * H_s
    
    # Bunch-Davies IC in chi -> u = chi * sqrt(aH)
    chi0 = np.exp(-1j*k*eta_s)/np.sqrt(2.0*k)
    chi_N0 = -1j*k*chi0/aH_s
    F0 = np.sqrt(aH_s)
    ur = chi0.real*F0; ui = chi0.imag*F0
    uNr = chi_N0.real*F0 + chi0.real*F0*0.5*(1.0-eps_s)
    uNi = chi_N0.imag*F0 + chi0.imag*F0*0.5*(1.0-eps_s)
    
    O2_func = SECTORS[sector]
    N_end = bg_interp.N_end
    
    N = N_start
    step_count = 0
    while N < N_max - 1e-10:
        # Adaptive step based on Q at current position
        a, H, eps, eps_p = bg_interp(N)
        Q_cur = np.sqrt(max(O2_func(k, a, H, eps, eps_p, m)/(a**2*H**2)
                            + 0.5*eps_p - 0.25*(1.0-eps)**2, 1e-60))
        dN = min(0.03, 0.5/max(Q_cur, 0.1))
        if N + dN > N_max:
            dN = N_max - N
        
        # Evaluate Q at midpoint for second-order accuracy
        N_mid = N + 0.5*dN
        a_mid, H_mid, eps_mid, eps_p_mid = bg_interp(N_mid)
        Q2_mid = max(O2_func(k, a_mid, H_mid, eps_mid, eps_p_mid, m)
                     / (a_mid**2*H_mid**2)
                     + 0.5*eps_p_mid - 0.25*(1.0-eps_mid)**2, 1e-60)
        Q_mid = np.sqrt(Q2_mid)
        
        cQ, sQ = np.cos(Q_mid*dN), np.sin(Q_mid*dN)
        ur_n = ur*cQ + uNr*sQ/Q_mid
        ui_n = ui*cQ + uNi*sQ/Q_mid
        uNr_n = -ur*Q_mid*sQ + uNr*cQ
        uNi_n = -ui*Q_mid*sQ + uNi*cQ
        
        ur, ui, uNr, uNi = ur_n, ui_n, uNr_n, uNi_n
        N += dN
        step_count += 1
    
    # Final state at N_max
    a_f, H_f, eps_f, eps_p_f = bg_interp(N_max)
    aH_f = a_f*H_f; F_f = np.sqrt(aH_f)
    chi_r = ur/F_f; chi_i = ui/F_f
    fp_f = 0.5*(1.0-eps_f)
    chi_N_r = uNr/F_f - ur*fp_f/F_f
    chi_N_i = uNi/F_f - ui*fp_f/F_f
    
    O2_f = O2_func(k, a_f, H_f, eps_f, eps_p_f, m)
    Om_f = np.sqrt(max(O2_f/(aH_f**2), 1e-60))
    
    beta_sq = 0.5*aH_f*Om_f*(chi_r**2+chi_i**2) + 0.5*aH_f/Om_f*(chi_N_r**2+chi_N_i**2) - 0.5
    wron = -2.0*aH_f*(chi_r*chi_N_i - chi_i*chi_N_r)
    
    # Adiabaticity check at N_max
    Q2_f = max(O2_f/(aH_f**2) + 0.5*eps_p_f - 0.25*(1.0-eps_f)**2, 1e-60)
    Q_f = np.sqrt(Q2_f)
    a_p = a_f*np.exp(0.01); H_p = H_f*np.exp(-eps_f*0.01)
    a_m = a_f*np.exp(-0.01); H_m = H_f*np.exp(eps_f*0.01)
    Q2_p = max(O2_func(k, a_p, H_p, eps_f, eps_p_f, m)/(a_p**2*H_p**2) + 0.5*eps_p_f - 0.25*(1.0-eps_f)**2, 1e-60)
    Q2_m = max(O2_func(k, a_m, H_m, eps_f, eps_p_f, m)/(a_m**2*H_m**2) + 0.5*eps_p_f - 0.25*(1.0-eps_f)**2, 1e-60)
    dQ_dN = abs(np.sqrt(Q2_p) - np.sqrt(Q2_m))/0.02
    A_val = dQ_dN/max(Q_f**2, 1e-60)
    
    reason = 'adiabatic' if A_val < A_threshold else 'not adiabatic'
    
    if abs(wron) > 1e-6:
        beta_sq_final = max((beta_sq + 0.5)/wron - 0.5, 0.0)
    else:
        beta_sq_final = max(beta_sq, 0.0)
    
    return {
        'k': k, 'm': m, 'sector': sector, 'beta_sq': beta_sq_final,
        'wronskian': wron, 'N_stop': N_max,
        'stop_reason': reason, 'Q_final': Q_f, 'A_final': A_val,
        'n_steps': step_count,
    }


# ============================================================
# Spectrum computation
# ============================================================

def compute_spectrum(bg, m_over_sqrt2_Hinf, sector='tensor_minimal',
                     k_min=0.01, k_max=500.0, n_k=80,
                     N_start_offset=3.0, N_extra=5.0,
                     A_threshold=0.01, N_buffer=0.7, verbose=True):
    """Compute |beta|^2 spectrum."""
    H_inf = bg['H_inf']; a_e = bg['a_e']; H_e = bg['H_e']
    m = m_over_sqrt2_Hinf * np.sqrt(2.0) * H_inf
    N_end = bg['N_end']; N_start = N_end - N_start_offset; N_max = N_end + N_extra
    bg_interp = BgInterp(bg)
    
    k_norm_arr = np.logspace(np.log10(k_min), np.log10(k_max), n_k)
    k_arr = k_norm_arr * a_e * H_e
    
    beta_sq_arr = np.zeros(n_k)
    stop_N_arr = np.zeros(n_k)
    reasons = []; wron_arr = np.zeros(n_k)
    
    if verbose:
        print(f"\n  {sector}  m/(√2 H_inf)={m_over_sqrt2_Hinf:.1f}  "
              f"k: {k_min:.2e}→{k_max:.2e}  n={n_k}")
    
    t0 = time.time()
    for i, (kv, kn) in enumerate(zip(k_arr, k_norm_arr)):
        res = integrate_mode_fast(kv, m, bg_interp, sector, N_start, N_max,
                                  A_threshold=A_threshold, N_buffer=N_buffer)
        beta_sq_arr[i] = res['beta_sq']; stop_N_arr[i] = res['N_stop']
        reasons.append(res['stop_reason']); wron_arr[i] = res['wronskian']
        
        if verbose and (i % 20 == 0 or i == n_k - 1):
            nk = kn**3/(2.0*np.pi**2)*res['beta_sq']
            wr_ok = "OK" if abs(res['wronskian']-1.0) < 0.2 else "BAD"
            print(f"    [{i:3d}/{n_k}] k={kn:.3e}  |β|²={res['beta_sq']:.3e}  "
                  f"n_k={nk:.3e}  N={res['N_stop']:.1f}({res['stop_reason']})  "
                  f"Wr={res['wronskian']:.4f} {wr_ok}  steps={res['n_steps']}  "
                  f"t={time.time()-t0:.1f}s")
    
    elapsed = time.time() - t0
    if verbose:
        print(f"    Done in {elapsed:.1f}s ({elapsed/n_k:.1f}s/mode)")
    
    nk_arr = k_norm_arr**3/(2.0*np.pi**2)*beta_sq_arr
    return {
        'k_norm': k_norm_arr, 'beta_sq': beta_sq_arr, 'nk': nk_arr,
        'm': m, 'm_over_sqrt2_Hinf': m_over_sqrt2_Hinf, 'sector': sector,
        'stop_N': stop_N_arr, 'stop_reason': reasons, 'wronskian': wron_arr,
    }


# ============================================================
# Summary / output
# ============================================================

def print_summary(spectra_dict):
    print("\n" + "=" * 70); print("SUMMARY"); print("=" * 70)
    for key, spec in spectra_dict.items():
        beta = spec['beta_sq']; kn = spec['k_norm']; nk = spec['nk']
        valid = beta > 1e-60
        if not valid.any():
            print(f"\n  {key}: no valid data"); continue
        idx = np.argmax(nk[valid]); wr = spec['wronskian']
        n_adiab = sum(1 for r in spec['stop_reason'] if r == 'adiabatic')
        n_tot = np.trapezoid(nk[valid], np.log(kn[valid]))
        print(f"\n  {key}:")
        print(f"    Peak k={kn[valid][idx]:.4f}  n_k={nk[valid][idx]:.4e}  |β|²={beta[valid][idx]:.4e}")
        print(f"    Adiabatic: {n_adiab}/{len(kn)}  Wr: μ={np.mean(wr):.4f} σ={np.std(wr):.4f}")
        print(f"    ∫ n_k d ln k = {n_tot:.3e}")


def save_output(spectra_dict, bg, outdir):
    os.makedirs(outdir, exist_ok=True)
    fname = os.path.join(outdir, 'cgpp_full_result.pkl')
    save_spec = {}
    for key, spec in spectra_dict.items():
        save_spec[key] = {k: spec[k] for k in
                          ['k_norm', 'beta_sq', 'nk', 'm', 'm_over_sqrt2_Hinf', 'sector']}
    save_bg = {k: bg[k] for k in ['N_end', 'a_e', 'H_e', 'H_inf']}
    with open(fname, 'wb') as f:
        pickle.dump({'spectra': save_spec, 'background': save_bg}, f)
    print(f"\nSaved to {fname}")


# ============================================================
# Quick / Full runs
# ============================================================

def run_quick():
    print("=" * 70); print("CGPP Full Numerics v2: Quick Verification"); print("=" * 70)
    bg = solve_background(N_extra=5.0)
    bg_int = BgInterp(bg)
    N_end = bg['N_end']; N_start = N_end - 3.0; N_max = N_end + 5.0
    H_inf = bg['H_inf']; a_e = bg['a_e']; H_e = bg['H_e']
    
    k_norms = [0.1, 1.0, 10.0, 100.0]
    m_val = 3.0*np.sqrt(2.0)*H_inf
    
    for sector in ['tensor_minimal', 'vector_minimal', 'tensor_nonmin', 'vector_nonmin']:
        print(f"\n[{sector}]  m/(√2 H_inf)=3.0:")
        for kn in k_norms:
            k = kn*a_e*H_e
            res = integrate_mode_fast(k, m_val, bg_int, sector, N_start, N_max)
            nk = kn**3/(2.0*np.pi**2)*res['beta_sq']
            print(f"  k/(a_eH_e)={kn:7.2f}  |β|²={res['beta_sq']:.4e}  "
                  f"n_k={nk:.4e}  Wr={res['wronskian']:.4f}  "
                  f"{res['stop_reason']}  steps={res['n_steps']}")
    print("\nDone.")


def run_full():
    print("=" * 70); print("CGPP Full Numerics v3: All Sectors"); print("=" * 70)
    bg = solve_background(N_extra=5.0)  # a/a_e ~ 148
    masses = [2.0, 3.0, 4.0, 5.0]
    sectors = ['tensor_minimal', 'vector_minimal', 'tensor_nonmin', 'vector_nonmin']
    
    spectra = {}
    for sector in sectors:
        for mv in masses:
            key = f"{sector}_m{mv:.1f}"
            try:
                spectra[key] = compute_spectrum(
                    bg, mv, sector=sector,
                    k_min=0.03, k_max=80.0, n_k=100,
                    N_start_offset=3.0, N_extra=5.0,
                    A_threshold=0.01, N_buffer=0.7, verbose=True)
            except Exception as e:
                print(f"    ERROR {key}: {e}")
                import traceback; traceback.print_exc()
    
    outdir = '/root/Agents/Spin2CGPP/Code/output'
    save_output(spectra, bg, outdir)
    print_summary(spectra)
    print("\nAll done.")


if __name__ == '__main__':
    if '--quick' in sys.argv or '-q' in sys.argv:
        run_quick()
    elif '--full' in sys.argv or '-f' in sys.argv:
        run_full()
    else:
        run_quick()

