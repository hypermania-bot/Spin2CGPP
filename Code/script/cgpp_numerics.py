#!/usr/bin/env python3
"""
CGPP Numerics: Direct numerical integration of mode equations in N-time.

Solves chi'' + (1-eps) chi' + Omega^2 chi = 0 in e-fold time N = ln(a),
with Bunch-Davies initial conditions.  Extracts the Bogoliubov coefficient
|beta|^2 using the adiabatic invariant when the mode has stabilized.

The integration can stop early when the adiabaticity parameter
A = |Omega'/Omega^2| falls below a threshold, indicating |beta|^2 has
converged.

Sectors implemented:
  - tensor_minimal : omega_k^2 = k^2 + a^2 m^2 - a^2 H^2 (2 - eps)
  - vector_minimal : omega_k^2 = k^2 + a^2 m^2 - f''/f,
                      f = a^2 / sqrt(k^2 + a^2 m^2)
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import os
import pickle
import time
import sys

# ============================================================
# Physical parameters (Planck units M_pl = 1)
# ============================================================
Mpl_GeV = 2.435e18
mphi_GeV = 4.14e12
mphi = mphi_GeV / Mpl_GeV
v_pl = 0.5

# ============================================================
# Hilltop inflation potential
# ============================================================

def V_hill(phi):
    x = phi / v_pl
    return mphi**2 * v_pl**2 / 72.0 * (1.0 - x**6)**2


def dV_hill(phi):
    x = phi / v_pl
    return -mphi**2 * v_pl / 6.0 * (1.0 - x**6) * x**5


# ============================================================
# Background evolution in e-fold time N
# ============================================================

def solve_background(N_extra=5.0):
    """
    Solve hilltop inflation background in e-fold time N.

    Integrates full hilltop ODE from N=0 (CMB) to N=N_end.
    Post-inflation background uses analytic continuation:
      eps(N) = 1  for N >= N_end  (kinetic-domination-like)
      H(N) = H_end * exp(-(N - N_end))
      a(N) = a_end * exp(N - N_end)

    This avoids the stiff oscillator phase and is always well-behaved.
    """
    bg_rtol = 1e-8
    bg_atol = 1e-12

    # phi0 giving N_cmb ~ 60 for this hilltop model
    phi0 = 0.04799

    psi0 = abs(dV_hill(phi0) / V_hill(phi0))
    H0 = np.sqrt(V_hill(phi0) / (3.0 - 0.5 * psi0**2))
    lnH0 = np.log(H0)

    def bg_ode(N, y):
        phi, psi, lnH, _eta = y
        H = np.exp(np.clip(float(lnH), -60, 10))
        eps = 0.5 * psi**2
        dpsi = -(3.0 - eps) * psi - dV_hill(phi) / H**2
        dlnH = -eps
        deta = 1.0 / (np.exp(N) * H)
        return [psi, dpsi, dlnH, deta]

    def ev_eps1(N, y):
        return 0.5 * y[1]**2 - 1.0
    ev_eps1.terminal = True
    ev_eps1.direction = 1

    sol = solve_ivp(bg_ode, [0, 80], [phi0, psi0, lnH0, 0.0],
                    events=[ev_eps1], method='RK45',
                    rtol=bg_rtol, atol=bg_atol)

    Ne = sol.t_events[0][0] if len(sol.t_events[0]) > 0 else sol.t[-1]

    # --- Analytic post-inflation continuation ---
    # For N >= N_end: use eps(N) that grows past 1 (matter-domination-like).
    # eps(N) = 1 + (eps_md - 1) * (1 - exp(-(N - N_end)/Delta))
    # with eps_md = 1.5 (matter domination) and Delta = 0.5 e-folds.
    #
    # This ensures eps > 1 for N > N_end, so aH decreases and k/(aH) grows,
    # allowing continued mode mixing and particle production.
    #
    # Numerical integration of d(ln H)/dN = -eps(N), a(N) = exp(N),
    # and d(eta)/dN = 1/(a H).

    # Build arrays for the inflationary part
    N_inf = sol.t
    phi_inf = sol.y[0]
    psi_inf = sol.y[1]
    lnH_inf = sol.y[2]
    H_inf_arr = np.exp(np.clip(lnH_inf, -60, 10))
    eps_inf_arr = 0.5 * psi_inf**2
    a_inf_arr = np.exp(N_inf)

    # End-of-inflation values
    idx_end = np.argmax(eps_inf_arr >= 1.0)
    if eps_inf_arr[idx_end] < 1.0:
        idx_end = len(eps_inf_arr) - 5
    N_end = N_inf[idx_end]
    a_e = a_inf_arr[idx_end]
    H_e = H_inf_arr[idx_end]

    # --- Eta: use dS formula during inflation to avoid catastrophic cancellation ---
    # During inflation: eta(N) = -1 / (a(N) * H(N))  [dS approximation, exact to O(eps)]
    # This avoids numerical issues from the integral being dominated by early times.
    eta_inf_dS = -1.0 / (a_inf_arr * H_inf_arr)

    # Post-inflation: eps transitions from 1 to eps_final
    eps_final = 1.5      # matter-domination-like
    Delta = 0.5           # transition width in e-folds

    N_post = np.linspace(N_end, N_end + N_extra,
                         max(40, int(N_extra * 20)))
    dN_fine = N_post[1] - N_post[0]
    a_post = np.exp(N_post)
    eps_post = 1.0 + (eps_final - 1.0) * (1.0 - np.exp(-(N_post - N_end) / Delta))

    # Integrate H and eta from end of inflation
    H_post = np.zeros(len(N_post))
    H_post[0] = H_e
    eta_post = np.zeros(len(N_post))
    eta_post[0] = -1.0 / (a_e * H_e)  # dS eta at N_end
    for i in range(1, len(N_post)):
        H_post[i] = H_post[i - 1] * np.exp(-eps_post[i - 1] * dN_fine)
        a_mid = a_post[i - 1]  # use beginning-of-step value for eta
        H_mid = H_post[i - 1]
        eta_post[i] = eta_post[i - 1] + dN_fine / (a_mid * H_mid)

    # Trim: exclude first point (N_end, already in inf array)
    N_post = N_post[1:]
    a_post = a_post[1:]
    H_post = H_post[1:]
    eps_post = eps_post[1:]
    eta_post = eta_post[1:]

    phi_post = np.full_like(N_post, phi_inf[idx_end])
    psi_post = np.full_like(N_post, np.sqrt(2.0 * eps_post))

    # Concatenate
    N_arr = np.concatenate([N_inf, N_post])
    a_arr = np.concatenate([a_inf_arr, a_post])
    H_arr = np.concatenate([H_inf_arr, H_post])
    eps_arr = np.concatenate([eps_inf_arr, eps_post])
    phi_arr = np.concatenate([phi_inf, phi_post])
    eta_arr = np.concatenate([eta_inf_dS, eta_post])

    eta_e = -1.0 / (a_e * H_e)

    # --- Add epsilon derivative eps' = d(eps)/dN ---
    eps_prime_arr = np.zeros_like(N_arr)
    eps_prime_arr[1:-1] = (eps_arr[2:] - eps_arr[:-2]) / (N_arr[2:] - N_arr[:-2])
    eps_prime_arr[0] = eps_prime_arr[1]
    eps_prime_arr[-1] = eps_prime_arr[-2]

    H_inf_est = np.mean(H_inf_arr[:max(10, len(H_inf_arr) // 20)])

    bg = {
        'N': N_arr,
        'a': a_arr,
        'H': H_arr,
        'eps': eps_arr,
        'eps_prime': eps_prime_arr,
        'phi': phi_arr,
        'eta': eta_arr,
        'N_end': N_end,
        'a_e': a_e,
        'H_e': H_e,
        'eta_e': eta_e,
        'H_inf': H_inf_est,
    }

    print(f"  N_end = {N_end:.2f}")
    print(f"  H_inf = {H_inf_est*Mpl_GeV:.3e} GeV")
    print(f"  H_e   = {H_e*Mpl_GeV:.3e} GeV  (ratio H_inf/H_e = {H_inf_est/H_e:.3f})")
    print(f"  m_phi / H_inf = {mphi/H_inf_est:.1f}")
    print(f"  eta_e = {bg['eta_e']:.3e}")
    return bg


# ============================================================
# Omega^2 functions
# ============================================================

def _omega2_tensor_minimal(k, a, H, eps, m):
    """omega_k^2 = k^2 + a^2 m^2 - a^2 H^2 (2 - eps)."""
    return k**2 + a**2 * (m**2 - H**2 * (2.0 - eps))


def _omega2_vector_minimal(k, a, H, eps, m):
    """
    omega_k^2 = k^2 + a^2 m^2 - f''/f,
    f = a^2 / sqrt(k^2 + a^2 m^2).

    Derived: f''/f = a^2 H^2 * P / D^2
    P = 6 k^4 + 5 k^2 a^2 m^2 + 2 a^4 m^4
        - eps (2 k^4 + 3 k^2 a^2 m^2 + a^4 m^4)
    D = k^2 + a^2 m^2
    """
    a2 = a**2
    k2 = k**2
    m2 = m**2
    D = k2 + a2 * m2
    a2H2 = a2 * H**2
    a2m2 = a2 * m2
    num = (6.0 * k2**2 + 5.0 * k2 * a2m2 + 2.0 * a2m2**2
           - eps * (2.0 * k2**2 + 3.0 * k2 * a2m2 + a2m2**2))
    fpp_over_f = a2H2 * num / D**2
    return k2 + a2m2 - fpp_over_f


def _omega2_tensor_nonmin(k, a, H, eps, m):
    """omega_k^2 = k^2 + a^2 m^2 + a^2 H^2 (1 - eps).  Eq. (4.15) with Lambda=0."""
    return k**2 + a**2 * (m**2 + H**2 * (1.0 - eps))


def _omega2_vector_nonmin(k, a, H, eps, m):
    """
    Non-minimal vector: omega_k^2 = (4 K_C M_C + K_C'^2 - 2 K_C K_C'')/(4 K_C^2)
    Evaluated via finite differences of K_C.
    """
    def _KC(av, Hv, ev):
        mu12 = m**2 + Hv**2 * (3.0 + ev)
        Dv = k**2 + av**2 * mu12
        if Dv <= 1e-300 or mu12 <= 0:
            return 1e-200
        return av**4 * k**2 * mu12 / Dv

    KC0 = _KC(a, H, eps)
    if KC0 <= 1e-200:
        return k**2

    dN = 0.001
    # Use nearby points for finite difference -- need eps values, assume eps varies slowly
    KCp = _KC(a * np.exp(dN), H * np.exp(-eps * dN), eps)
    KCm = _KC(a * np.exp(-dN), H * np.exp(eps * dN), eps)

    aH = a * H
    KCprime = aH * (KCp - KCm) / (2.0 * dN)
    KCpp = aH**2 * ((KCp - 2.0*KC0 + KCm)/dN**2 + (1.0 - eps)*(KCp - KCm)/(2.0*dN))

    MC0 = a**4 * k**2 * (m**2 + H**2*(3.0 - 2.0*eps))

    num = 4.0*KC0*MC0 + KCprime**2 - 2.0*KC0*KCpp
    den = 4.0*KC0**2
    if den <= 1e-200:
        return k**2
    return num / den


def _omega2_scalar_nonmin(k, a, H, eps, m):
    """
    Non-minimal scalar: omega_k^2 = (4 K_B M_B + K_B'^2 - 2 K_B K_B'')/(4 K_B^2)
    Evaluated via finite differences.
    Uses the full K_B/M_B polynomial formulas.
    """
    dN_fd = 0.001
    KB0, MB0 = _compute_KB_MB(k, a, H, eps, m)
    if KB0 <= 1e-200:
        return k**2

    KBp, _ = _compute_KB_MB(k, a*np.exp(dN_fd), H*np.exp(-eps*dN_fd), eps, m)
    KBm, _ = _compute_KB_MB(k, a*np.exp(-dN_fd), H*np.exp(eps*dN_fd), eps, m)

    L0 = np.log(KB0)
    Lp = np.log(max(KBp, 1e-200)); Lm = np.log(max(KBm, 1e-200))
    L_dN = (Lp - Lm) / (2.0*dN_fd)
    L_d2N = (Lp - 2.0*L0 + Lm) / dN_fd**2

    a2H2 = a**2 * H**2
    Omega2 = MB0/(a2H2*KB0) - 0.5*L_d2N - 0.5*(1.0-eps)*L_dN - 0.25*L_dN**2
    return max(Omega2*a2H2, 1e-200)


def _compute_KB_MB(k, a, H, eps, m):
    """Compute K_B and M_B for non-minimal scalar. See paper Eq. (4.32)-(4.33)."""
    # dotH = dH/dt = H dH/dN = -H^2 eps  (where d/dN gives d/dt = H d/dN)
    # ddotH = d^2 H/dt^2 = H^3(2 eps^2)  (approximation, ignoring eps')
    dotH = -H**2 * eps
    ddotH = H**3 * 2.0 * eps**2  # approximation

    a2, a4, a6 = a**2, a**4, a**6
    k2, k4, k6 = k**2, k**4, k**6
    k8, k10 = k**8, k**10
    H2, H4, H6 = H**2, H**4, H**6
    m2, m4, m6 = m**2, m**4, m**6

    m2pH2 = m2 + H2
    m2p3H2 = m2 + 3.0*H2
    m2p3H2m_dH = m2p3H2 - dotH
    m2p3H2p2dH = m2p3H2 + 2.0*dotH

    P = (4.0*(m2p3H2 + 3.0*dotH)*k4
         + 12.0*a2*(m2pH2*m2p3H2 + 2.0*m2pH2*dotH - dotH**2)*k2
         + 9.0*a4*m2pH2*m2p3H2m_dH*m2p3H2p2dH)
    if P <= 1e-200:
        return 1e-200, 1e-200

    KB = a4/P * (-4.0*dotH**2*k6 + 3.0*a2*m2pH2*m2p3H2m_dH*m2p3H2p2dH*k4)

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

    MB = a6/P**2 * (c10*k10 + c8*k8 + c6*k6 + c4*k4)
    return KB, MB


OMEGA2_FUNC = {
    'tensor_minimal': _omega2_tensor_minimal,
    'vector_minimal': _omega2_vector_minimal,
    'tensor_nonmin': _omega2_tensor_nonmin,
    'vector_nonmin': _omega2_vector_nonmin,
    'scalar_nonmin': _omega2_scalar_nonmin,
}


def compute_omega2(k, a, H, eps, m, sector):
    return OMEGA2_FUNC[sector](k, a, H, eps, m)


# ============================================================
# Background evaluator
# ============================================================

class BackgroundEval:
    """Fast evaluation of a, H, eps, eps_prime, eta at arbitrary N."""

    def __init__(self, bg):
        self.N_arr = bg['N']
        self.a_arr = bg['a']
        self.H_arr = bg['H']
        self.eps_arr = bg['eps']
        self.eps_prime_arr = bg['eps_prime']
        self.eta_arr = bg['eta']
        self.N_end = bg['N_end']

    def __call__(self, N):
        """Return (a, H, eps, eps_prime) interpolated at N."""
        i = np.searchsorted(self.N_arr, N)
        if i <= 0:
            return (float(self.a_arr[0]), float(self.H_arr[0]),
                    float(self.eps_arr[0]), float(self.eps_prime_arr[0]))
        if i >= len(self.N_arr):
            return (float(self.a_arr[-1]), float(self.H_arr[-1]),
                    float(self.eps_arr[-1]), float(self.eps_prime_arr[-1]))
        t = (N - self.N_arr[i - 1]) / (self.N_arr[i] - self.N_arr[i - 1])
        a = self.a_arr[i - 1] + t * (self.a_arr[i] - self.a_arr[i - 1])
        H = self.H_arr[i - 1] + t * (self.H_arr[i] - self.H_arr[i - 1])
        eps = self.eps_arr[i - 1] + t * (self.eps_arr[i] - self.eps_arr[i - 1])
        eps_p = self.eps_prime_arr[i - 1] + t * (self.eps_prime_arr[i] -
                                                  self.eps_prime_arr[i - 1])
        return float(a), float(H), float(eps), float(eps_p)

    def eta_at(self, N):
        """Conformal time at N."""
        i = np.searchsorted(self.N_arr, N)
        if i <= 0:
            return float(self.eta_arr[0])
        if i >= len(self.N_arr):
            return float(self.eta_arr[-1])
        t = (N - self.N_arr[i - 1]) / (self.N_arr[i] - self.N_arr[i - 1])
        return float(self.eta_arr[i - 1] + t *
                     (self.eta_arr[i] - self.eta_arr[i - 1]))


# ============================================================
# Adiabaticity helper
# ============================================================

def _compute_adiabaticity(k, a, H, eps, m, sector):
    """Compute A = |dOmega/dN| / Omega^2 using finite differences."""
    O2 = compute_omega2(k, a, H, eps, m, sector)
    aH2 = a**2 * H**2
    O2_dimless = max(O2 / aH2, 1e-60)
    Omega = np.sqrt(O2_dimless)

    dN = 1e-5
    a_p = a * np.exp(dN)
    H_p = H * np.exp(-eps * dN)
    O2_p = compute_omega2(k, a_p, H_p, eps, m, sector)
    a_m = a * np.exp(-dN)
    H_m = H * np.exp(eps * dN)
    O2_m = compute_omega2(k, a_m, H_m, eps, m, sector)
    dO2_dN = (O2_p - O2_m) / (2.0 * dN)
    dO2dim_dN = (dO2_dN - 2.0 * (1.0 - eps) * O2) / aH2
    dOmega_dN = dO2dim_dN / (2.0 * Omega)
    return float(abs(dOmega_dN) / max(Omega**2, 1e-60))


# ============================================================
# Mode integration
# ============================================================

def integrate_mode(k, m, bg_eval, sector, N_start, N_max,
                   A_threshold=0.01, N_buffer=0.7,
                   rtol=1e-10, atol=1e-14):
    """
    Integrate one k-mode from N_start to N_max.

    Uses the Hamiltonian variable u = chi * sqrt(aH), which satisfies
    u'' + Q^2 u = 0  (no friction term), where
    Q^2 = Omega^2 + eps'/2 - (1-eps)^2/4.

    This transformation gives much better Wronskian conservation
    than integrating chi directly in N-time.

    After integration, the physical chi is recovered as chi = u / sqrt(aH),
    and |beta|^2 is computed from the adiabatic invariant.
    """
    a_start, H_start, eps_start, _ = bg_eval(N_start)
    eta_start = bg_eval.eta_at(N_start)
    aH_start = a_start * H_start

    # Bunch-Davies in chi variables
    chi0 = np.exp(-1j * k * eta_start) / np.sqrt(2.0 * k)
    chi_N0 = -1j * k * chi0 / aH_start

    # Transform to u = chi * sqrt(aH)
    F0 = np.sqrt(aH_start)
    u0 = chi0 * F0
    # u' = chi' * F + chi * F'
    # F' = d(sqrt(aH))/dN = (1/(2F)) * d(aH)/dN = (1-eps) * F / 2
    F_prime_factor = 0.5 * (1.0 - eps_start)
    u_N0 = chi_N0 * F0 + chi0 * F0 * F_prime_factor

    y0 = np.array([u0.real, u0.imag, u_N0.real, u_N0.imag])

    # Storage for tracking (beta only; adiabaticity computed once at end)
    traj_N = []
    traj_beta = []
    traj_wron = []

    def mode_ode(N, y):
        u_r, u_i, u_N_r, u_N_i = y
        a, H, eps, eps_prime = bg_eval(N)
        O2 = compute_omega2(k, a, H, eps, m, sector)
        O2_dimless = O2 / (a**2 * H**2)
        # Q^2 = Omega^2 + eps'/2 - (1-eps)^2/4
        Q2 = max(O2_dimless + 0.5 * eps_prime - 0.25 * (1.0 - eps)**2, 1e-60)
        u_NN_r = -Q2 * u_r
        u_NN_i = -Q2 * u_i
        return [u_N_r, u_N_i, u_NN_r, u_NN_i]

    def record_state(N, y):
        """Record beta^2 and Wronskian at current state."""
        u_r, u_i, u_N_r, u_N_i = y
        a, H, eps, _ = bg_eval(N)
        aH = a * H
        F = np.sqrt(aH)

        chi_r = u_r / F
        chi_i = u_i / F
        Fp_over_F = 0.5 * (1.0 - eps)
        chi_N_r = u_N_r / F - u_r * Fp_over_F / F
        chi_N_i = u_N_i / F - u_i * Fp_over_F / F

        O2 = compute_omega2(k, a, H, eps, m, sector)
        O2_dimless = max(O2 / (aH**2), 1e-60)
        Omega = np.sqrt(O2_dimless)

        chi_abs2 = chi_r**2 + chi_i**2
        chi_N_abs2 = chi_N_r**2 + chi_N_i**2
        beta_sq = 0.5 * aH * Omega * chi_abs2 \
            + 0.5 * aH / Omega * chi_N_abs2 \
            - 0.5

        wron = -2.0 * aH * (chi_r * chi_N_i - chi_i * chi_N_r)

        traj_N.append(N)
        traj_beta.append(beta_sq)
        traj_wron.append(wron)

    # Integrate with specified output grid to decouple recording from solver steps
    N_end_integ = bg_eval.N_end
    N_grid_inf = np.linspace(N_start, N_end_integ - 0.5,
                              max(15, int((N_end_integ - 0.5 - N_start) * 5)))
    N_grid_trans = np.linspace(N_end_integ - 0.5, N_max,
                                max(30, int((N_max - N_end_integ + 0.5) * 20)))
    N_grid = np.unique(np.concatenate([N_grid_inf, N_grid_trans]))
    N_grid = N_grid[(N_grid >= N_start) & (N_grid <= N_max)]

    sol = solve_ivp(mode_ode, [N_start, N_max], y0,
                    method='RK45', rtol=rtol, atol=atol,
                    max_step=0.2, t_eval=N_grid)

    # Record at output times
    for i, N_val in enumerate(sol.t):
        record_state(N_val, sol.y[:, i])

    # Convert to arrays for post-processing
    traj_N = np.array(traj_N)
    traj_beta = np.array(traj_beta)
    traj_wron = np.array(traj_wron)

    # Find stop time: first time past N_end + N_buffer with A < threshold
    N_end = bg_eval.N_end
    mask_late = traj_N >= N_end + N_buffer
    if mask_late.any():
        late_indices = np.where(mask_late)[0]
        idx_stop = None
        for idx in late_indices:
            a, H, eps, _ = bg_eval(traj_N[idx])
            A_val = _compute_adiabaticity(k, a, H, eps, m, sector)
            if A_val < A_threshold:
                idx_stop = idx
                break
        if idx_stop is not None:
            stop_reason = 'adiabatic'
        else:
            idx_stop = len(traj_N) - 1
            stop_reason = 'reached N_max'
    else:
        idx_stop = len(traj_N) - 1
        stop_reason = 'reached N_max'

    N_stop = traj_N[idx_stop]
    wron_final = traj_wron[idx_stop]
    beta_raw = float(traj_beta[idx_stop])

    # Wronskian renormalization: |beta|^2_ren = (beta_raw + 0.5) / W - 0.5
    if abs(wron_final) > 1e-6:
        beta_sq_final = max((beta_raw + 0.5) / wron_final - 0.5, 0.0)
    else:
        beta_sq_final = max(beta_raw, 0.0)

    return {
        'k': k,
        'm': m,
        'sector': sector,
        'N_start': N_start,
        'N_stop': N_stop,
        'N_end': N_end,
        'beta_sq': beta_sq_final,
        'wronskian': wron_final,
        'stop_reason': stop_reason,
        'traj_N': traj_N,
        'traj_beta': traj_beta,
        'traj_wron': traj_wron,
    }


# ============================================================
# Spectrum computation
# ============================================================

def compute_spectrum(bg, m_over_sqrt2_Hinf, sector='tensor_minimal',
                     k_min=0.01, k_max=500.0, n_k=100,
                     N_start_offset=3.0, N_extra=5.0,
                     A_threshold=0.01, N_buffer=0.7,
                     verbose=True):
    """
    Compute |beta|^2 spectrum for a given sector and mass.

    k_norm = k / (a_e H_e) is the dimensionless wavenumber.
    """
    H_inf = bg['H_inf']
    a_e = bg['a_e']
    H_e = bg['H_e']
    m = m_over_sqrt2_Hinf * np.sqrt(2.0) * H_inf

    N_end = bg['N_end']
    N_start = N_end - N_start_offset
    N_max = N_end + N_extra

    bg_eval = BackgroundEval(bg)

    k_norm_arr = np.logspace(np.log10(k_min), np.log10(k_max), n_k)
    k_arr = k_norm_arr * a_e * H_e

    beta_sq_arr = np.zeros(n_k)
    stop_N_arr = np.zeros(n_k)
    stop_reason_arr = ['' for _ in range(n_k)]
    wron_arr = np.zeros(n_k)

    if verbose:
        print(f"\n  Sector: {sector}")
        print(f"    m/(sqrt(2) H_inf) = {m_over_sqrt2_Hinf:.2f}  "
              f"(m/H_inf = {m/H_inf:.3f})")
        print(f"    k range: {k_min:.2e} - {k_max:.2e}  ({n_k} modes)")
        print(f"    N: start={N_start:.2f}  end={N_end:.2f}  max={N_max:.2f}")
        print(f"    A_threshold={A_threshold}  N_buffer={N_buffer}")

    t0 = time.time()
    for i, (kv, k_norm) in enumerate(zip(k_arr, k_norm_arr)):
        res = integrate_mode(kv, m, bg_eval, sector,
                             N_start, N_max,
                             A_threshold=A_threshold,
                             N_buffer=N_buffer)
        beta_sq_arr[i] = res['beta_sq']
        stop_N_arr[i] = res['N_stop']
        stop_reason_arr[i] = res['stop_reason']
        wron_arr[i] = res['wronskian']

        if verbose and (i % 25 == 0 or i == n_k - 1):
            nk = k_norm**3 / (2.0 * np.pi**2) * res['beta_sq']
            wr_ok = "OK" if abs(res['wronskian'] - 1.0) < 0.2 else "BAD"
            elapsed = time.time() - t0
            print(f"    [{i:3d}/{n_k}] k={k_norm:.3e}  "
                  f"|beta|^2={res['beta_sq']:.3e}  "
                  f"n_k={nk:.3e}  "
                  f"stop={res['N_stop']:.2f}({res['stop_reason']})  "
                  f"Wr={res['wronskian']:.4f} {wr_ok}  "
                  f"t={elapsed:.1f}s")

    t_elapsed = time.time() - t0
    if verbose:
        print(f"    Completed in {t_elapsed:.1f}s "
              f"({t_elapsed/n_k:.1f}s per mode)")

    nk_arr = k_norm_arr**3 / (2.0 * np.pi**2) * beta_sq_arr

    return {
        'k_norm': k_norm_arr,
        'beta_sq': beta_sq_arr,
        'nk': nk_arr,
        'm': m,
        'm_over_sqrt2_Hinf': m_over_sqrt2_Hinf,
        'sector': sector,
        'stop_N': stop_N_arr,
        'stop_reason': stop_reason_arr,
        'wronskian': wron_arr,
    }


# ============================================================
# Summary and output
# ============================================================

def print_summary(spectra_dict, bg):
    """Print summary statistics."""
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    for key, spec in spectra_dict.items():
        beta = spec['beta_sq']
        kn = spec['k_norm']
        nk = spec['nk']

        valid = beta > 1e-50
        if not valid.any():
            print(f"\n  {key}: no valid data points")
            continue

        idx_peak = np.argmax(nk[valid])
        reasons = [spec['stop_reason'][j] for j in range(len(kn))]
        n_stopped = sum(1 for r in reasons if r == 'adiabatic')
        wrons = spec['wronskian']

        print(f"\n  {key}:")
        print(f"    Peak: k/(a_e H_e) = {kn[valid][idx_peak]:.4f}")
        print(f"    Peak n_k           = {nk[valid][idx_peak]:.4e}")
        print(f"    Peak |beta|^2      = {beta[valid][idx_peak]:.4e}")
        print(f"    Stopped adiabatic:   {n_stopped}/{len(kn)}")
        print(f"    Wronskian stats:     mean={np.mean(wrons):.4f}  "
              f"std={np.std(wrons):.4f}  "
              f"min={np.min(wrons):.4f}  max={np.max(wrons):.4f}")

        # Integrated number density per polarization: n = ∫ n_k d(ln k)
        n_tot = np.trapezoid(nk[valid], np.log(kn[valid]))
        print(f"    a^3 n/(a_e H_e)^3  = {n_tot:.3e}")


def save_output(spectra_dict, bg, outdir):
    """Save spectra and background summary to pickle."""
    os.makedirs(outdir, exist_ok=True)
    fname = os.path.join(outdir, 'cgpp_numerics_result.pkl')

    save_spec = {}
    for key, spec in spectra_dict.items():
        save_spec[key] = {
            'k_norm': spec['k_norm'],
            'beta_sq': spec['beta_sq'],
            'nk': spec['nk'],
            'm': spec['m'],
            'm_over_sqrt2_Hinf': spec['m_over_sqrt2_Hinf'],
            'sector': spec['sector'],
        }

    save_bg = {k: bg[k] for k in ['N_end', 'a_e', 'H_e', 'H_inf', 'eta_e']}

    with open(fname, 'wb') as f:
        pickle.dump({'spectra': save_spec, 'background': save_bg}, f)
    print(f"\nSaved to {fname}")


# ============================================================
# Quick test function
# ============================================================

def run_quick(rtol=1e-9, atol=1e-14):
    """Quick test with a single k-mode to verify correctness."""
    print("=" * 60)
    print("CGPP Numerics: Quick verification")
    print("=" * 60)

    print("\n[1] Background...")
    bg = solve_background(N_extra=5.0)

    bg_eval = BackgroundEval(bg)
    N_end = bg['N_end']
    N_start = N_end - 3.0
    N_max = N_end + 5.0

    m_over_sqrt2 = 3.0
    m = m_over_sqrt2 * np.sqrt(2.0) * bg['H_inf']
    a_e = bg['a_e']
    H_e = bg['H_e']

    # Test a few k values
    k_norms = [0.1, 1.0, 10.0, 100.0]
    print(f"\n[2] Testing m/(sqrt(2) H_inf) = {m_over_sqrt2} "
          f"(tensor minimal):")

    for k_norm in k_norms:
        k = k_norm * a_e * H_e
        res = integrate_mode(k, m, bg_eval, 'tensor_minimal',
                             N_start, N_max,
                             A_threshold=0.01, N_buffer=0.7,
                             rtol=rtol, atol=atol)
        nk = k_norm**3 / (2.0 * np.pi**2) * res['beta_sq']
        print(f"  k/(a_e H_e) = {k_norm:7.2f}  "
              f"|beta|^2 = {res['beta_sq']:.4e}  "
              f"n_k = {nk:.4e}  "
              f"stop = {res['N_stop']:.2f} ({res['stop_reason']})  "
              f"Wr = {res['wronskian']:.4f}")

        # Verify Wronskian stability
        wr_arr = res['traj_wron']
        wr_drift = (wr_arr[-1] - wr_arr[0]) / max(abs(wr_arr[0]), 1e-30)
        if abs(wr_drift) > 1e-2:
            print(f"    *** WARNING: Wronskian drift = {wr_drift:.2e} ***")

    print("\n[3] Testing vector minimal sector:")
    for k_norm in k_norms:
        k = k_norm * a_e * H_e
        res = integrate_mode(k, m, bg_eval, 'vector_minimal',
                             N_start, N_max,
                             A_threshold=0.01, N_buffer=0.7,
                             rtol=rtol, atol=atol)
        nk = k_norm**3 / (2.0 * np.pi**2) * res['beta_sq']
        print(f"  k/(a_e H_e) = {k_norm:7.2f}  "
              f"|beta|^2 = {res['beta_sq']:.4e}  "
              f"n_k = {nk:.4e}  "
              f"stop = {res['N_stop']:.2f} ({res['stop_reason']})  "
              f"Wr = {res['wronskian']:.4f}")

    print("\nDone.")


# ============================================================
# Full spectrum
# ============================================================

def run_full():
    """Compute full spectra for all sectors."""
    print("=" * 60)
    print("CGPP Numerics: Full spectrum")
    print("=" * 60)

    print("\n[1] Solving background...")
    bg = solve_background(N_extra=5.0)

    print("\n[2] Computing spectra...")
    masses = [2.0, 3.0, 4.0, 5.0]
    sectors = ['tensor_minimal', 'vector_minimal', 'tensor_nonmin', 'vector_nonmin']

    spectra = {}
    for sector in sectors:
        for mv in masses:
            key = f"{sector}_m{mv:.1f}"
            print(f"\n  --- {key} ---")
            try:
                spectra[key] = compute_spectrum(
                    bg, mv, sector=sector,
                    k_min=0.03, k_max=50.0, n_k=40,
                    N_start_offset=3.0, N_extra=5.0,
                    A_threshold=0.01, N_buffer=0.7,
                    verbose=True)
            except Exception as e:
                print(f"    ERROR: {e}")
                import traceback; traceback.print_exc()

    # Save first (robust against post-processing errors)
    outdir = '/root/Agents/Spin2CGPP/Code/output'
    save_output(spectra, bg, outdir)

    print_summary(spectra, bg)

    print("\nAll done.")


# ============================================================
# Main
# ============================================================

if __name__ == '__main__':
    if '--quick' in sys.argv or '-q' in sys.argv:
        run_quick()
    else:
        run_full()
