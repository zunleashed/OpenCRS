"""Ben-Tal Layer 1 ventilatory module (Ishbulatov 2022 + Ben-Tal & Smith 2008).

Replaces ONLY OpenCRS's ventilatory-drive block: the VAflow equation and
the breath-timing t1/t2 spline lookup in All_derivatives_njit.py.
Cardiovascular core, Chiari gas exchange and lung mechanics are untouched.

Primary sources (read from the papers, not from secondary summaries):
  - Ishbulatov et al. 2022, Mathematics 10:1088, sec 2.1.3, eqs 22-26,
    Table 1 (all values below verified against the MDPI PDF, 2026-11-15).
  - Ben-Tal & Smith 2008, J Theor Biol 251:480-497 (PMC2440661), eqs 11, 13,
    sec 4.1.3 (open loop: g_n = 13.2, K_n = 0.5 -> 5 s period) and sec 4.2.1
    Figs 16/18 (P vs PI controller on a 7.5% inspired-CO2 step).

Equations (Ishbulatov numbering):
  (22) v_c(t)  = k1c*tanh((p_ce(t - th_vc) - k2c)/k3c) + k4c
                central CO2 receptors; p_ce = end-capillary / brain PCO2.
                Here p_ce is OpenCRS PmbCO2 (brain PCO2), which already
                includes the blood-to-brain transport lag, so no explicit
                th_vc delay buffer is built (th_vc = 5 s, published).
  (23) v_o(t)  = k1o*tanh((p_oe(t) - k2o)/k3o) + k4o   (carotid O2)
                AS PRINTED the paper's eq 23 is piecewise: the sigmoidal
                branch for p_oe <= 100 mmHg and a linear branch
                k5o*p_oe + k4o above 100. That printed form is
                DISCONTINUOUS at 100 mmHg (sigmoidal: -0.004; linear:
                +0.491; jump +0.495, i.e. +30 in g_tilde via k_og = 60)
                and the >100 branch is physiologically backwards
                (hyperoxia driving ventilation). Layer 1 uses the
                sigmoidal branch for ALL PO2. With published k2o = -50
                the receptor is numerically silent (v_o ~= -0.004) across
                the whole physiological range, exactly as in the paper's
                own operating regime; k_og / k_oK are therefore
                unidentifiable from Protocol B / FiCO2 (documented).
  (24) g_tilde = g0 + k_cg*v_c + k_og*v_o + k_ex_g*max(0, MRV)
                frequency drive. k_ex_g is NEW (no paper contains an
                exercise feedforward); MRV is OpenCRS's existing lagged
                metabolic ratio (tau_MRV = 50 s), i.e. tau_ex = 50 s.
  (25) K       = K0 + k_cK*v_c + k_oK*v_o + I_k + k_ex_K*max(0, MRV)
                amplitude drive. k_ex_K is NEW, as above.
  (26) I_k     = trapezoidal integral of chemoreceptor activity.
                Implemented as the ODE dI_k/dt = k1_Ik*v_c + k2_Ik*v_o
                MINUS a leak -I_k/tau_Ik. The leak is a Layer 1
                adaptation, added on prototype evidence: the faithful
                cumulative integral winds up without bound under a
                sustained CO2 stimulus (prototype: I_k = +21 after
                660 s of FiCO2 8%, VT 5.1 L -- gate V1 fail), while
                tau_Ik = 60 s passes V1 (FiCO2 8%: RR 27, VT 2.45 L,
                VE 67 L/min; PI beats P on the 7.5% step, 51.6 vs 46.3).
                FITTED value is tau_Ik = 30 s (fit 1, Protocol B, set-C
                base): it halves the I_k windup that otherwise drives
                VT past VC at Z5. Trade-off, measured: at tau_Ik = 30 s
                with the fitted gains the FiCO2 8% response falls to
                VE 20.8 L/min (target ~51) -- see the eval-repo README,
                "Ben-Tal Layer 1 ventilatory module".

Map (option (a) in the scope doc) from (g_tilde, K) to (BF, VT):
  BF = BF_rest + s_BF*(g_tilde - g_tilde_rest)
  VT = VT_rest + s_VT*(K - K_rest)
  where g_tilde_rest / K_rest are evaluated at the rest point
  (PmbCO2 = 40, PamO2 = 100, MRV = 0, I_k = 0) so rest is exact by
  construction under any fitted gains. Scales are fixed by the
  operating-range argument, NOT fitted: g_tilde 16 -> ~39 gives
  BF 12 -> ~53/min; K 0.575 -> ~14 gives VT 0.5 -> ~2.5 L.
  VAflow follows by dead-space inversion:
  VAflow = BF*(VT - V0_dead) / (1 + GV_dead*BF), clamped to
  [0.04, VAflow_max]. The work-of-breathing spline is used ONLY for
  the inspiratory fraction t1/(t1+t2) at the computed VAflow
  (a spline-only map would collapse the two DOF into one).

Ceiling: VAflow_max = 3.5 L/s (fork hard ceiling, from the MEASURED
unclamped Z5 drive of 2.8-3.5 L/s). The scope doc's 2.2-2.5 L/s MVV
ceiling is superseded -- see the eval-repo README ceiling note.

State: I_k is a proper ODE state (83rd state, appended to the
respiratory-control group), initialised to 0.
"""
import math

import numpy as np
from numba import njit


@njit
def bental_layer1_drive(PmbCO2, PamO2, MRV, Ik,
                        k1c, k2c, k3c, k4c,
                        k1o, k2o, k3o, k4o,
                        g0, k_cg, k_og, K0, k_cK, k_oK,
                        k1_Ik, k2_Ik, tau_Ik,
                        k_ex_g, k_ex_K,
                        BF_rest, VT_rest, s_BF, s_VT,
                        GV_dead, V0_dead, VAflow_max):
    # Ishbulatov 2022 eq 22: central CO2 receptor (sigmoidal).
    # Input is OpenCRS PmbCO2 (brain PCO2); the 5 s th_vc blood-to-brain
    # delay is already in that quantity, so no delay buffer is built.
    v_c = k1c * math.tanh((PmbCO2 - k2c) / k3c) + k4c
    # Ishbulatov 2022 eq 23: carotid O2 receptor. Sigmoidal branch for
    # ALL PO2 -- see the module docstring: the printed piecewise form is
    # discontinuous at 100 mmHg (+0.495) and physiologically backwards.
    # With published k2o = -50 this is ~= -0.004 everywhere (silent).
    v_o = k1o * math.tanh((PamO2 - k2o) / k3o) + k4o
    # Exercise feedforward (NEW -- no Ben-Tal paper contains one):
    # MRV is OpenCRS's lagged metabolic ratio (tau_MRV = 50 s).
    ff = MRV if MRV > 0.0 else 0.0
    # Ishbulatov 2022 eq 24: frequency drive + feedforward.
    g_tilde = g0 + k_cg * v_c + k_og * v_o + k_ex_g * ff
    # Ishbulatov 2022 eq 25: amplitude drive + integral + feedforward.
    K = K0 + k_cK * v_c + k_oK * v_o + Ik + k_ex_K * ff
    # Ishbulatov 2022 eq 26 as an ODE, plus the Layer 1 anti-windup leak
    # (tau_Ik = 60 s; faithful integral = tau_Ik -> inf).
    dIk_dt = k1_Ik * v_c + k2_Ik * v_o - Ik / tau_Ik
    # Rest operating point (PmbCO2 = 40, PamO2 = 100, MRV = 0, Ik = 0),
    # recomputed from the current gains so rest stays exact when fitted.
    v_c_rest = k1c * math.tanh((40.0 - k2c) / k3c) + k4c
    v_o_rest = k1o * math.tanh((100.0 - k2o) / k3o) + k4o
    g_rest = g0 + k_cg * v_c_rest + k_og * v_o_rest
    K_rest = K0 + k_cK * v_c_rest + k_oK * v_o_rest
    # Option (a) map: g_tilde -> BF, K -> VT (scales fixed, not fitted).
    BF = BF_rest + s_BF * (g_tilde - g_rest)
    VT = VT_rest + s_VT * (K - K_rest)
    if BF < 0.05:
        BF = 0.05
    if VT < 0.10:
        VT = 0.10
    # Dead-space inversion: VT = (VAflow + BF*VD)/BF,
    # VD = GV_dead*VAflow + V0_dead.
    VAflow = BF * (VT - V0_dead) / (1.0 + GV_dead * BF)
    if VAflow < 0.04:
        VAflow = 0.04
    if VAflow > VAflow_max:
        VAflow = VAflow_max
    return v_c, v_o, g_tilde, K, dIk_dt, BF, VT, VAflow


# Published values (Ishbulatov 2022 Table 1) + Layer 1 design constants.
# Mirrored in Entire_system/fixed_params.py (the fork's parameter source).
BENTAL_DEFAULTS = dict(
    # Central CO2 receptor, eq 22.
    k1c=0.585, k2c=61.0, k3c=20.7, k4c=0.442,
    # Carotid O2 receptor, eq 23 (sigmoidal branch only; k5o/k6o unused).
    k1o=-0.5, k2o=-50.0, k3o=20.0, k4o=0.496,
    # Frequency / amplitude drives, eqs 24-25 (published gains).
    g0=16.0, k_cg=8.0, k_og=60.0,
    K0=0.75, k_cK=10.0, k_oK=7.5,
    # Integral weights, eq 26.
    k1_Ik=0.1, k2_Ik=0.1,
    # Layer 1 adaptations / new terms (see module docstring).
    tau_Ik=30.0,          # s, anti-windup leak (fit 1; V1 at 60 s)
    k_ex_g=0.8, k_ex_K=0.5,  # FITTED fit 1 (Protocol B, set-C base)
    # Option (a) map scales (fixed by operating-range argument).
    BF_rest=0.20, VT_rest=0.50, s_BF=0.03, s_VT=0.15,
)
# The six first-fit knobs (scope doc sec 4): k_ex_g, k_ex_K, k_cg, k_og,
# k_cK, k_oK. Everything else held at the published values above.
FIT_KNOBS = ("k_ex_g", "k_ex_K", "k_cg", "k_og", "k_cK", "k_oK")
