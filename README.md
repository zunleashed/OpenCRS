# OpenCRS

> **Fork note (zunleashed/OpenCRS, Zone 5 fork).** This is a fork of
> [Sheng-Ya/OpenCRS](https://github.com/Sheng-Ya/OpenCRS) at upstream commit
> `f0c3c4612f310017f07d5e67353d7a770034c396`, MIT licence (see `LICENSE`).
> Upstream is unchanged except for exactly four changes: the two to raise
> the ventilatory ceiling (below), plus the Ben-Tal Layer 1 ventilatory
> module (change 3), plus the dynamic muscle-pressure fix (change 4).
>
> 1. **VAflow clamp parameterised** — `Entire_system/DGSM_Union_Paper_gas/All_derivatives_njit.py:242`:
>    `VAflow = min(max(VAflow, 0.04), 1.6)` → `min(max(VAflow, 0.04), VAflow_max)`,
>    with new parameter `VAflow_max` (default **3.5 L/s**, added to
>    `Entire_system/fixed_params.py` and to
>    `Entire_system/DGSM_Union_Paper_gas/fixed_params.py`; threaded through
>    `Samples_for_DGSM_Union.py` as the last entry of `Input_Parameters` and
>    unpacked as the last entry in `All_derivatives_njit.njit_compatible`).
>    3.5 L/s was set by measurement: the unclamped ventilatory drive at 95%
>    VO2max was 2.8–3.5 L/s, so 3.5 covers it with a small margin. It is NOT
>    a physiological ceiling — a maximal-ventilation soft ceiling (~2.2–2.5
>    L/s) is a later refinement, not this change.
> 2. **Breath-timing spline grid extended** — `Entire_system/DGSM_Union_Paper_gas/Samples_for_DGSM_Union.py:124`
>    (`minimise_breathing`): `np.linspace(0.04, 1.6, 200)` →
>    `np.linspace(0.04, VAflow_max, 440)`, keeping the original ~0.0078 L/s
>    knot spacing; `VAflow_max` is passed in from `simulate_cpu`.
> 3. **Ben-Tal Layer 1 ventilatory module** —
>    `Entire_system/DGSM_Union_Paper_gas/Bental_Layer1.py` (new file) replaces
>    only the ventilatory-drive block in `All_derivatives_njit.py` (the
>    VAflow equation + breath-timing t1/t2 splines). Sigmoidal chemoreceptors
>    (Ishbulatov et al. 2022, *Mathematics* 10:1088, eqs 22–26) driving
>    separate frequency/amplitude outputs, plus a new exercise-feedforward
>    term (`k_ex`, no paper contains one) driven by OpenCRS's lagged
>    metabolic ratio MRV. State `Ik` added (83-state RHS). Two documented
>    adaptations of the published equations: the carotid-O2 receptor uses
>    the sigmoidal branch of eq 23 for all PO2 (the printed piecewise form
>    is discontinuous at PO2 = 100 mmHg and physiologically backwards
>    above it), and `Ik` has an anti-windup leak `−Ik/tau_Ik` (the faithful
>    cumulative integral winds up without bound under sustained CO2).
>    Fitted values (refit 2026-11-15, fixed plant, set-C base — best of 6
>    jointly scored candidates, joint score 0.239): `k_cg` 20.0,
>    `k_cK` 35.0, `k_ex_g` 0.8, `k_ex_K` 1.2, `tau_Ik` 30 s
>    (published Ishbulatov gains were 20 / 25 / — / — / —; fit 1 was
>    8 / 10 / 0.8 / 0.5 / 30). **Gate V2 is still NOT passed** — the
>    refit fixes FiCO2 8% (VE 50.4 vs ~51, VT 2.20 exact) and HR
>    grading (112→187 vs ACSM 129→181) and VE at Z5 (153 vs ~139),
>    but the paper test still fails (VE 32.3 vs 59.9) and the
>    VE/RR/VT shape is wrong (deep-slow breathing: VT 3.89 L at Z5
>    vs ~2.32, RR 39 vs ~60). Full results, sensitivity and the
>    coupling verdict: evaluation repo README, section "Layer 1 refit
>    on fixed plant".
>    Diagnostic stores `PAO2_store`, `PACO2_store`, `Pd5O2_store`,
>    `Pd5CO2_store` (alveolar + dead-space gas) were added alongside the
>    Layer 1 diagnostics for the V2 diagnosis; they are read-only
>    diagnostics, not model changes.
> 4. **Dynamic muscle pressure in the breath optimiser** —
>    `Resp_Control_Breath_Optimiser.py` (`compute_constants`, in both
>    copies of the file — the fork ships it at the repo root and under
>    `Entire_system/DGSM_Union_Paper_gas/`, and the harness imports the
>    root copy): the muscle-pressure quadratic coefficient `a2` is now
>    solved so the true solution of the equation of motion
>    `R_rs*dV/dt + E_rs*V = P_musc(t)`, `V(0) = 0`, reaches `Vt1` at `t1`
>    (closed form, `a2 = Vt1*R_rs/(I2 - t1^2*I0)`), instead of the old
>    quasi-static `Pt1 = P_ao + E_rs*Vt1`. The old choice ignored the
>    resistive drop, so the ODE solution fell short of `Vt1` by `d1`
>    (25% at Z5) and the `V(t)` formula patched the endpoint with a
>    spurious constant — the `dV/dt` fed to gas exchange was only 0.38x
>    the true volume rate at Z5, the dead-space cascade flushed 2.6x
>    slow, and effective alveolar ventilation saturated ~48 L/min
>    (PaCO2 → 91 mmHg). With the fix `d1 = 0` and `dV/dt` is the true
>    derivative of `V(t)` at rest, paper exercise and Z5. Effect (same
>    Layer 1 fit-1 knobs): Protocol B PaCO2 43.8–46.7 mmHg every zone
>    (was 54–91), VA_eff tracks nominal VA (92.9 vs 88.1 L/min at Z5,
>    was 47.6 vs 172.6) — but VE/HR now under-shoot (VE 106.6, HR 138
>    at Z5) because the fitted knobs were tuned against the broken
>    plant and leaned on the pathological PaCO2 signal for drive;
>    refit done 2026-11-15 (best set: `k_cg` 20, `k_cK` 35,
>    `k_ex_g` 0.8, `k_ex_K` 1.2 — see the refit note above). Full
>    results: evaluation repo README, section "Plant alveolar-VA cap
>    fix" and section "Layer 1 refit on fixed plant".
>
> Everything else — model equations, parameters, calibration — is upstream
> as of `f0c3c46`. Evaluation of this fork (Protocol B, paper test, FiCO2
> 8%): https://github.com/zunleashed/metabolisim-exercise-eval, README
> section "OpenCRS fork — raised VAflow clamp".

A closed-loop cardiorespiratory model: a lumped-parameter cardiovascular system
coupled to respiratory mechanics and gas exchange, under baroreflex and
chemoreflex control, simulated at rest and during exercise.

This repository holds the model implementation and the sensitivity, history
matching and Bayesian calibration analyses built on it.

---

## Where the paper's code is

> ### **[`Entire_system/DGSM_Union_Paper_gas/`](Entire_system/DGSM_Union_Paper_gas/)**
>
> All code, configuration and analysis scripts for the accompanying manuscript
> live in this one folder, which has its own detailed
> [README](Entire_system/DGSM_Union_Paper_gas/README.md) covering the
> environment, the pipeline stage by stage, and which script produces which
> figure.

Start there. The rest of this repository is the development history behind it
and is **not** needed to reproduce the paper.

---

## Repository layout

| Path | What it is |
|---|---|
| `Entire_system/DGSM_Union_Paper_gas/` | **The manuscript's analysis.** Union (rest + exercise) DGSM, 8-wave history matching, NUTS calibration, MAP refinement, figures |
| `Entire_system/` (root files) | The full-system model: cardiovascular, respiratory, gas-exchange and control equations (`All_derivatives_njit.py`, `All_Cardiovascular_system.py`, `All_Respiratory_controller.py`, `All_Gas_exchange.py`), plus earlier rest-only and exercise-only calibration pipelines |
| `Entire_system/DGSM_Union_Paper_Final_50/`, `DGSM_Union_Paper_Longer_Delay/`, `DGSM_Union/`, `DGSM_Pericardium/` | Earlier union-calibration variants, kept for provenance |
| `Entire_system/DGSM_Exercise_Paper/`, `DGSM_Exercise_Paper_Final_20/` | Exercise-only calibration variants |
| `Entire_system/Union_Calibration/`, `Plot_abstract/` | Supporting calibration and figure scripts |
| Root `.py` files | The original single-file model implementation that the `Entire_system` version grew out of |
| `Just Resp/`, `Reduced/` | Respiratory-only and reduced-order model experiments |

The `DGSM_*` folders are deliberate snapshots rather than branches: each one
pins the model and parameter ranges used for that run, so results stay
reproducible. They share filenames but are not interchangeable.

---

## Quick start

```bash
git clone https://github.com/Sheng-Ya/OpenCRS.git
cd OpenCRS/Entire_system/DGSM_Union_Paper_gas
```

Then follow that folder's README. It lists pinned dependencies (Python 3.12,
numpy, scipy, numba, SALib, torch, gpytorch, pyro, autoemulate), the seven
pipeline stages with their commands, and the figure-to-script map.

---

## Data

The repository contains source code only. The numerical data behind the
published results — the DGSM design and model outputs, the history-matching
NROY region, the fitted Gaussian-process emulators and the MCMC posterior
samples — is archived on Zenodo under CC BY 4.0:

> ### **https://doi.org/10.5281/zenodo.22732021**

The deposit is a 2.7 GB minimal dataset that reproduces the main reported
results without re-running the simulator, documented by a `MANIFEST.txt` and
verifiable against `CHECKSUMS.sha256`. The complete set of intermediate arrays
the workflow produces is roughly 19 GB.

See the [Data section](Entire_system/DGSM_Union_Paper_gas/README.md#data) of the
paper folder's README for the file-by-file breakdown and for the ~2.8 GB subset
that reproduces the main results.

---

## License

MIT — see [`LICENSE`](LICENSE). You may use, modify and redistribute this code,
including commercially, provided the copyright notice is retained. It is
supplied without warranty.

This is a research model. It is not validated for clinical use and must not be
used to inform patient care.

## Citation

If you use this code, please cite the accompanying manuscript.

## Contact

Please open an issue for questions about reproduction.
