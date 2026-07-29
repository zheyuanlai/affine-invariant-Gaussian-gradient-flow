"""Scalar natural gradient schemes on the sharp bump train, and the stepsize arms.

The two scalar maps are eq. (bump-R-map) and eq. (bump-KL-map) of Appendix C,

    Riemannian:  m+ = m - dt c b,   c+ = c exp(dt [1 - c A])
    KL/Bregman:  m+ = m - dt c b,   c+ = (1 + dt) c / (1 + dt c A)

written directly in the scalar form the appendix uses rather than through the
dense ``discretization_step`` of the stepsize group (identical maps at ``d = 1``,
but without a per-step eigendecomposition, which matters at ``Omega(kappa^2)``
iterations).

Three stepsize arms, all with the same ``gamma``:

``theory``   ``dt = gamma / kappa``      -- the globally certified fixed step of
             ``thm:glob-Riem``/``thm:glob-KL``, the one ``thm:sharp-disc`` targets.
``const``    ``dt = gamma``              -- a ``kappa``-independent fixed step, outside
             the certified range; the regime the manuscript's remark after
             ``thm:sharp-disc`` leaves open.
``relcurv``  ``dt = gamma / (c_n A_n)``  -- the adaptive relative-curvature rule that
             same remark names ("a method whose stability test is based on the
             current relative curvature could accept a step of order one").

Runs stop at a relative gap tolerance, on non-finite state, or at ``max_steps``.
"""
from __future__ import annotations

import math

SCHEMES = ["riemannian", "kl"]
ARMS = ["theory", "const", "relcurv"]

# Bures--Wasserstein forward--backward step (manuscript eq. BW-step; Lambert et al.,
# arXiv:2304.05398), scalar form. Its certified stepsize is eta <= 1/beta, which
# involves beta ALONE -- no covariance factor. That is the structural contrast with
# Fisher--Rao, whose certified step 1/(beta lambda_max) carries the covariance scale.
BW_SCHEMES = ["bures_wasserstein"]


def bw_step(m, c, b, A, eta):
    """One Bures--Wasserstein forward--backward step.

    ``m+ = m - eta b`` -- note the mean update has NO covariance preconditioner, so
    its multiplier is ``1 - eta * (secant slope)`` and the secant is bounded by
    ``beta``. With ``eta <= 1/beta`` the multiplier lies in ``[0, 1)`` at every
    admissible state, which is exactly why the two-cycle cannot form.
    """
    m_next = m - eta * b
    c_tilde = (1.0 - eta * A) ** 2 * c
    c_next = 0.5 * (c_tilde + 2.0 * eta + math.sqrt(c_tilde * (c_tilde + 4.0 * eta)))
    return m_next, c_next


def nominal_dt(arm, gamma, kappa):
    """Nominal (trajectory) stepsize of an arm -- what the retuned train is built for.

    ``relcurv`` has ``dt = gamma/(cA)`` and the construction keeps ``cA = 1 + o(1)``
    along the shadowed trajectory, so its nominal step is ``gamma``, same as ``const``.
    """
    if arm == "theory":
        return gamma / kappa
    if arm in ("const", "relcurv"):
        return gamma
    raise ValueError(f"unknown arm '{arm}' (known: {ARMS})")


def step_size(arm, gamma, kappa, c, A):
    """Actual stepsize used at the current state.

    ``relcurv`` is ``gamma / max(c A, 1)``: a *stability test*, so it only ever
    shrinks the order-one step, never inflates it. Without the cap the rule
    inflates ``dt`` to ``gamma/(cA)`` wherever the relative curvature is small
    (``cA = 1/kappa`` once the mean leaves the bump train), which makes the
    Riemannian covariance map ``c exp(dt[1 - cA])`` overshoot by ``e^{gamma kappa}``
    -- an artifact of the uncapped rule, not of the construction.
    """
    if arm == "relcurv":
        return gamma / max(c * A, 1.0)
    return nominal_dt(arm, gamma, kappa)


def scheme_step(scheme, m, c, b, A, dt):
    """One step of the named scalar scheme."""
    m_next = m - dt * c * b
    if scheme == "riemannian":
        c_next = c * math.exp(dt * (1.0 - c * A))
    elif scheme == "kl":
        c_next = (1.0 + dt) * c / (1.0 + dt * c * A)
    else:
        raise ValueError(f"unknown scheme '{scheme}' (known: {SCHEMES})")
    return m_next, c_next


def simulate_two_cycle(target, scheme, dt, max_steps, tol_rel=1e-6,
                       max_saved_rows=400, perturb=0.0):
    """Run one scheme on the two-cycle target from ``(M, 1/p)``.

    ``scheme`` is ``riemannian``, ``kl`` or ``bures_wasserstein``; ``dt`` is the fixed
    stepsize (``gamma`` for the Fisher--Rao schemes, ``eta`` for BW). ``perturb``
    multiplies the initial mean by ``1 + perturb``, which probes whether the orbit is
    attracting or a knife edge.

    Returns ``(records, summary)``. The summary reports the cycle diagnostics: the
    sign pattern of the mean, the amplitude drift ``|m_n|/M``, and whether the gap
    ever falls below ``tol_rel`` of its initial value.
    """
    m = target.M * (1.0 + perturb)
    c = target.c
    gap0 = target.gap(m, c)
    save_every = max(1, max_steps // max_saved_rows)

    records, n_tol = [], -1
    amp_min, amp_max = abs(m) / target.M, abs(m) / target.M
    gap_min, gap_max = gap0, gap0
    sign_flips = 0
    prev_sign = math.copysign(1.0, m)
    status = "ok"
    gap = gap0

    for n in range(max_steps):
        b, A = target.b_A(m, c)
        if n % save_every == 0:
            records.append({"n": n, "m": m, "c": c, "gap": gap,
                            "rel_gap": gap / gap0, "amplitude": abs(m) / target.M,
                            "cA": c * A})
        if scheme == "bures_wasserstein":
            m, c = bw_step(m, c, b, A, dt)
        else:
            m, c = scheme_step(scheme, m, c, b, A, dt)
        if not (math.isfinite(m) and math.isfinite(c) and c > 0.0):
            status = "failed"
            break
        gap = target.gap(m, c)
        if not math.isfinite(gap):
            status = "failed"
            break
        amp = abs(m) / target.M
        amp_min, amp_max = min(amp_min, amp), max(amp_max, amp)
        gap_min, gap_max = min(gap_min, gap), max(gap_max, gap)
        s = math.copysign(1.0, m) if m != 0.0 else prev_sign
        if s != prev_sign:
            sign_flips += 1
        prev_sign = s
        if n_tol < 0 and gap <= tol_rel * gap0:
            n_tol = n + 1
            break

    n_done = len(records) and records[-1]["n"]
    records.append({"n": max_steps if n_tol < 0 else n_tol, "m": m, "c": c, "gap": gap,
                    "rel_gap": gap / gap0, "amplitude": abs(m) / target.M, "cA": float("nan")})
    converged = n_tol > 0
    if status == "ok" and not converged:
        status = "cycling" if gap_min > 0.5 * gap0 else "max_steps"
    return records, {
        "kappa": target.kappa, "scheme": scheme, "gamma": target.gamma, "dt": dt,
        "perturb": perturb,
        "p": target.p, "r": target.r, "q": target.q, "c_cycle": target.c, "M": target.M,
        "initial_gap": gap0, "final_gap": gap, "final_rel_gap": gap / gap0,
        "min_rel_gap": gap_min / gap0, "max_rel_gap": gap_max / gap0,
        "amplitude_min": amp_min, "amplitude_max": amp_max,
        "n_steps": max_steps if n_tol < 0 else n_tol,
        "sign_flips": sign_flips, "n_tol": n_tol, "tol_rel": tol_rel,
        "converged": int(converged), "status": status,
    }


def simulate(train, scheme, arm, gamma, max_steps, tol_rel=1e-6,
             max_saved_rows=400):
    """Run one ``(train, scheme, arm)`` from ``(m, c) = (x_0, 1/kappa)``.

    Returns ``(records, summary)``. ``records`` is the decimated trajectory;
    ``summary`` reports the iteration counts to a constant-factor and to a
    ``tol_rel`` relative reduction of the gap, the shadowing diagnostics against
    the ideal centers ``x_j`` (meaningful while ``n <= train.N``), and stability
    flags. Iteration counts are computed from the full trajectory.
    """
    kappa = train.kappa
    m, c = train.x0, train.c_kappa
    gap0 = train.gap(m, c)
    if not (gap0 > 0.0 and math.isfinite(gap0)):
        raise ValueError(f"non-positive initial gap {gap0!r}")

    save_every = max(1, max_steps // max_saved_rows)
    records = []
    n_half = n_tol = -1
    max_energy_increase = 0.0
    shadow_m = shadow_q = 0.0       # sup |m_n - x_n| / w  and  sup |kappa c_n - 1|
    cA_min, cA_max = math.inf, -math.inf
    c_min, c_max = c, c
    dt_sum = 0.0
    status = "ok"
    prev_gap = gap0
    n = 0

    while n < max_steps:
        b, A = train.b_A(m, c)
        cA = c * A
        cA_min, cA_max = min(cA_min, cA), max(cA_max, cA)
        if n <= train.N:
            shadow_m = max(shadow_m, abs(m - train.centers[n]) / train.w)
            shadow_q = max(shadow_q, abs(kappa * c - 1.0))
        if n % save_every == 0:
            records.append({"n": n, "m": m, "c": c, "gap": prev_gap,
                            "rel_gap": prev_gap / gap0, "cA": cA, "A": A})
        dt = step_size(arm, gamma, kappa, c, A)
        dt_sum += dt
        m, c = scheme_step(scheme, m, c, b, A, dt)
        n += 1
        if not (math.isfinite(m) and math.isfinite(c) and c > 0.0):
            status = "failed"
            break
        c_min, c_max = min(c_min, c), max(c_max, c)
        gap = train.gap(m, c)
        if not math.isfinite(gap):
            status = "failed"
            break
        max_energy_increase = max(max_energy_increase, gap - prev_gap)
        prev_gap = gap
        if n_half < 0 and gap <= 0.5 * gap0:
            n_half = n
        if n_tol < 0 and gap <= tol_rel * gap0:
            n_tol = n
            break

    records.append({"n": n, "m": m, "c": c, "gap": prev_gap,
                    "rel_gap": prev_gap / gap0, "cA": cA_max, "A": float("nan")})
    if status == "ok" and n_tol < 0:
        status = "max_steps"
    return records, {
        "kappa": kappa, "scheme": scheme, "arm": arm, "gamma": gamma,
        "s_train": train.s, "N_train": train.N, "x0": train.x0,
        "dt_nominal": nominal_dt(arm, gamma, kappa),
        "dt_mean": dt_sum / max(n, 1),
        "initial_gap": gap0, "final_gap": prev_gap,
        "final_rel_gap": prev_gap / gap0,
        "n_steps": n,
        "n_half": n_half, "n_tol": n_tol, "tol_rel": tol_rel,
        "n_half_over_N_train": (n_half / train.N) if (n_half > 0 and train.N > 0) else float("nan"),
        "shadow_mean_err_over_w": shadow_m, "shadow_cov_err": shadow_q,
        "cA_min": cA_min, "cA_max": cA_max,
        "c_min": c_min, "c_max": c_max,
        "max_energy_increase": max_energy_increase,
        "monotone": int(max_energy_increase <= 1e-9 * gap0),
        "status": status,
    }
