"""Symbolic checks for the algebraic Itô cores in the FA audit.

This script verifies scalar coefficient identities only.  Matrix inner
products such as |Dv|^2 are represented by independent scalar symbols.
"""

import sympy as sp


def verify_normalized_alignment_drift() -> None:
    z, r = sp.symbols("z r", nonzero=True)
    norm_dv2, norm_dc2, dv_dc = sp.symbols(
        "norm_dv2 norm_dc2 dv_dc", real=True
    )
    norm_d2, c_ac = sp.symbols("norm_d2 c_ac", real=True)

    phi = z**2 / r
    drift_r = norm_d2 - 2 * c_ac
    bracket_z = norm_dv2
    bracket_r = 4 * norm_dc2
    bracket_zr = 2 * dv_dc

    ito_drift = (
        sp.diff(phi, r) * drift_r
        + sp.Rational(1, 2) * sp.diff(phi, z, 2) * bracket_z
        + sp.diff(phi, z, r) * bracket_zr
        + sp.Rational(1, 2) * sp.diff(phi, r, 2) * bracket_r
    )
    claimed = (
        norm_dv2 / r
        - 4 * z * dv_dc / r**2
        + 4 * z**2 * norm_dc2 / r**3
        - z**2 * norm_d2 / r**2
        + 2 * z**2 * c_ac / r**2
    )
    assert sp.simplify(ito_drift - claimed) == 0


def verify_cone_radial_angular_drift() -> None:
    q, r = sp.symbols("q r", nonzero=True)
    z_norm, d = sp.symbols("Z d", real=True)
    eta, p, h, k = sp.symbols("eta p h k", real=True)

    projection = q**2 / r
    drift_q = 2 * h * q + eta * r * p + k * r * q
    drift_r = d + 2 * h * r
    bracket_q = z_norm
    bracket_r = 4 * r
    bracket_qr = 2 * q

    ito_projection = (
        sp.diff(projection, q) * drift_q
        + sp.diff(projection, r) * drift_r
        + sp.Rational(1, 2)
        * sp.diff(projection, q, 2)
        * bracket_q
        + sp.diff(projection, q, r) * bracket_qr
        + sp.Rational(1, 2)
        * sp.diff(projection, r, 2)
        * bracket_r
    )

    p_visible = projection
    h_angular = z_norm - p_visible
    claimed_projection = (
        2 * eta * p * q
        + 2 * h * p_visible
        + 2 * k * q**2
        + (h_angular - (d - 1) * p_visible) / r
    )
    assert sp.simplify(ito_projection - claimed_projection) == 0

    drift_z_norm = 2 * eta * p * q + 2 * h * z_norm + 2 * k * q**2
    claimed_angular = (
        2 * h * h_angular
        - (h_angular - (d - 1) * p_visible) / r
    )
    assert sp.simplify(drift_z_norm - ito_projection - claimed_angular) == 0


if __name__ == "__main__":
    verify_normalized_alignment_drift()
    verify_cone_radial_angular_drift()
    print("All FA Itô algebra checks passed.")

