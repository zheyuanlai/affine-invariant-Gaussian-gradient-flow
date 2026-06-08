"""Display-only plotting helpers and style constants for the WFR experiment.

Matplotlib only. None of these mutate or re-save experiment data; they only
control how the committed CSV values are displayed in the report figures.
"""
from __future__ import annotations

import numpy as np

# Method display order, colors, labels, and line styles.
METHOD_ORDER = ["w_only", "fr_only", "wfr_fixed", "wfr_theory", "wfr_adaptive"]
METHOD_STYLE = {
    "w_only":      {"color": "#9467bd", "label": "W-only",       "ls": ":"},
    "fr_only":     {"color": "#d62728", "label": "FR-only",      "ls": "--"},
    "wfr_fixed":   {"color": "#1f77b4", "label": "WFR-fixed",    "ls": "-"},
    "wfr_theory":  {"color": "#7f7f7f", "label": "WFR-theory",   "ls": "-."},
    "wfr_adaptive": {"color": "#2ca02c", "label": "WFR-adaptive", "ls": "-"},
}
# Methods shown in the four-curve phase-separation figures (Figs 1-2).
PHASE_METHODS = ["w_only", "fr_only", "wfr_fixed", "wfr_adaptive"]

TARGET_TITLE = {
    "gaussian": "Gaussian posterior",
    "smooth_log_cosh": "smooth log-cosh posterior",
}

# Display floor for the objective-gap log axes (data is never modified).
GAP_FLOOR = 1e-13


def clip_gap(gap, floor=GAP_FLOOR):
    """Clip a gap array to a positive display floor for log-scale plotting."""
    g = np.asarray(gap, dtype=np.float64)
    return np.clip(g, floor, None)
