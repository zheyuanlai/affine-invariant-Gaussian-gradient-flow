"""Guard: active README / report prose must not present obsolete theory.

The obsolete formulas of the pre-revision KL proof -- the cubic KL stepsize
penalty ``max{1, lambda_max^3/(2 lambda_min^3)}`` and the clipped
relative-smoothness constant ``beta*max(lambda_plus, lambda_plus^4/lambda_minus^3)``
-- must not appear in the current README or affected reports unless the same
paragraph explicitly flags them as deprecated / obsolete / historical.

Basename is globally unique across tests/ (flat default import mode).
"""
import os
import re

import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Files whose *current* prose is governed by the new theory.
TEXT_FILES = [
    os.path.join(_ROOT, "README.md"),
    os.path.join(_ROOT, "reports",
                 "natural_gradient_discretization_stepsize_report.tex"),
    os.path.join(_ROOT, "reports",
                 "natural_gradient_nonconvex_instability_report.tex"),
]

# Obsolete formula fragments, matched in either ASCII (README) or TeX (reports).
OBSOLETE_PATTERNS = [
    re.compile(r"lambda_plus\s*\^\s*\{?\s*4", re.IGNORECASE),     # lambda_plus^4
    re.compile(r"lambda_minus\s*\^\s*\{?\s*3", re.IGNORECASE),    # lambda_minus^3
    re.compile(r"lambda_max\s*\^\s*\{?\s*3", re.IGNORECASE),      # lambda_max^3
    re.compile(r"lambda_min\s*\^\s*\{?\s*3", re.IGNORECASE),      # lambda_min^3
    re.compile(r"\\lambda_\+\s*\^\s*\{?\s*4"),                    # \lambda_+^4
    re.compile(r"\\lambda_-\s*\^\s*\{?\s*3"),                     # \lambda_-^3
    re.compile(r"\\lambda_\{\\max\}\s*\^\s*\{?\s*3"),             # \lambda_{\max}^3
    re.compile(r"\\lambda_\{\\min\}\s*\^\s*\{?\s*3"),             # \lambda_{\min}^3
]

ALLOW = re.compile(r"deprecat|obsolet|historical", re.IGNORECASE)


def _paragraphs(text):
    """Split into blank-line separated paragraphs (keeps each self-contained)."""
    return re.split(r"\n\s*\n", text)


@pytest.mark.parametrize("path", TEXT_FILES)
def test_no_active_obsolete_theory_text(path):
    assert os.path.exists(path), path
    with open(path, "r", encoding="utf-8") as f:
        text = f.read()
    offenders = []
    for para in _paragraphs(text):
        hit = next((p.pattern for p in OBSOLETE_PATTERNS if p.search(para)), None)
        if hit and not ALLOW.search(para):
            offenders.append((hit, para.strip()[:200]))
    assert not offenders, (
        f"{os.path.basename(path)} has active obsolete-theory text "
        f"(not flagged deprecated/obsolete/historical):\n"
        + "\n".join(f"  [{h}] {snip}" for h, snip in offenders))
