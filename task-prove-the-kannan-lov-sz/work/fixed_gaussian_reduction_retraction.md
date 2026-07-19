# Retraction: the claimed fixed-Gaussian reverse reduction

## Mandatory correction

Let \(a\) be the minimum input spectral gap and \(b\) the spectral gap of
the normalized convolution.  The verified forward inequality is

\[
 b(1-b)\le 4(b-a).
\]

Its correct rearrangement is

\[
 b-a\ge \frac{b(1-b)}4,
 \qquad
 \boxed{a\le \frac{b(3+b)}4}.
\]

It is **not** true that this argument gives
\(a\ge b(3+b)/4\).  Therefore a universal lower bound on the gap after
fixed unit-Gaussian regularization does not transfer back to the original
law by the two-law convolution theorem.

The following claims in the draft state are retracted wherever they
occur:

1. that the two-law inequality gives
   \(\lambda_1(\mu)\ge c_0(3+c_0)/4\) from
   \(\lambda_1(\mathcal G\mu)\ge c_0\);
2. that KLS is thereby quantitatively equivalent to KLS on the fixed
   Gaussian-regularized subclass;
3. that it is enough, on the basis of this inequality, to prove a gap for
   isotropic potentials satisfying \(0\preceq D^2V\preceq2I\).

The posterior calculation

\[
 0\preceq D^2V_{\mathcal G\mu}\preceq2I
\]

remains correct as a structural fact about the Gaussian-smoothed output.
Likewise, the unequal-law forward amplification theorem remains correct.
Neither supplies the missing reverse-smoothing estimate.

This retraction supersedes the fixed-Gaussian reduction passages in
`two_law_convolution_regularization.md` and in any generated research note
until those passages are physically replaced and the PDF is rebuilt.
