# Thin shell about every center

Let `mu` be isotropic and log-concave on `R^n`, and let `X` have law `mu`.
The dimension-free thin-shell theorem, in its equivalent quadratic-radius
form, states that

\[
 \operatorname {Var}(|X|^2)\le C_0 n.                    \tag{1}
\]

The equivalence between (1) and the usual bounded thin-shell parameter uses
the standard log-concave radial moment estimates; its constants are
universal.

For every deterministic `z in R^n`, put

\[
 Z=|X-z|^2,
 \qquad m=\mathbb EZ=n+|z|^2.                             \tag{2}
\]

Isotropy and (1) give

\[
 \begin{split}
 \operatorname {Var}Z
 &=\operatorname {Var}(|X|^2-2\langle z,X\rangle)\\
 &\le2\operatorname {Var}(|X|^2)
       +8\operatorname {Var}\langle z,X\rangle\\
 &\le 2C_0n+8|z|^2
 \le C_1(n+|z|^2).                                      \tag{3}
 \end{split}
\]

For every nonnegative random variable `Z` with finite second moment and
`m=EZ>0`,

\[
 \operatorname {Var}(\sqrt Z)
 \le\mathbb E(\sqrt Z-\sqrt m)^2
 =\mathbb E\left[{(Z-m)^2\over(\sqrt Z+\sqrt m)^2}\right]
 \le{\operatorname {Var}Z\over m}.                       \tag{4}
\]

Combining (2)--(4) proves the uniform translated thin-shell estimate

\[
 \boxed{\quad
 \sup_{z\in\mathbb R^n}\operatorname {Var}|X-z|\le C_1.
 \quad}                                                   \tag{5}
\]

Consequently, if `g:[0,infinity)->R` is 1-Lipschitz, then for an independent
copy `X'`,

\[
 \operatorname {Var}(g(|X-z|))
 ={1\over2}\mathbb E
   [g(|X-z|)-g(|X'-z|)]^2
 \le\operatorname {Var}|X-z|le C_1.                    \tag{6}
\]

The same conclusion for first absolute moments follows by Cauchy--Schwarz.
Thus the thin-shell obstruction excludes not only radial witnesses centered
at the barycenter, but every translated radial/concurrent-ray witness, with
one universal constant.
