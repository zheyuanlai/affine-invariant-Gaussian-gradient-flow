# Nonuniqueness of Stein kernels: a sharp trace-budget sanity check

Even for the standard Gaussian, a Hilbert--Schmidt Stein budget does not control
the operator error of an arbitrary Stein kernel.

Let \(\gamma_2(x)=(2\pi)^{-1}e^{-|x|^2/2}\), and choose a nonzero
\(\psi\in C_c^\infty(\mathbb R^2)\). Define the symmetric matrix field
\[
 B_\psi(x)=\gamma_2(x)^{-1}
 \begin{pmatrix}
  \partial_{22}\psi(x)&-\partial_{12}\psi(x)\\
  -\partial_{12}\psi(x)&\partial_{11}\psi(x)
 \end{pmatrix}. \tag{1}
\]
Its weighted row divergences vanish:
\[
 \partial_1(\gamma_2 B_{11})+\partial_2(\gamma_2 B_{12})
 =\partial_{122}\psi-\partial_{212}\psi=0,
\]
\[
 \partial_1(\gamma_2 B_{21})+\partial_2(\gamma_2 B_{22})
 =-\partial_{112}\psi+\partial_{211}\psi=0. \tag{2}
\]
Therefore, for every smooth compactly supported \(u\),
\[
 E_{\gamma_2}[B_\psi\nabla u]=0. \tag{3}
\]
Hence \(\tau_L=I+L B_\psi\) is a symmetric Stein kernel of \(\gamma_2\) for
every scalar \(L\). Because \(\psi\) is compactly supported, all polynomial
moments of \(B_\psi\) are finite, but
\[
 \|E[(\tau_L-I)^2]\|_{\rm op}=L^2
 \|E[B_\psi^2]\|_{\rm op}\longrightarrow\infty. \tag{4}
\]
Tensoring with \(I_{n-2}\) gives the same example in every dimension.

Thus the usual existence of a Stein kernel, even one with finite (or
trace-controlled) Hilbert--Schmidt discrepancy, says nothing about the
operator budget needed in the ordered-tensor ANOVA proof. The canonical
Gaussian kernel \(I\) is of course optimal; the example only rules out
arguments that do not specify and control the chosen kernel.

