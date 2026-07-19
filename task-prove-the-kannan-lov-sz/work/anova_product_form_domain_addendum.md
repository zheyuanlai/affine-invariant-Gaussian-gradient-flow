# Form-domain addendum for the product ANOVA proof

In the product proof, the affine reduction is made before the jet estimates:
\(f\) is assumed orthogonal to constants and all coordinates. For a general
\(f\in L^2_0(\nu)\), choose bounded smooth \(f_j\to f\) in \(L^2(\nu)\), and
then replace each \(f_j\) by
\[
 \widetilde f_j=f_j-\sum_{i=1}^n
       E_\nu[S_i f_j(S)]\,S_i .
\]
(The constant projection is zero because \(Ef_j\) may first be subtracted.)
Each \(\widetilde f_j\) is exactly affine-orthogonal, and
\[
 R(\widetilde f_j)=R(f_j),\qquad
 \operatorname {dist}(\widetilde f_j,\operatorname {Aff})
 =\operatorname {dist}(f_j,\operatorname {Aff}).
\]
The coefficient vectors converge because \(S_i\in L^2(\nu)\), so
\(\widetilde f_j\to f-\Pi_{\operatorname{Aff}}f\). Apply the finite-jet estimate
to \(\widetilde f_j\), then pass by \(L^2\) convergence of the lifts and
conditional expectations. This supplies the orthogonality needed in the
zeroth-jet estimate \(V_0+3M_1\le12V_1\) and avoids an unjustified use of
that estimate for arbitrary approximants.

