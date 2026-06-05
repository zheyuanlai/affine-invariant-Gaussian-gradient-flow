\section{Local convergence near equilibrium}
We proceed to prove the local convergence rate of Gaussian natural gradient flow \eqref{ODEs:NG} near equilibrium under Assumption \ref{assump:logconcave-smooth}, where we denote $\kappa = \frac{\beta}{\alpha}$. Based on the evidence from numerical experiments, we expect that the local convergence rate Gaussian natural gradient flow does not rely on initialization as in Theorem \ref{thm:glob-conv}.

We work on the equilibrium-whitened coordinates without loss of generality
\begin{equation}\label{eq:whitened}
    z = C_\star^{-1/2} (\theta - m_\star)
\end{equation}
such that the optimality is $a_\star = \mathcal{N} (0, I_{N_\theta})$ and $z \sim \mathcal{N} (0, I_{N_\theta})$, in which case we have $\alpha \leq 1 \leq \beta$.
This gives the stationary conditions $\E_{\mathcal{N} (0, I_{N_\theta})} [\nabla V(z)] = 0$ and $\E_{\mathcal{N} (0, I_{N_\theta})} [\nabla^2 V(z)] = I_{N_\theta}$. Hence for any $a = (m, C) \in \Aspace$, we can parameterize $(m, C) = (u, I + X)$ with $u \in \R^{N_\theta}$ and $X \in \Sym (N_\theta)$.

\begin{definition}\label{def:local-operators}
    For $u \in \R^{N_\theta}$ and $X \in \Sym (N_\theta)$, define the operators,
    \begin{align*}
        T [X]_k &:= \E[\Tr (X \nabla^2 V(z)) z_k] = \sum_{i, j} \E[\partial^3_{ijk} V(z)] X_{ij}, \\
        T^* [u]_{ij} &:= \E[(u^\top z) \nabla^2 V_{ij} (z)], \\
        H [X]_{ij} &:= \E[\nabla^2 V_{ij} (z) (z^\top Xz - \Tr X)] = \sum_{k, l} \E[\partial^4_{ijkl} V(z)] X_{kl},
    \end{align*}
    where we have applied Stein's identity in $\mathcal{T}$ and $\mathcal{H}$. For the diagonal mode, we define
    \begin{equation*}
        \widetilde{T}_{ki} := \E[\nabla^2 V_{ii} (z) z_k], \qquad G_{ij} := \E[\nabla^2 V_{ii} (z) z_j^2].
    \end{equation*}
\end{definition}

We first show the tail-cut bound of $G_{ij}$.

\begin{lemma}[Tail-cut bound]\label{lem:tailcut}
    Under Assumption \ref{assump:logconcave-smooth} and \eqref{eq:whitened}, for any $i, j$, we have
    \begin{equation*}
        \alpha \leq G_{ij} \leq 4 + \frac{4}{\sqrt{\pi}} ( 1 + \log \kappa).
    \end{equation*}
\end{lemma}

\begin{proof}
    Since $\alpha I \preceq \nabla^2 V \preceq \beta I$, then $\alpha \leq \nabla^2 V_{ii} (z) \leq \beta$, hence $G_{ij} = \E[\nabla^2 V_{ii} (z) z_j^2] \geq \alpha \E[z_j^2] = \alpha$ since $z_j \sim \N (0, 1)$. This proves the lower bound. For the upper bound, if $\kappa \leq e^{1/2}$, then $G_{ij} \leq \beta \leq \kappa \leq e^{1/2} < 4 + \frac{4}{\sqrt{\pi}} ( 1 + \log \kappa)$, so we consider the case where $\kappa > e^{1/2}$ and choose $t_0 = \sqrt{2 \log \kappa}> 1$. We split the expectation into the bulk region and the tail region:
    \begin{equation*}
        G_{ij} = \E[\nabla^2 V_{ii} (z) z_j^2] = \E[\nabla^2 V_{ii} (z) z_j^2 \mathbf{1}_{|z_j| < t_0}] + \E[\nabla^2 V_{ii} (z) z_j^2 \mathbf{1}_{|z_j| \geq t_0}].
    \end{equation*}
    For the bulk part,
    \begin{equation*}
        \E[\nabla^2 V_{ii} (z) z_j^2 \mathbf{1}_{|z_j| < t_0}] \leq t_0^2 \E[\nabla^2 V_{ii} (z)] = 2 \log \kappa.
    \end{equation*}
    For the tail part, we define $\varphi (t) = \frac{1}{\sqrt{2\pi}} e^{-t^2 / 2}$, then
    \begin{equation*}
        \E[\nabla^2 V_{ii} (z) z_j^2 \mathbf{1}_{|z_j| \geq t_0}] \leq \kappa \E[z_j^2 \mathbf{1}_{|z_j| \geq t_0}] = 2 \kappa \int_{t_0}^\infty t^2 \varphi(t) \dd t,
    \end{equation*}
    using integration by parts and Mills ratio,
    \begin{equation*}
        \int_{t_0}^\infty t^2 \varphi(t) \dd t = t_0 \varphi (t_0) + \int_{t_0}^\infty \varphi(t) \dd t \leq t_0 \varphi (t_0) + \frac{\varphi (t_0)}{t_0} \leq 2t_0 \varphi(t_0) = 2\sqrt{\frac{\log \kappa}{\pi}} \kappa^{-1}.
    \end{equation*}
    Hence,
    \begin{equation*}
        G_{ij} \leq 2 \log \kappa + \frac{4}{\sqrt{\pi}} \sqrt{\log \kappa} \leq 4 + \frac{4}{\sqrt{\pi}} (1 + \log \kappa).
    \end{equation*}
\end{proof}

We then prove an identity which characterizes the ``perturbation'' near equilibrium.

\begin{lemma}\label{lem:bochner}
    Under Assumption \ref{assump:logconcave-smooth} and \eqref{eq:whitened}, for any $u \in \R^{N_\theta}$, $X \in \Sym (N_\theta)$, and $\eta(z) := u + Xz$, 
    \begin{equation*}
        \E[\eta(z)^\top \nabla^2 V(z) \eta (z)] = \|u\|^2 + 2 u^\top T[X] + \|X\|_F^2 + \Tr(X H[X]).
    \end{equation*}
    Specifically when $X = \diag \lambda$,
    \begin{equation*}
        \E[\eta(z)^\top \nabla^2 V(z) \eta (z)] = \|u\|^2 + 2 u^\top \widetilde{T} \lambda + \|\lambda\|^2 + \lambda^\top (G - \mathbf{1} \mathbf{1}^\top) \lambda.
    \end{equation*}
\end{lemma}

\begin{proof}
    We expand 
    \begin{equation*}
        \eta^\top \nabla^2 V \eta = u^\top \nabla^2 V u + 2u^\top \nabla^2 V (Xz) + (Xz)^\top \nabla^2 V (Xz)
    \end{equation*}
    and analyze the expectation for each term.
    For the constant term, since $\E[\nabla^2 V] = I$, then $\E[u^\top \nabla^2 V u] = u^\top \E[\nabla^2 V] u = \|u\|^2$. For the cross term, by Stein's Lemma and the definition of $T$, 
    \begin{equation*}
        \E[u^\top \nabla^2 V (Xz)] = \sum_{a, b, c} u_a X_{bc} \E[\nabla^2 V_{ab} (z) z_c] = \sum_{a, b, c} u_a X_{bc} \E[\partial_a \partial_b \partial_c V (z)] = u^\top T[X].
    \end{equation*}
    For the quadratic term, we apply Stein's Lemma twice to obtain 
    \begin{align*}
        \E[(Xz)^\top \nabla^2 V (Xz)] &= \sum_{a, b, c, d} X_{ac} X_{bd} \E[\nabla^2 V_{ab} (z) z_c z_d] \\
        &= \sum_{a, b, c, d} X_{ac} X_{bd} (\E[\partial_a \partial_b \partial_c \partial_d V (z)] + \delta_{cd} \delta_{ab}) = \Tr(X H[X]) + \|X\|_F^2.
    \end{align*}
    Adding these three terms up completes the proof.
\end{proof}

Based on Lemma \ref{lem:bochner}, we aim to prove the spectral upper bound of $\widetilde{T}^\top \widetilde{T}$.

\begin{lemma}\label{lem:matrix-schur}
    Under Assumption \ref{assump:logconcave-smooth} and \eqref{eq:whitened}, 
    \begin{equation*}
        \widetilde{T}^\top \widetilde{T} \preceq I + (G - \mathbf{1} \mathbf{1}^\top).
    \end{equation*}
\end{lemma}

\begin{proof}
    We choose $\eta (z) = u + \diag (\lambda) z$ in Lemma \ref{lem:bochner}, since $\nabla^2 V(z) \succeq 0$, then
    \begin{equation*}
        0 \leq  \|u\|^2 + 2 u^\top \widetilde{T} \lambda + \|\lambda\|^2 + \lambda^\top (G - \mathbf{1} \mathbf{1}^\top) \lambda
    \end{equation*}
    for any $u \in \R^{N_\theta}$ and $\lambda \in \R^{N_\theta}$, which is equivalent to 
    \begin{equation*}
        \begin{pmatrix}
            I & \widetilde{T} \\
            \widetilde{T} & I + (G - \mathbf{1} \mathbf{1}^\top)
        \end{pmatrix} \succeq 0.
    \end{equation*}
    As $I \succeq 0$, then we have $I + (G - \mathbf{1} \mathbf{1}^\top) - \widetilde{T}^\top \widetilde{T} \succeq 0$ by Schur complement criterion, which completes the proof.
\end{proof}

Based on these Lemmas, we prove the following operator bound.

\begin{lemma}[Operator bound]\label{lem:operator-bound}
    Under Assumption \ref{assump:logconcave-smooth} and \eqref{eq:whitened}, 
    \begin{equation*}
        \lambda_{\max} (G - \mathbf{1} \mathbf{1}^\top) \leq \Gamma := \min \left\{\beta - 1, N_\theta \left(3 + \frac{4}{\sqrt{\pi}} (1 + \log \kappa) \right) \right\}.
    \end{equation*}
\end{lemma}

\begin{proof}
    We prove two spectral upper bounds of $G - \mathbf{1} \mathbf{1}^\top$ separately. On the one hand, we take $\eta(z) = \diag (\lambda) z$ as in Lemma \ref{lem:bochner}, this gives
    \begin{equation*}
        \E[\eta(z)^\top \nabla^2 V (z) \eta (z)] = \|\lambda\|^2 + \lambda^\top (G - \mathbf{1} \mathbf{1}^\top) \lambda
    \end{equation*}
    and concurrently
    \begin{equation*}
        \E[\eta(z)^\top \nabla^2 V (z) \eta (z)] \leq \beta \E \|\diag (\lambda) z\|^2 = \beta \|\diag \lambda\|_F^2 = \beta \|\lambda\|^2.
    \end{equation*}
    Collectively, $\lambda^\top (G - \mathbf{1} \mathbf{1}^\top) \lambda \leq (\beta - 1) \|\lambda\|^2$, hence $\lambda_{\max} (G - \mathbf{1} \mathbf{1}^\top) \leq \beta - 1$.

    On the other hand, by Lemma \ref{lem:tailcut}, we obtain
    \begin{equation*}
        |G_{ij} - 1| \leq 3 + \frac{4}{\sqrt{\pi}} (1 + \log \kappa),
    \end{equation*}
    this gives
    \begin{align*}
        \lambda^\top (G - \mathbf{1} \mathbf{1}^\top) \lambda &= \sum_{i,j} \lambda_i \lambda_j (G_{ij} - 1) \leq \sum_{i, j} |\lambda_i| \cdot |\lambda_j| \left( 3 + \frac{4}{\sqrt{\pi}} (1 + \log \kappa) \right) \\
        &\leq \left( \sum_{i} |\lambda_i| \right)^2 \left( 3 + \frac{4}{\sqrt{\pi}} (1 + \log \kappa) \right) \leq  \left( 3 + \frac{4}{\sqrt{\pi}} (1 + \log \kappa) \right) N_\theta \|\lambda\|^2,
    \end{align*}
    this gives $\lambda_{\max} (G - \mathbf{1} \mathbf{1}^\top ) \leq \left( 3 + \frac{4}{\sqrt{\pi}} (1 + \log \kappa) \right) N_\theta$. Hence collectively,
    \begin{equation*}
        \lambda_{\max} (G - \mathbf{1} \mathbf{1}^\top) \leq \min \left\{\beta - 1, N_\theta \left(3 + \frac{4}{\sqrt{\pi}} (1 + \log \kappa) \right) \right\}.
    \end{equation*}
\end{proof}

We then compute the linearized Jacobian matrix of Gaussian natural gradient flow \eqref{ODEs:NG} and bound its largest eigenvalue.

\begin{proposition}\label{prop:Jac-eig}
    We consider the Gaussian natural gradient flow \eqref{ODEs:NG} under Assumption \ref{assump:logconcave-smooth} and \eqref{eq:whitened}. For $u = m \in \R^{N_\theta}$ and $X = C - I \in \Sym (N_\theta)$, the linearized Jacobian matrix at equilibrium is given by
    \begin{equation*}
        - J_\star (u, X) = \left( u + \frac{1}{2} T[X], X + T^*[u] + \frac{1}{2} H[X] \right).
    \end{equation*}
    Moreover, the largest eigenvalue $\lambda_{\star, \max}$ of $J_\star$ satisfies
    \begin{equation*}
        -\lambda_{\star, \max} \geq \frac{1}{4 + \Gamma} = \frac{1}{4 + \min \left\{\beta - 1, N_\theta \left(3 + \frac{4}{\sqrt{\pi}} (1 + \log \kappa) \right) \right\}},
    \end{equation*}
    and the smallest eigenvalue $\lambda_{\star, \min}$ of $J_\star$ satisfies
    \begin{equation*}
        -\lambda_{\star, \min} \leq \beta + \frac{1 - \alpha}{2}.
    \end{equation*}
\end{proposition}

\begin{proof}
    We apply Bonnet's identity and Price's identity, i.e.,
    \begin{equation*}
        \partial_m \E_{\N (m, C)} [f] = \E_{\N(m, C)} [\nabla f], \qquad \partial_C \E_{\N(m, C)} [f] = \frac{1}{2} \E_{\N(m, C)} [\nabla^2 f],
    \end{equation*}
    by Taylor's expansion near $a_\star$, we obtain
    \begin{align*}
        \E_{\N (u, I + X)} [\nabla V] &= u + \frac{1}{2} T[X] + O (\| (u, X) \|^2), \\
        \E_{\N (u, I + X)} [\nabla^2 V] &= I + T^*[u] + \frac{1}{2} H[X] + O(\|(u, X)\|^2).
    \end{align*}
    By substituting $\E_{\N (u, I + X)} [\nabla V]$ and $\E_{\N (u, I + X)} [\nabla^2 V]$ into \eqref{ODEs:NG}, we obtain the local ODEs
    \begin{equation*}
        \begin{cases}
            \dot u &= -u - \frac{1}{2} T[X] + O(\|(u, X)\|_\star^2), \\
            \dot X &= -X - T^*[u] - \frac{1}{2} H[X] + O(\|(u, X)\|_\star^2).
        \end{cases}
    \end{equation*}
    Therefore the linearized Jacobian matrix is 
    \begin{equation*}
        - J_\star (u, X) = \left( u + \frac{1}{2} T[X], X + T^*[u] + \frac{1}{2} H[X] \right).
    \end{equation*}

    We proceed to show that the operator $J_\star$ is self-adjoint with respect to $\langle \cdot, \cdot \rangle_\star^\mathrm{FR}$. For $(u_1, X_1), (u_2, X_2) \in \R^{N_\theta} \times \Sym (N_\theta)$, we calculate
    \begin{align*}
        \langle -J_\star (u_1, X_1), (u_2, X_2) \rangle_\star^\mathrm{FR} &= u_1^\top u_2 + \frac{1}{2} T[X_1]^\top u_2 + \frac{1}{2} \Tr(X_1 X_2) + \frac{1}{2} \Tr(T^*[u_1] X_2) + \frac{1}{4} \Tr (H[X_1] X_2),
    \end{align*}
    since $u^\top T[X] = \Tr(T^*[u] X)$ and $\Tr(H[X_1] X_2) = \Tr(X_1 H[X_2])$, then the expression is symmetric, which justifies the self-adjointness.

    Thus by the spectral theorem, we obtain
    \begin{equation*}
        -\lambda_{\star, \max} = \min_{(u, X) \neq 0} \frac{\langle -J_\star (u, X), (u, X) \rangle_\star^\mathrm{FR}}{\langle (u, X), (u, X) \rangle_\star^\mathrm{FR}} = \min_{(u, X) \ne 0} \frac{Q(u, X)}{\|u\|^2 + \frac{1}{2} \|X\|_F^2},
    \end{equation*}
    where $Q(u, X) := \|u\|^2 + u^\top T[X] + \frac{1}{2} \|X\|_F^2 + \frac{1}{4} \Tr(X H[X])$. 
    Since the standard Gaussian measure and the Hessian bounds are invariant under orthogonal changes of basis, we work in orthogonal basis $X = \diag \lambda$, hence
    \begin{equation*}
        Q(u, \diag \lambda) = \|u\|^2 + u^\top \widetilde{T} \lambda + \frac{1}{2} \|\lambda\|^2 + \frac{1}{4} \lambda^\top (G - \mathbf{1} \mathbf{1}^\top) \lambda.
    \end{equation*}
    We hence only need to show $Q(u, \diag \lambda) \geq \frac{1}{4 + \Gamma} (\|u\|^2 + \frac{1}{2} \|\lambda\|^2)$, which is equivalent to 
    \begin{equation*}
        \left( 1 - \nu \right) \|u\|^2 + u^\top \widetilde{T} \lambda + \frac{1 - \nu}{2} \|\lambda\|^2 + \frac{1}{4} \lambda^\top (G - \mathbf{1} \mathbf{1}^\top) \lambda \geq 0,
    \end{equation*}
    where we denote $\nu = \frac{1}{4 + \Gamma}$ and $0 < \nu \leq \frac{1}{4}$.
    It suffices to show that
    \begin{equation*}
        \frac{1 - \nu}{2} I + \frac{1}{4} (G - \mathbf{1} \mathbf{1}^\top) - \frac{\widetilde{T}^\top \widetilde{T}}{4(1 - \nu)} \succeq 0,
    \end{equation*}
    from Lemma \ref{lem:matrix-schur} and Lemma \ref{lem:operator-bound}, a sufficient condition is
    \begin{equation*}
        \frac{1 - \nu}{2} - \frac{1 + \nu \Gamma}{4 (1 - \nu)} \geq 0 \Longleftrightarrow f(\nu) := 1 - \nu (4 + \Gamma) + 2\nu^2 \geq 0.
    \end{equation*}
    The smaller root of $f$ is 
    \begin{equation*}
        \nu_1 = \frac{(4 + \Gamma) - \sqrt{(4 + \Gamma)^2 - 8}}{4} = \frac{2}{4 + \Gamma + \sqrt{(4 + \Gamma)^2 - 8}} \geq \frac{1}{4 + \Gamma},
    \end{equation*}
    hence $f (\frac{1}{4 + \Gamma}) \geq 0$, which completes the proof of the largest eigenvalue bound. 

    We proceed to bound the smallest eigenvalue of $J_\star$. From Lemma \ref{lem:bochner} and Assumption \ref{assump:logconcave-smooth}, we have
    \begin{equation*}
        \E[(u + Xz)^\top \nabla^2 V(z) (u + Xz)] = \|u\|^2 + 2 u^\top T[X] + \|X\|_F^2 + \Tr (XH[X]) \leq \beta (\|u\|^2 + \|X\|_F^2),
    \end{equation*}
    also by taking $u = 0$, we obtain
    \begin{equation*}
        \|X\|_F^2 + \Tr(X H[X]) \geq \alpha \|X\|_F^2 \Longrightarrow \Tr(XH[X]) \geq (\alpha - 1) \|X\|_F^2.
    \end{equation*}
    We then upper bound $Q(u, X)$ by $\|(u, X)\|_\star^2$:
    \begin{align*}
        Q(u, X) &= \|u\|^2 + u^\top T[X] + \frac{1}{2} \|X\|_F^2 + \frac{1}{4} \Tr (XH[X]) \\
        &= \frac{1}{2} \E[(u + Xz)^\top \nabla^2 V(z) (u + Xz)] + \frac{1}{2} \|u\|^2 - \frac{1}{4} \Tr (XH[X]) \\
        &\leq \frac{\beta + 1}{2} \|u\|^2 + \frac{2\beta + 1 - \alpha}{4} \|X\|_F^2 \\
        &\leq \max\left\{ \frac{\beta + 1}{2}, \frac{2\beta + 1 - \alpha}{2} \right\} \|(u, X)\|_\star^2 \leq \left( \beta + \frac{1 - \alpha}{2} \right) \|(u, X)\|_\star^2,
    \end{align*}
    where the last inequality comes from $\alpha \leq 1 \leq \beta$ following \eqref{eq:whitened}. Then by spectral theorem, 
    \begin{equation*}
        -\lambda_{\star, \min} \leq \beta + \frac{1 - \alpha}{2}.
    \end{equation*}
\end{proof}

We provide local convergence rate for the continuous dynamics \eqref{ODEs:NG}. Before that, we first provide explicit second-order bound for the natural gradient vector field.

\begin{lemma}\label{lem:2ord-bound-NG}
    Under Assumption \ref{assump:logconcave-smooth} and \eqref{eq:whitened}, we define the vector field of \eqref{ODEs:NG} as
    \begin{equation*}
        F (m, C) = (F_1 (m, C), F_2 (m, C)) := (-C \E_{\N(m, C)} [\nabla V], C - C \E_{\N (m, C)} [\nabla^2 V] C),
    \end{equation*}
    then on the neighborhood $\mathcal{S}_\star = \left\{ (m, C) : \frac{1}{2} I \preceq C \preceq \frac{3}{2} I \right\}$, the second differential of vector field $F$ satisfies
    \begin{equation*}
        \| D^2 F(a) [\zeta, \zeta] \|_\star \leq 135 \beta N_\theta^{5/2} \|\zeta\|_\star^2,
    \end{equation*}
    where $a = (m, C) \in \mathcal{S}_\star$ and the tangent direction $\zeta = (v, Y) \in \R^{N_\theta} \times \Sym(N_\theta)$.
\end{lemma}

\begin{proof}
    We consider the line $a(s) = (m + sv, C + sY)$, since $\|\zeta\|_\star^2 = \|v\|^2 + \frac{1}{2} \|Y\|_F^2$, then we have $\|v\| \leq \|\zeta\|_\star$ and $\|Y\|_\mathrm{op} \leq \|Y\|_F \leq \sqrt{2} \|\zeta\|_\star$. For $\theta \in \N (m, C)$ with $(m, C) \in \mathcal{S}_\star$, we denote $w := C^{-1} (\theta - m) = C^{-1/2} z$ for $z \in \N(0, I_{N_\theta})$. Using the standard Gaussian bounds, we obtain the bounds of moments of $w$,
    \begin{equation*}
        \E\|w\| \leq \sqrt{2} N_\theta, \quad \E\|w\|^2 \leq 2 N_\theta, \quad \E\|w\|^3 \leq 2\sqrt{2} N_\theta \sqrt{N_\theta + 2}, \quad \E\|w\|^4 \leq 4 N_\theta (N_\theta + 2).
    \end{equation*}

    For any smooth function $f$, we apply Bonnet--Price differentiation to obtain
    \begin{equation*}
        \frac{\dd}{\dd s} \E_{\N(m + sv, C + sY)} [f] = \langle v, \nabla_m \E[f] \rangle + \langle Y, \nabla_C \E[f] \rangle_F = \E \left[ v^\top \nabla f + \frac{1}{2} \Tr(Y \nabla^2 f) \right],
    \end{equation*}
    which is equivalent to 
    \begin{equation*}
        \frac{\dd}{\dd s} \E[f] = \E[\mathcal{L}_\zeta f], \qquad \mathcal{L}_\zeta = \sum_p v_p \partial_p + \frac{1}{2} \sum_{p,q} Y_{pq} \partial_p \partial_q,
    \end{equation*}
    for the second derivative, we have
    \begin{equation*}
        \frac{\dd^2}{\dd s^2} \E[f] = \E[\mathcal{L}_\zeta^2 f], \qquad \mathcal{L}_\zeta^2 = \sum_{p, q} v_p v_q \partial_p \partial_q + \sum_{p, q, r} v_p Y_{qr} \partial_p \partial_q \partial_r + \frac{1}{4} \sum_{p, q, r, s} Y_{pq} Y_{rs} \partial_p \partial_q \partial_r \partial_s.
    \end{equation*}
    
    We thus compute the first and second derivative of $\E[\nabla V]$ and $\E[\nabla^2 V]$ using Gaussian integration by parts:
\begin{comment}
    \begin{align*}
        \frac{\dd}{\dd s} \E[\nabla V]_i &= \E\left[ \sum_k v_k \partial_k (\partial_i V) + \frac{1}{2} \sum_{k, l} Y_{kl} \partial_k \partial_l (\partial_i V) \right] = \E \left[ [\nabla^2 V v]_i + \frac{1}{2} [\nabla^2 V (Yw)]_i \right], \\
        \frac{\dd }{\dd s} \E[\nabla^2 V]_{ij} &= \E\left[ \sum_k v_k \partial_k (\nabla^2 V_{ij}) + \frac{1}{2} \sum_{k, l} Y_{kl} \partial_k \partial_l (\nabla^2 V_{ij}) \right] \\
        &= \E\left[ \nabla^2 V_{ij} \sum_k v_k w_k + \frac{1}{2} \nabla^2 V_{ij} \sum_{k,l} Y_{kl} (w_k w_l - C_{kl}^{-1}) \right], \\
        \frac{\dd^2}{\dd s^2} \E[\nabla V]_i &= \E\left[ \sum_{p, q} v_p v_q \partial_p \partial_q (\partial_i V) + \sum_{p, q, r} v_p Y_{qr} \partial_p \partial_q \partial_r (\partial_i V) + \frac{1}{4} \sum_{p, q, r, s} Y_{pq} Y_{rs} \partial_p \partial_q \partial_r \partial_s (\partial_i V) \right] \\
        &= \E\left[ \sum_q \nabla^2 V_{iq} v_q (v^\top w) + \sum_q \nabla^2 V_{iq} ((v^\top w) (Yw)_q - (YC^{-1} v)_q) \right.\\
        &\qquad \left.+ \frac{1}{4} \sum_q \nabla^2 V_{iq} ((Yw)_q (w^\top Y w) - 2 (Y C^{-1} Y w)_q - (Yw)_q \Tr (YC^{-1})) \right], \\
        \frac{\dd^2 }{\dd s^2} \E[\nabla^2 V]_{ij} &= \E\left[ \sum_{p, q} v_p v_q \partial_p \partial_q (\nabla^2 V_{ij}) + \sum_{p, q, r} v_p Y_{qr} \partial_p \partial_q \partial_r (\nabla^2 V_{ij}) + \frac{1}{4} \sum_{p, q, r, s} Y_{pq} Y_{rs} \partial_p \partial_q \partial_r \partial_s (\nabla
        ^2 V_{ij}) \right] \\
        &= \E\left[ \nabla^2 V_{ij} \sum_{p, q} v_p v_q (w_p w_q - C_{pq}^{-1}) + \nabla^2 V_{ij} \sum_{p, q, r} v_p Y_{qr} (w_p w_q w_r - C_{pq}^{-1} w_r - C_{pr}^{-1} w_q - C_{qr}^{-1} w_p) \right. \\ 
        &\qquad \left. + \frac{1}{4} \nabla^2 V_{ij} \left( (w^\top Y w)^2 - 2\Tr (YC^{-1}) (w^\top Y w) - 4w^\top YC^{-1} Yw + \Tr (YC^{-1})^2 + 2\Tr (C^{-1} Y C^{-1} Y) \right) \right],
    \end{align*}
    these give
\end{comment}
    
    \begin{align*}
        \frac{\dd}{\dd s} \E[\nabla V] &= \E\left[ \nabla^2 V \left( v + \frac{1}{2} Yw \right) \right], \\
        \frac{\dd}{\dd s} \E[\nabla^2 V] &= \E\left[ \nabla^2V \left( v^\top w + \frac{1}{2} w^\top Y w - \frac{1}{2} \Tr (YC^{-1}) \right) \right],\\
        \frac{\dd^2}{\dd s^2} \E[\nabla V] &= \E\left[ \nabla^2 V \left( (v^\top w) v + (v^\top w) Yw - YC^{-1} v + \frac{1}{4} (w^\top Yw - \Tr (YC^{-1})) Yw - \frac{1}{2} YC^{-1} Yw \right) \right], \\
        \frac{\dd^2}{\dd s^2} \E[\nabla^2 V] &= \E\left[\nabla^2 V \left((v^\top w)^2 - v^\top C^{-1} v + (v^\top w)(w^\top Yw) - 2v^\top C^{-1}Yw - (v^\top w)\Tr(YC^{-1}) + \frac{1}{4}\Pi \right)\right],
    \end{align*}
    where $\Pi = (w^\top Yw)^2 - 2\Tr(YC^{-1})(w^\top Yw) - 4 w^\top YC^{-1}Yw + \Tr(YC^{-1})^2 + 2\Tr(C^{-1}YC^{-1}Y)$.
    Using $\|\nabla^2 V\|_\op \le \beta$, $\|\nabla^2 V\|_F \le \sqrt{N_\theta} \beta$, we bound each term:
    \begin{align*}
        \left\| \frac{\dd}{\dd s} \E[\nabla V] \right\| &\le \beta (1 + \sqrt{N_\theta}) \|\zeta\|_\star, \\
        \left\| \frac{\dd}{\dd s} \E[\nabla^2 V] \right\|_F &\le \sqrt{2}\beta (N_\theta^{3/2} + 2N_\theta) \|\zeta\|_\star, \\
        \left\| \frac{\dd^2}{\dd s^2} \E[\nabla V] \right\| &\le \sqrt{2}\beta \left(N_\theta \sqrt{N_\theta + 2} + 3N_\theta + 3\sqrt{N_\theta} + 2\right) \|\zeta\|_\star^2, \\
        \left\| \frac{\dd^2}{\dd s^2} \E[\nabla^2 V] \right\|_F &\le \sqrt{N_\theta} \beta \left(2N_\theta^2 + 4N_\theta \sqrt{N_\theta + 2} + 4N_\theta^{3/2} + 18N_\theta + 8\sqrt{N_\theta} + 8\right) \|\zeta\|_\star^2.
    \end{align*}
    Differentiating $F_1$ and $F_2$ twice along $a(s)$ and substituting the bounds above gives:
    \begin{align*}
        \| D^2 F_1 [\zeta, \zeta] \| &\le \frac{3}{2} \left\| \frac{\dd^2}{\dd s^2} \E[\nabla V] \right\| + 2\sqrt{2} \|\zeta\|_\star \left\| \frac{\dd}{\dd s} \E[\nabla V] \right\|, \\
        \| D^2 F_2 [\zeta, \zeta] \|_F &\le \frac{9}{4} \left\| \frac{\dd^2}{\dd s^2} \E[\nabla^2 V] \right\|_F + 6\sqrt{2} \|\zeta\|_\star \left\| \frac{\dd}{\dd s} \E[\nabla^2 V] \right\|_F + 4\beta \|\zeta\|_\star^2.
    \end{align*}
    By the definition of the Fisher--Rao norm, $\|D^2 F [\zeta, \zeta]\|_\star \le \|D^2 F_1 [\zeta, \zeta]\| + \frac{1}{\sqrt{2}} \|D^2 F_2 [\zeta, \zeta]\|_F$. Summing the bounds and evaluating the coefficients using $\sqrt{N_\theta + 2} \le \sqrt{3}\sqrt{N_\theta}$ ($N_\theta \geq 1$) yields:
    \begin{equation*}
        \|D^2 F [\zeta, \zeta]\|_\star \le 135 \beta N_\theta^{5/2} \|\zeta\|_\star^2.
    \end{equation*}
\end{proof}

\begin{theorem}\label{thm:loc-conv}
    Under Assumption \ref{assump:logconcave-smooth} and \eqref{eq:whitened}, we denote \begin{equation*}
        \delta := \min \left\{ \frac{1}{2\sqrt{2}}, \frac{1}{135 \beta N_\theta^{5/2} (4 + \Gamma)} \right\}, \qquad \Gamma = \min \left\{ \beta - 1, N_\theta \left( 3 + \frac{4}{\sqrt{\pi}} (1 + \log \kappa) \right) \right\}.
    \end{equation*}
    If the initialization satisfies $\|a_0 - a_\star\|_\star \leq \delta$, then the local convergence rate of the Gaussian natural gradient flow \eqref{ODEs:NG} reads
    \begin{equation*}
        \|a_t - a_\star\|_\star^2 \leq \exp \left(- \frac{t}{4 + \Gamma}\right) \|a_0 - a_\star \|_\star^2, \qquad t\geq 0.
    \end{equation*}
\end{theorem}

\begin{proof}
    We denote $\xi_t = a_t - a_\star$. Since the vector field of \eqref{ODEs:NG} satisfies $F(a_\star) = 0$ with $DF(a_\star) = J_\star$, we apply Taylor's theorem to obtain
    \begin{equation*}
        F(a_\star + \xi) = J_\star \xi + R(\xi), \qquad R(\xi) = \int_0^1 (1-s) D^2 F(a_\star + s\xi)[\xi, \xi] \dd s.
    \end{equation*}
    By Lemma \ref{lem:2ord-bound-NG}, since $\|\xi\|_\star \leq \delta \leq \frac{1}{2\sqrt{2}}$, then the trajectory remains in the convex neighborhood $\mathcal{S}_\star$, hence the remainder satisfies $\|R(\xi)\|_\star \le \frac{135}{2} \beta N_\theta^{5/2} \|\xi\|_\star^2$.
    
    We proceed to compute and bound the derivative of the squared norm 
    \begin{align*}
        \frac{\dd}{\dd t} \|\xi_t\|_\star^2 &= 2 \langle \xi_t, J_\star \xi_t \rangle_\star + 2 \langle \xi_t, R(\xi_t) \rangle_\star \\
        &\leq -\frac{2}{4 + \Gamma} \|\xi_t\|_\star^2 + 135 \beta N_\theta^{5/2} \|\xi_t\|_\star^3 \leq -\frac{1}{4 + \Gamma} \|\xi_t\|_\star^2,
    \end{align*}
    where the first inequality comes from Proposition \ref{prop:Jac-eig} and Lemma \ref{lem:2ord-bound-NG} while the second inequality comes from $\|\xi_t\|_\star \leq \delta \leq \frac{1}{135\beta N_\theta^{5/2} (4 + \Gamma)}$. This completes the proof by Gronwall's inequality.
\end{proof}

We then show the local convergence rate of the discrete algorithm with Riemannian distance \eqref{upd:Riem}.

\begin{theorem}
    Under Assumption \ref{assump:logconcave-smooth} and \eqref{eq:whitened}, we define
    \begin{align*}
        \delta_\mathrm{R} (\Delta t) := \min \left\{ \frac{1}{2\sqrt{2}}, \frac{1}{(4 + \Gamma) \cdot 148 \beta^2 N_\theta^{5/2}} \right\}, \quad \Gamma = \min \left\{ \beta - 1, N_\theta \left( 3 + \frac{4}{\sqrt{\pi}} (1 + \log \kappa) \right) \right\}.
        %C_\mathrm{R} &:= 135 \beta N_\theta^{5/2} + 2C_0 + \frac{C_0^2}{\sqrt{2}}, \qquad C_0 := \sqrt{2} + \frac{\beta - \alpha}{2} + \beta N_\theta
    \end{align*}
    Consider the discrete update with Riemannian distance \eqref{upd:Riem} with stepsize $0 < \Delta t \leq \left( \beta + \frac{1 - \alpha}{2} \right)^{-1}$, if the initialization satisfies $\|a_0 - a_\star\|_\star \leq \delta_\mathrm{R} (\Delta t)$, then the local convergence rate is
    \begin{equation*}
        \|a_N - a_\star\|_\star \leq \left( 1 - \frac{\Delta t}{2 (4 + \Gamma)} \right)^N \|a_0 - a_\star\|_\star.
    \end{equation*}
\end{theorem}

\begin{proof}
    We denote $\xi_n := a_n - a_\star$. By Proposition \ref{prop:upd-map-Riem}, for the discrete update \eqref{upd:Riem}, we have
    \begin{equation*}
        a_{n+1} = \Exp_{a_n} (-\Delta t \grad \mathcal{E} (a_n)) = \Exp_{a_n} (\Delta t F(a_n)).
    \end{equation*}
    The update of error admits the following expansion
    \begin{equation*}
        \xi_{n+1} = (I + \Delta t J_\star) \xi_n + R_{\mathrm{R}, \Delta t} (\xi_n).
    \end{equation*}
    By Proposition \ref{prop:Jac-eig}, the linear part satisfies \begin{equation*}
        \| (I + \Delta t J_\star) \xi \|_\star \leq \left( 1 - \frac{\Delta t}{4 + \Gamma} \right) \|\xi_n\|_\star
    \end{equation*}
    since the stepsize satisfies $0 < \Delta t \leq \left( \beta + \frac{1 - \alpha}{2} \right)^{-1}$.

    We proceed to look into the nonlinear term $R_{\mathrm{R}, \Delta t}$. From Lemma \ref{lem:2ord-bound-NG}, the vector field part contributes
    \begin{equation*}
        \Delta t \left\| \int_0^1 (1 - s) D^2 F(a_\star + s \xi_n) [\xi_n, \xi_n] \dd s \right\|_\star \leq \frac{\Delta t}{2} 135 \beta N_\theta^{5/2} \|\xi_n\|_\star.
    \end{equation*}
    The remaining part comes from the exponential covariance increment in the covariance component. We denote
    \begin{equation*}
        P(a) = I + C^{1/2} \E_{\rho_a} [\nabla_\theta \nabla_\theta \log \rho_\post] C^{1/2}
    \end{equation*}
    and bound the first variation of $P$ at $a_\star$ for $\xi = (u, X)$:
    \begin{equation*}
        DP(a_\star) [\xi] = -\left( X + T^*[u] + \frac{1}{2} H[X] \right) \Longrightarrow \| DP(a_\star) [\xi]\|_F \leq \left( \sqrt{2} + \frac{\beta - \alpha}{2} + \beta N_\theta \right) \|\xi\|_\star.
    \end{equation*}
    We expand the covariance update $C_{n+1} = C_n^{1/2} \exp(\Delta t P(a_n)) C_n^{1/2}$, since \begin{equation*}
        \exp (\Delta t P) = I + \Delta t P + \frac{\Delta t^2}{2} P^2 + O(\Delta t^3 \|P\|^3),
    \end{equation*} 
    we obtain
    \begin{align*}
        \|R_{\mathrm{R}, \Delta t} (\xi_n)\|_\star &\leq \frac{\Delta t}{2} 135 \beta N_\theta^{5/2} \|\xi_n\|_\star^2 + \frac{\Delta t}{2} \left(2 \left( \sqrt{2} + \frac{\beta - \alpha}{2} + \beta N_\theta \right)  + \frac{1}{\sqrt{2}} \left( \sqrt{2} + \frac{\beta - \alpha}{2} + \beta N_\theta \right)^2 \Delta t  \right) \|\xi_n\|_\star^2 \\
        &\leq \frac{\Delta t}{2} \left( 135 \beta N_\theta^{5/2} + 6\beta N_\theta + \frac{9 \beta^2 N_\theta^2}{\sqrt{2}} \right) \|\xi_n\|_\star^2 \leq \frac{\Delta t}{2} \cdot 148 \beta^2 N_\theta^{5/2} \|\xi_n\|_\star^2,
    \end{align*}
    where the continuous-field component from the first inequality comes from $\|\xi_n\|_\star \leq \delta_\mathrm{R} (\Delta t) \leq \frac{1}{2\sqrt{2}}$, which ensures $a_n \in \mathcal{S}_\star$. 
    
    Therefore
    \begin{align*}
        \|\xi_{n+1}\|_\star &\leq \| (I + \Delta t J_\star) \xi_n \|_\star + \|R_{\mathrm{R}, \Delta t} (\xi)\|_\star \\
        &\leq \left( 1 - \frac{\Delta t}{4 + \Gamma} \right) \|\xi_n\|_\star + \frac{\Delta t}{2} 148 \beta^2 N_\theta^{5/2} \|\xi_n\|_\star^2 \\
        &\leq \left( 1 - \frac{\Delta t}{2(4 + \Gamma)} \right) \|\xi_n\|_\star,
    \end{align*}
    where the last inequality comes from $\|\xi_n\|_\star \leq \delta_\mathrm{R} (\Delta t) \leq \frac{1}{(4 + \Gamma) 148 \beta^2 N_\theta^{5/2}}$. This completes the proof by iterating on $n$.
\end{proof}

We proceed to proceed to prove the local convergence rate of discretization scheme with KL divergence \eqref{upd:KL}.

\begin{theorem}
    Under Assumption \ref{assump:logconcave-smooth} and \eqref{eq:whitened}, we define
    \begin{align*}
        \delta_\mathrm{KL} (\Delta t) := \min \left\{ \frac{1}{2\sqrt{2}}, \frac{1}{(4 + 2 \Delta t + \Gamma) \cdot 216 \beta^2 N_\theta^{5/2}} \right\},\quad
        \Gamma = \min \left\{ \beta - 1, N_\theta \left( 3 + \frac{4}{\sqrt{\pi}} (1 + \log \kappa) \right) \right\}.
    \end{align*}
    Consider the discrete update with KL divergence \eqref{upd:KL} with the following local norm
    \begin{equation*}
        \|(u, X)\|_{\star, \Delta t}^2 = \|u\|^2 + \frac{1 + \Delta t}{2} \|X\|_F^2, \qquad (u, X) \in \R^{N_\theta} \times \Sym(N_\theta),
    \end{equation*}
    if the stepsize satisfies $0 < \Delta t \leq \left( \beta + \frac{1 - \alpha}{2} \right)^{-1}$ and the initialization satisfies $\|a_0 - a_\star\|_{\star, \Delta t} \leq \delta_{\mathrm{KL}} (\Delta t)$, then the local convergence guarantee reads
    \begin{equation*}
        \|a_N - a_\star\|_{\star, \Delta t} \leq \left( 1 - \frac{\Delta t}{2 (4 + 2 \Delta t + \Gamma)} \right)^N \|a_0 - a_\star\|_{\star, \Delta t}.
    \end{equation*}
\end{theorem}

\begin{proof}
    We write $\xi_n:= a_n - a_\star = (u_n, X_n)$, and characterize \eqref{upd:KL} into linear part and non-linear part
    \begin{equation*}
        \xi_{n+1} = M_{\Delta t, \star}^\KL (\xi_n) + R_{\Delta t}^\KL (\xi_n).
    \end{equation*}
    We first calculate the linearized map using the similar technique in Proposition \ref{prop:Jac-eig}:
    \begin{equation*}
        M_{\Delta t, \star}^\KL (u, X) = \left( (1 - \Delta t) u - \frac{\Delta t}{2} T[X],\ \frac{1}{1 + \Delta t} \left( X - \Delta t T^*[u] - \frac{\Delta t}{2} H[X] \right) \right),
    \end{equation*}
    which is self-adjoint with respect to the weighted inner product at equilibrium because $$\langle M_{\Delta t, \star}^\KL (u_1, X_1), (u_2, X_2) \rangle_{\star, \Delta t} = \langle (u_1, X_1), M_{\Delta t, \star}^\KL (u_2, X_2) \rangle_{\star, \Delta t},$$ where the inner product at equilibrium is defined by
    \begin{equation*}
        \langle (u_1, X_1), (u_2, X_2) \rangle_{\star, \Delta t} := u_1^\top u_2 + \frac{1 + \Delta t}{2} \Tr (X_1 X_2).
    \end{equation*}
    We observe that
    \begin{equation*}
        \|(u, X)\|_{\star, \Delta t}^2 - \langle M_{\Delta t, \star}^\KL (u, X), (u, X) \rangle_{\star, \Delta t} = \Delta t Q(u, X),
    \end{equation*}
    hence the largest eigenvalue $\rho_\star^\KL$ of the linearized Jacobian matrix $M_{\Delta t, \star}^\KL$ satisfies
    \begin{align*}
        1 - \rho_\star^\KL &= \Delta t \min_{(u, X) \neq 0} \frac{Q(u, X)}{\|u\|^2 + \frac{1 + \Delta t}{2} \|X\|_F^2} \geq \frac{\Delta t}{4 + 2\Delta t + \Gamma},
    \end{align*}
    which is because $g(\nu) := 1 - \nu (4 + 2 \Delta t + \Gamma) + 2 (1 + \Delta t) \nu^2 \geq 0$ at $\nu = \frac{1}{4 + 2 \Delta t + \Gamma}$.
    We thus have
    \begin{equation*}
        \|M_{\Delta t, \star}^\KL (\xi)\|_{\star, \Delta t} \leq \left( 1 - \frac{\Delta t}{4 + 2\Delta t + \Gamma} \right) \|\xi\|_{\star, \Delta t}
    \end{equation*}
    since $0 < \Delta t \leq \left( \beta + \frac{1 - \alpha}{2} \right)^{-1}$ and hence $\langle M_{\Delta t, \star}^\KL (\xi), \xi \rangle_{\star, \Delta t}\geq 0$.

    We proceed to control the nonlinear part. Along $a = a_\star + s\xi \in \R^{N_\theta} \times $, we compute
    \begin{equation*}
        \left. \frac{\dd}{\dd s} \E_{\rho_{a + s \xi}}[\nabla^2 V] \ \right|_{s=0} = T^*[u] + \frac{1}{2} H[X] := A_1 (\xi)
    \end{equation*}
    and bound \begin{equation*}
        \left\| A_1 (\xi) \right\|_F \leq \left( \frac{\beta - \alpha}{2} + \sqrt{2} \beta N_\theta \right) \|\xi\|_\star.
    \end{equation*}
    From Lemma \ref{lem:2ord-bound-NG}, 
    \begin{equation*}
        \| D^2 F[\zeta, \zeta]\|_{\star, \Delta t} \leq \sqrt{1 + \Delta t} \| D^2 F[\zeta, \zeta]\|_{\star} \leq 135 \sqrt{1 + \Delta t} \beta N_\theta^{5/2} \|\zeta\|_{\star}^2 \leq 135\sqrt{2} \beta N_\theta^{5/2} \|\zeta\|_{\star, \Delta t}^2.
    \end{equation*}
    By expanding the covariance update $C_{n+1} = (1 + \Delta t) (C_n + \Delta t \E_{\rho_{a_n}}[\nabla^2 V])^{-1}$ and combining the second-order correction bound with continuous-time high-order remainder, we obtain
    \begin{align*}
        \|R_{\Delta t}^\KL (\xi_n) \|_{\star, \Delta t} &\leq \frac{\Delta t}{2} \left( 135 \sqrt{2} \beta N_\theta^{5/2} + 2\sqrt{2} + 8\left( \frac{\beta - \alpha}{2} + \sqrt{2} \beta N_\theta \right) + \sqrt{2} \left( \frac{\beta - \alpha}{2} + \sqrt{2} \beta N_\theta \right)^2 \right) \|\xi_n\|_{\star, \Delta t}^2 \\
        &\leq \frac{\Delta t}{2} \left( 135 \sqrt{2} \beta N_\theta^{5/2} + 2\sqrt{2} + 16 \beta N_\theta + 4\sqrt{2} \beta^2 N_\theta^2 \right) \|\xi_n\|_{\star, \Delta t}^2 \leq \frac{\Delta t}{2} 216 \beta^2 N_\theta^{5/2} \|\xi_n\|_{\star, \Delta t}^2,
    \end{align*}
    where the vector field component from the first inequality comes from Lemma \ref{lem:2ord-bound-NG} and $\|\xi_n\|_\star \leq \delta_\KL (\Delta t) \leq \frac{1}{2\sqrt{2}}$, which ensures that $a_n \in \mathcal{S}_\star$.

    Therefore, as the stepsize is chosen to satisfy $0 < \Delta t \leq \left( \beta + \frac{1 - \alpha}{2} \right)^{-1}$, then
    \begin{align*}
        \|\xi_{n+1} \|_{\star, \Delta t} &\leq \|M_{\Delta t, \star}^\KL (\xi)\|_{\star, \Delta t} + \|R_{\Delta t}^\KL (\xi_n) \|_{\star, \Delta t} \\
        &\leq \left( 1 - \frac{\Delta t}{4 + 2\Delta t + \Gamma} \right) \|\xi\|_{\star, \Delta t} + \frac{\Delta t}{2} 216 \beta^2 N_\theta^{5/2} \|\xi_n\|_{\star, \Delta t}^2 \\
        &\leq \left( 1 - \frac{\Delta t}{2(4 + 2\Delta t + \Gamma)} \right) \|\xi\|_{\star, \Delta t},
    \end{align*}
    where the last inequality comes from $\|\xi_n\|_{\star, \Delta t} \leq \delta_\KL (\Delta t) \leq \frac{1}{(4 + 2 \Delta t + \Gamma) \cdot 216 \beta^2 N_\theta^{5/2}}$. This completes the proof by iterating on $n$.
\end{proof}

