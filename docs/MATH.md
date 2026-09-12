# Mathematics — wavefield

## Derivation sketch of the stencil

Replace $\partial_{tt}u$ and $\nabla^2 u$ by centered second differences of order $O(\Delta t^2)$ and $O(\Delta x^2)$. Collecting terms yields the leapfrog update used in all three language ports.

## CFL from von Neumann analysis

Assume a plane-wave mode $u^n_{i,j} = \xi^n e^{I(k_x i\Delta x + k_y j\Delta y)}$. Requiring $|\xi| \le 1$ for all wavevectors yields the isotropic bound $c\Delta t/\Delta x \le 1/\sqrt{2}$.

## Energy proxy

A discrete energy that is approximately conserved under Dirichlet boundaries (exact conservation for the continuous Hamiltonian; discrete form is a monitor, not a theorem):

$$
\mathcal{E}^n \approx \sum_{i,j}\Bigl[\bigl(u^n_{i,j}-u^{n-1}_{i,j}\bigr)^2 + c^2\bigl((\delta_x u)^2 + (\delta_y u)^2\bigr)\Bigr]
$$

Mur boundaries intentionally **drain** $\mathcal{E}$ as waves exit.
