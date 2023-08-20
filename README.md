# Implementation of TASE-RK-methods (by Dajana Conte, Giovanni Pagano, Beatrice Paternoster, Department of Mathematics, University of Salerno, Italy) 
TASE-RK-methods [Calvo, Montijano, Rández, “A note on the stability of time-accurate and highly-stable explicit operators for stiff differential equations”, J. Comput. Phys. 436 (2021), 110316] are linearly implicit numerical schemes for solving stiff initial value problems of the type $y'(t)=f(t,y(t))$,  $y(t_0)=y_0$, $t \in [t_0,t_e]$, $f:\mathbb{R}\times \mathbb{R}^d \rightarrow \mathbb{R}^d$ (1).

Such numerical schemes require the solution of s*p (s = number of stages, p = order of the method) linear systems per step involving the Jacobian matrix of the problem. The proposed implementation allows two options:
- update the Jacobian of the problem at each step;
- fix a constant approximation of the Jacobian for all the time integration.

In the latter case, the computational cost of the method drops drastically, while the order of consistency is preserved. However, if the chosen Jacobian approximation 'is not good', stability issues may be encountered.

The codes reported here for the implementation of TASE-RK methods are written in Matlab language. 
These codes are explained in detail in the manuscript “A MATLAB implementation of TASE-RK methods”, by Dajana Conte, Giovanni Pagano, Beatrice Paternoster (Department of Mathematics, University of Salerno, Italy), accepted for publication in the Journal of Approximation Software (JAS).

# Input and output arguments

Input arguments

• N - integer scalar

Number of (equally spaced) discrete time intervals into which the user decides to subdivide the continuous time grid $[t_0,t_e]$.

• tspan - double array

Row vector of length two containing the first and last grid points $t_0$, $t_e$, respectively.

• y0 - double array

Column vector with the initial condition $y_0 \in \mathbb{R}^d$.

• Fun - function handle

Function which returns the vector field $f$, evaluated at $(t,y)$, of the problem (1) that the user wants to solve; $t \in  \mathbb{R}$, $y \in  \mathbb{R}^d$, constitute the input arguments of Fun, and the column vector $f(t,y) \in  \mathbb{R}^d$ is the output.

• Jac - function handle

Function which returns the Jacobian matrix $J_f$ of the problem (1) that the user wants to solve, evaluated at a point $(t,y)$, or a suitable fixed approximation of $J_f$; in the first case, $t \in \mathbb{R}$, $y \in \mathbb{R}^d$, constitute the input arguments of Jac, and the matrix $f_y(t,y) \in \mathbb{R}^{d,d}$ is the output.

• Method - integer array

Array with the TASE-RK methods to apply; in particular:

- 20 corresponds to the TASE-RK with s = p = 2, using the midpoint rule as under-
lying explicit RK;

- 30 corresponds to the TASE-RK with s = p = 3, using the Ralston’s method as
underlying explicit RK;

- 40 corresponds to the TASE-RK with s = p = 4, using the Kutta’s method as un-
derlying explicit RK.

For example, if we want to use methods 20 and 30, then Method=[20,30]. We will
later show an example where we apply all the TASE-RK methods using the same main
program.
• jacup - integer scalar
Parameter which is equal to 0 if the user wants to use constant Jn, 1 otherwise.
Output arguments
• yT - double array

Column vector of length d with the numerical solution computed by the chosen TASE-
RK method at last grid point te.

• y - double array
Matrix of size d ×(N +1) having in column n+1 the numerical solution yn compute by
the chosen TASE-RK method.
• t - double array
Row vector of length N+1 with all the discrete grid points {tn =t0+nh;n = 0,...,N;tN =
te}.

D. Conte, G. Pagano, B. Paternoster 8/21

A MATLAB implementation of TASE-RK methods

• CPUtime - double scalar
Total CPU time in seconds taken by the chosen TASE-RK method.

# Test case

# References
