# Authors
The authors are: Dajana Conte, Giovanni Pagano, Beatrice Paternoster (Department of Mathematics, University of Salerno, Italy).

Some information about the authors:
- Dajana Conte webpage https://docenti.unisa.it/020280/home; email dajconte@unisa.it;
- Giovanni Pagano (corresponding) webpage https://www.researchgate.net/profile/Giovanni-Pagano-2; email gpagano@unisa.it;
- Beatrice Paternoster webpage https://docenti.unisa.it/000793/home; email beapat@unisa.it.

# Implementation of TASE-RK methods  
TASE-RK (Runge-Kutta) methods [https://www.sciencedirect.com/science/article/pii/S0021999121002114] are linearly implicit numerical schemes for solving stiff initial value problems of the type $y'(t)=f(t,y(t))$,  $y(t_0)=y_0$, $t \in [t_0,t_e]$, $f:\mathbb{R}\times \mathbb{R}^d \rightarrow \mathbb{R}^d$ (1).

Such numerical schemes require the solution of $s \times p$ ($s$ = number of stages, $p$ = order of the method) linear systems per step involving the Jacobian matrix of the problem. The proposed implementation allows two options:
- update the Jacobian of the problem at each step;
- fix a constant approximation of the Jacobian for all the time integration.

In the latter case, the computational cost of the method drops drastically, while the order of consistency is preserved. However, if the chosen Jacobian approximation 'is not good', stability issues may be encountered.

The codes reported here for the implementation and application of TASE-RK methods are written in Matlab language. 
These codes are explained in detail in the manuscript “A MATLAB implementation of TASE-RK methods”, by Dajana Conte, Giovanni Pagano, Beatrice Paternoster (Department of Mathematics, University of Salerno, Italy), accepted for publication in the Journal of Approximation Software.

# Input and output arguments
The function TASERK.m applies the chosen TASE-RK method to solve a selected initial value problem. We describe the input and output arguments of TASERK.m below.

• Input arguments

- N - integer scalar

Number of (equally spaced) discrete time intervals into which the user decides to subdivide the continuous time grid $[t_0,t_e]$.

- tspan - double array

Row vector of length two containing the first and last grid points $t_0$, $t_e$, respectively.

- y0 - double array

Column vector with the initial condition $y_0 \in \mathbb{R}^d$.

- Fun - function handle

Function which returns the vector field $f$, evaluated at $(t,y)$, of the problem (1) that the user wants to solve; $t \in  \mathbb{R}$, $y \in  \mathbb{R}^d$, constitute the input arguments of Fun, and the column vector $f(t,y) \in  \mathbb{R}^d$ is the output.

- Jac - function handle

Function which returns the Jacobian matrix $J_f$ of the problem (1) that the user wants to solve, evaluated at a point $(t,y)$, or a suitable fixed approximation of $J_f$; in the first case, $t \in \mathbb{R}$, $y \in \mathbb{R}^d$, constitute the input arguments of Jac, and the matrix $f_y(t,y) \in \mathbb{R}^{d,d}$ is the output.

- Method - integer array

Array with the TASE-RK methods to apply; in particular:

• 20 corresponds to the TASE-RK with $s = p = 2$, using the midpoint rule as underlying explicit RK;

• 30 corresponds to the TASE-RK with $s = p = 3$, using the Ralston’s method as underlying explicit RK;

• 40 corresponds to the TASE-RK with $s = p = 4$, using the Kutta’s method as underlying explicit RK.

For example, if we want to use methods 20 and 30, then Method=[20,30]. The main program exampleBurgers.m, e.g., applies all the TASE-RK methods for the solution of the Burgers' equation.

- jacup - integer scalar
  
Parameter which is equal to 0 if the user wants 'constant Jacobian', 1 otherwise.

• Output arguments

- yT - double array

Column vector of length $d$ with the numerical solution computed by the chosen TASE-RK method at last grid point $t_e$.

- y - double array
  
Matrix of size $d \times (N +1)$ having in column $n+1$ the numerical solution $y_n$ computed by the chosen TASE-RK method.

- t - double array
  
Row vector of length $N+1$ with all the discrete grid points $[t_n =t_0+n h;n = 0,\dots,N;t_N =t_e]$.

- CPUtime - double scalar
  
Total CPU time in seconds taken by the chosen TASE-RK method.

# Test case
Here, we test the code for an example case. In particular, we show the application of the function TASERK.m in solving the well known Euler’s problem [https://en.wikipedia.org/wiki/Euler%27s_equations_(rigid_body_dynamics)]. To do this, just run the main code exampleEuler.m. This code uses the functions funEuler.m, jacEuler.m, jacEulerfix.m, and of course TASERK.m.

As you can see from line 3 of exampleEuler.m, in this case we apply only the TASE-RK method of order $p=s=4$. Furthermore, from line 5, note that we fix a constant approximation of the Jacobian (i.e. the exact Jacobian fixed at the first time grid point $t_0$), thus deciding not to update it at each step. From line 12 to line 16 we set the initial conditions, and from line 18 to line 21 we compute a reference solution in order to be able to calculate the TASE-RK error. In line 23 we apply the selected TASE-RK method.

The outputs are as follows:

errT_TRK =
   3.3776e-08;
   
CPUtimeTRK =
   6.4062e-01.

The first corresponds to the relative error committed by the TASE-RK method with respect to the reference solution using the infinity norm. The second corresponds to the total CPU time used by the TASE-RK method.

Similarly, running the main code exampleBurgers.m, which uses the functions funBurgers.m, jacBurgers.m, jacBurgersfix.m, and of course TASERK.m, it is possible to solve the famous Burgers' equation [https://en.wikipedia.org/wiki/Burgers%27_equation] with all the TASE-RK methods (i.e. of order $p=s=2,3,4$), after a spatial semi-discretization of the problem (method of lines) done by means of finite differences of order four.


# References
- M. Calvo, J. I. Montijano, L. Rández, “A note on the stability of time-accurate and highly-stable explicit operators for stiff differential equations”, J. Comput. Phys. 436 (2021), 110316.
- D. Conte, G. Pagano, B. Paternoster, “A MATLAB implementation of TASE-RK methods”, J. Approx. Softw. (2023), to appear.
