# Implementation of TASE-RK-methods (by Dajana Conte, Giovanni Pagano, Beatrice Paternoster, Department of Mathematics, University of Salerno, Italy) 
TASE-RK-methods [Calvo, Montijano, Rández, “A note on the stability of time-accurate and highly-stable explicit operators for stiff differential equations”, J. Comput. Phys. 436 (2021), 110316] are linearly implicit numerical schemes for solving stiff initial value problems. 
Such numerical schemes require the solution of s*p (s = number of stages, p = order of the method) linear systems per step involving the Jacobian matrix of the problem. The proposed implementation allows two options:
- update the Jacobian of the problem at each step;
- fix a constant approximation of the Jacobian for all the time integration.

In the latter case, the cost of the method drops drastically, while the order of consistency is preserved. However, if the chosen Jacobian approximation is not 'good', stability issues may be encountered.

The codes reported here for the implementation of TASE-RK methods are written in Matlab language. 
These codes are explained in detail in the manuscript “A MATLAB implementation of TASE-RK methods”, by Dajana Conte, Giovanni Pagano, Beatrice Paternoster (Department of Mathematics, University of Salerno, Italy) accepted for publication in the Journal of Approximation Software (JAS).

# Input and output arguments

# Test case

# References
