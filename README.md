# Racket_NumericalMethods
A collection of Numerical Methods implemented in Racket. Includes Root Finders, ODE Solvers, Linear System Solvers, etc. as well as Applications/Demos.
The different methods are mostly split into individual files.

The "Experimenting" files include demos to show the different methods through examples.

Restricted_3_Body.rkt has a few examples of Restricted Three Body Problems. These are designed to test Numerical ODE's, since each set of initial conditions is calculated to a very high accuracy to be precisely periodic (returns to start). Thus, one may look at the solutions from particular ODE solvers and compare their accuracy based on how close they return to the start.
