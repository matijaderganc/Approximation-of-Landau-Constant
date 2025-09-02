# Approximation-of-Landau-Constant
Programiranje 2 project
The goal of the project is to calculate the Landau constant to the best precision so far. More about the constant and the algorithm can be found in article https://lmcs.episciences.org/1189/pdf. Put shortly, the constant gives us the largest circle to fit inside the image of the unit circle with any normalised holomorphic function.

So far, we have implemented Dyadic numbers for efficient computation, a psi function to represent holomorphic functions and a struct of holomorphic functions. We will then take images of the unit complex disk and calculate the area of the largest disk inside the image using an EDT algorithm. We have also implemented epsilon-covering grids that will help us calculate this. The computations will be parallelized with Tokyo to save time. 

- [x] implement psi function as a recursive function (Luka)
- [x] implement sequence of m_i so they are bounded by some constant? (Luka)
- [x] implement series to represent holomorphic functions in rust (Luka)
- [x] Add 'fast power' method for ComplexDyadic (Luka)
- [x] indefinite integral of a holomorphic function (Luka)
- [x] struct for holomorphic functions (derivative, integral, power series representation) (Luka)
- [x] epsilon-covering grid (Matija)
- [x] represent psi function with complex functions
- [x] parallelization of distance calculations
- [ ] main loop of words
- [ ] Dashboard (Luka)
