# Approximation-of-Landau-Constant
Programiranje 2 project
The goal of the project is to calculate the Landau constant to the best precision so far. More about the constant and the algorithm can be found in article https://lmcs.episciences.org/1189/pdf.

So far, we have implemented Dyadic numbers for efficient computation, a psi function to represent holomorphic functions and a struct of holomorphic functions. We will then take images of the unit complex disk and calculate the area of the largest disk inside the image. We have also implemented epsilon-covering grids that will help us calculate this. 


- [x] fix dyadic operations implementations (Matija) 
- [x] implement psi function as a recursive function (Luka)
- [x] implement sequence of m_i so they are bounded by some constant? (Luka)
- [x] implement series to represent holomorphic functions in rust (Luka)
- [ ] indefinite integral of a holomorphic function (Luka)
- [ ] struct for holomorphic functions (derivative, integral, power series representation) (Luka)
- [x] epsilon-covering grid (Matija)
- [ ] represent psi function with complex functions
