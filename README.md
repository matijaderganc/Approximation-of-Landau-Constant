# Approximation-of-Landau-Constant
Goal of the project is to compute the Landau's constant for any accuracy bound. The implemented algorithm is based on the paper [here](https://lmcs.episciences.org/1189/pdf).

Kaj bi bilo potrebno implementirati?
Torej, funkicja psi neskončno deluje po naslednjem principu:
imamo m1, m2, ... pozitivna diadična števila, tako da potenčna vrsta s temi koeficienti konvergira na enotskem disku. (od kje dobimo ta števila še ne vem)

imamo t1, t2, .. zaporedje, kjer se vsako naravno število pojavi neskončnokrat (je en primer v članku, nevem če ga bomo res uporabili)

nato pač vsako kompleksno v potenčni vrsti število predstavimo s funkcijo psi, pri čemer je različen tako začeten interval kot beseda

nato narediti zaporedje b1, b2, ... za katerega držijo dovolj lepe stvari

od nekje dobit cifro e

### TODO

- [x] fix dyadic operations implementations (Matija) 
- [ ] implement psi function as a recursive function (Luka)
- [x] implement sequence of m_i so they are bounded by some constant? (Luka)
- [x] implement series to represent holomorphic functions in rust (Luka)
- [ ] indefinite integral of a holomorphic function (Luka)
- [ ] struct for holomorphic functions (derivative, integral, power series representation) (Luka)
- [ ] epsilon-covering grid (Luka)
- [ ] represent psi function with complex functions
- [ ] Get value of constant e 

  
