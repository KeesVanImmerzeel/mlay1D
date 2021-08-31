
# Analytic steady-state, one-dimensional, leaky multi-aquifer groundwater model

C.H. van Immerzeel

31/8/2021

![Head plot](https://user-images.githubusercontent.com/16401251/127302532-25676075-f24f-4cff-a928-3d9f7801a28a.png)
![Lateral flux plot](https://user-images.githubusercontent.com/16401251/127303441-a00f6c81-5490-4db0-ae7e-936b49cd2a37.png)
![Seepage plot](https://user-images.githubusercontent.com/16401251/127303673-8e56bb3c-8060-4e89-99ef-34c1bd45474f.png)

## Link to app
<https://sweco.shinyapps.io/mlay1D/>

## Source code
R-source code of the app:

<https://github.com/KeesVanImmerzeel/mlay1D>

Here, also Matlab/Octave and Python routines of the numerical solution is available.

## Numerical stability
Exotic values for kD and c, as well as large distances between intersection points may lead to numerical instability. In that case, no results or plots are presented. Instead, an error message appears.

In order to fix this, you can try to divide the sections in parts using the field "f" in the Control panel. Typical values range from f=2 to f=100. However, doing so increases the calculation time. Try increasing f-values until numerical stability (a result) is reached.

## References
- Olsthoorn, T.N. (2000). Eendimensionale stroming in meerlagensystemen. Stromingen, jaargang 6, nummer 2, blz. 13-24.
  <https://edepot.wur.nl/10049>
