
# Analytic steady-state, one-dimensional, leaky multi-aquifer groundwater model

C.H. van Immerzeel

2/9/2021

![Head plot](https://user-images.githubusercontent.com/16401251/131862426-cd709642-6294-474a-b7f2-b94bcc074b7f.png)
![Lateral flux plot](https://user-images.githubusercontent.com/16401251/131862604-399d08f8-b54b-4362-916c-2400ff77e806.png)
![Seepage plot](https://user-images.githubusercontent.com/16401251/131862832-a8e9f988-0015-436d-b2d3-8b4b1874d526.png)

## Link to the app
<https://sweco.shinyapps.io/mlay1D/>

## Source code
R-source code of the app:

<https://github.com/KeesVanImmerzeel/mlay1D>

An example is included if spreadsheet data is to be used as input.

Also Matlab/Octave and Python routines of the numerical solution is made available.

## Numerical stability
Exotic values for kD and c, as well as large distances between intersection points may lead to numerical instability. In that case, no results or plots are presented. Instead, an error message appears.

In order to fix this, you can try to divide the sections in parts using the field "f" in the Control panel. Typical values range from f=2 to f=100. However, doing so increases the calculation time. Try increasing f-values until numerical stability (a result) is reached.


## References
- Olsthoorn, T.N. (2000). Eendimensionale stroming in meerlagensystemen. Stromingen, jaargang 6, nummer 2, blz. 13-24.
  <https://edepot.wur.nl/10049>
