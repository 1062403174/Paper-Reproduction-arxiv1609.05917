# Reproducing Figure 3 of *The Radial Acceleration Relation in Rotationally Supported Galaxies*

This repository reproduces Figure 3 of the paper *The Radial Acceleration Relation in Rotationally Supported Galaxies* using publicly available SPARC data.
The goal of this project is not to recreate the figure pixel by pixel, but to recover its main physical structure: the tight relation between the observed centripetal acceleration and the acceleration predicted by the baryonic mass distribution.

## 1. Paper Summary
Using the SPARC galaxy sample, the original paper found an extremely tight empirical relation between the observed centripetal acceleration, \(g_{\rm obs}\), and the baryon-predicted acceleration, \(g_{\rm bar}\).
The key result is that galaxy rotation curves are not arbitrarily shaped by dark matter alone. Instead, the baryonic distribution strongly constrains the observed acceleration, producing a remarkably narrow radial acceleration relation (RAR).


## 2. Reproduction Result
The main reproduction figure is provided in this repository.


## 3. How to Read the Figure
Figure 3 has two panels.
The upper panel shows the relation between the observed centripetal acceleration \(g_{\rm obs}\) and the acceleration predicted from the baryonic distribution \(g_{\rm bar}\). The black solid line is the empirical fit to all unbinned data points. The lower panel shows the residuals of each data point with respect to that fitted curve.
The most important feature in the upper panel is the overall shape. At low acceleration, the data clearly deviate from the line \(y=x\). As \(g_{\rm bar}\) increases, the relation gradually approaches the Newtonian limit and becomes close to \(y=x\) in the high-acceleration regime. This implies that baryons alone are not sufficient to explain the total gravitational acceleration in the low-acceleration regime, where an additional gravitational component is needed. In the high-acceleration regime, baryons become sufficient to explain most of the acceleration, so the relation returns toward the one-to-one Newtonian expectation. This is one reason why the result is so important: it suggests that dark matter cannot shape galaxy rotation curves arbitrarily, and that the baryonic distribution strongly constrains the final rotation pattern.
The lower panel shows the residuals of each data point relative to the black fitted curve in the upper panel. The residuals do not scatter widely at random. Instead, they form a narrow band centered around zero, indicating that the correspondence between \(g_{\rm obs}\) and \(g_{\rm bar}\) is highly stable, and that most of the remaining deviation can be attributed to observational uncertainties and model-parameter uncertainties.


## 4. Differences from the Original Figure
This reproduction captures the main structure of the original figure, but it is not identical to the paper’s version.
- **Sample selection is not fully identical.**  
  Although this reproduction uses the same SPARC data source, the original Figure 3 was based on a filtered sample of 153 galaxies and 2693 data points. I did not fully replicate all of the paper’s exclusion criteria, so the final sample size and composition differ.
- **The mass-to-light ratio and baryonic component handling may differ slightly.**  
  In Figure 3, the paper adopted fixed near-infrared \(M/L\) values for the stellar disk and bulge. Any small difference in these parameters, in the component-scaling implementation, or in the fitting details can shift the fitted characteristic acceleration \(g_\dagger\).  
  The paper reports \(g_\dagger = [1.20 \pm 0.02\,(\mathrm{random}) \pm 0.24\,(\mathrm{systematic})] \times 10^{-10}\,\mathrm{m\,s^{-2}}\), while my reproduction gives \(g_\dagger \approx 1.07 \times 10^{-10}\,\mathrm{m\,s^{-2}}\).
Here, \(g_\dagger\) is the characteristic acceleration scale in the RAR fit. It sets the transition between the low-acceleration regime and the Newtonian regime. When \(g_{\rm bar} \gg g_\dagger\), the relation approaches \(g_{\rm obs} \approx g_{\rm bar}\).


## 5. Scientific Significance
Before this result, it was natural to assume that baryonic mass and dark matter content could vary relatively independently from galaxy to galaxy. In that picture, one might expect galaxies with little baryonic mass but a large dark matter contribution, where the total acceleration is dominated by dark matter.
This paper shows that the real Universe is much less arbitrary than that. Even if dark matter still exists, the observed rotational acceleration is strongly locked to the baryonic distribution. The result implies that any successful galaxy-formation model must explain why baryons and dark matter are so tightly linked, and why the scatter in this relation is so small.


## 6. Data and Code
### Data source
The data were downloaded from the official SPARC website:
- SPARC homepage: https://astroweb.case.edu/SPARC/
- RAR data file: `RAR.mrt`
### Reproduction workflow
1. Load the SPARC RAR data file.
2. Compute \(g_{\rm obs}\) and \(g_{\rm bar}\).
3. Fit the empirical RAR relation with an orthogonal-distance regression (ODR) model.
4. Plot the upper density map, binned means, fit curve, residual panel, and residual histogram.
### Key parameter choices
- `U_DISK = 0.50`
- `U_BUL = 0.70`
### Code
The main reproduction script is provided in this repository.

## 7. Repository Structure
```text
.
├── data/            # Downloaded SPARC data
├── figures/         # Original and reproduced figures
├── code/            # Reproduction scripts
├── paper/           # Original paper screenshots / notes
└── README.md
