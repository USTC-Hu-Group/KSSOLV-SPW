# Correction notice

## 1. Translation-phase sign in SPW construction

**Date: 2026-09-20**

### Summary

We identified and corrected a sign error in the translation-phase factor in a previously uploaded version of `genSPWs_allirep_k.p`.

**This error affected only the code version that was uploaded to the public repository. The calculations reported in the manuscript were performed using the correct implementation. Therefore, the numerical results and conclusions presented in the manuscript are unaffected by this correction.**

The current version of the repository has been updated and is consistent with the implementation used to generate the results reported in the manuscript.

### Code error

For the symmetry-operation convention adopted in the implementation, a translation by $\boldsymbol{\tau}$ acts on a plane wave as

$$
e^{-i\mathbf{k}\cdot\boldsymbol{\tau}}
e^{i\mathbf{k}\cdot\mathbf{r}}.
$$

In the previously uploaded version of `genSPWs_allirep_k.p`, the translation-phase factor was implemented with the opposite sign,
$$
e^{+i\mathbf{k}\cdot\boldsymbol{\tau}},
$$

because of a mismatch between the symmetry-operation convention used in the implementation and the convention used when interpreting the stored BCS irreducible-representation data.

This was an implementation/version mismatch in the uploaded code and has now been corrected.

### Potential impact

The incorrect sign can produce plane-wave combinations that do not transform according to the intended irreducible representations. Consequently, calculations performed using the affected uploaded version may produce incorrect numerical results.

Checks of basis dimension and orthonormality alone do not necessarily detect this problem. Users who performed calculations using the affected version should rerun potentially affected cases with the corrected implementation.

Again, this issue does **not** affect the numerical results reported in the associated manuscript, which were obtained using the correct implementation.