# CrossSections
Collate, fit, and plot dimuon cross-sections from Whizard 3.0.3

## python3 CrossSections.py
Read in data files and fit
i)  ALR vs ECM and 
ii) unpolarized cross-section (in units of R) vs ECM.

Result is plots of ALR, xsLR, xsRL, xsU, Rmumu vs ECM 
with the superimposed parametrizations. (xsLR, xsRL, xsU models 
are inferred from the ALR and Rmumu fits using the QED point-like 
cross-section).

Current fits are a polynomial with 3 free parameters to ALR, and 
a fit to Rmumu using a BW and polynomial terms with 4 free parameters.
See mymodels.py, fitconfig,py and fitting.py for details.

The uncertainties in the Whizard cross-section have been scaled up a little 
using the chi-squared amongst the various iterations. I have seen some 
occurrences of the Whizard cross-section integration being quite discrepant, so 
I think there is still an expectation of some outliers and not particularly 
well modeled uncertainties. The run tests, while nominally less 
powerful (than chi-squared) look encouraging.

ALR fit.
Chi-squared / dof = 38.0/24, pvalue = 3.4%
Run-test p-value = 88%

Rmumu fit.
Chi-squared / dof = 65.2/23, pvalue = small
Run-test p-value = 37%

## CrossSections.png
The 6 figures with superimposed models. Note only Figure 4 and Figure 6 are 
directly fitted.

## CrossSections.log
Logfile from running the code
