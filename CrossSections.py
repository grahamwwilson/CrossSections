# Plot dimuon cross-sections calculated with Whizard 3.0.3
# Simulations include ISR from both beams and were carried out using a scan of 
# center-of-mass energy from sqrt(s) = 125.0 to 255.0 GeV in 5 GeV steps.
from matplotlib import pyplot as plt
import numpy as np
import pylab
import plotfn                    # Various customized plotting functions
import fitting                   # Put all the fit details in here

ECM = pylab.genfromtxt('Data/ecm.dat',usecols=(0),unpack=True)
xsLR, dxsLR, scaleLR = pylab.genfromtxt('Data/mumu_LR.dat',usecols=(2,3,7),unpack=True) # e- (L) e+ (R) cross-sections in fb
xsRL, dxsRL, scaleRL = pylab.genfromtxt('Data/mumu_RL.dat',usecols=(2,3,7),unpack=True) # e- (R) e+ (L) cross-sections in fb

print(ECM)
print(xsLR)
print(dxsLR)
print(xsRL)
print(dxsRL)

# Convert from fb to pb
xsLR  = 0.001*xsLR
dxsLR = 0.001*dxsLR*np.sqrt(scaleLR)
xsRL  = 0.001*xsRL
dxsRL = 0.001*dxsRL*np.sqrt(scaleRL)

ALR = (xsLR - xsRL)/(xsLR + xsRL)     # Left-right asymmetry
xsU = 0.25*(xsLR + xsRL)
dxsU = 0.25*np.sqrt(dxsLR**2 + dxsRL**2)
# Writing A = (x-y)/(x+y) and doing error propagation
dAdx =  2.0*xsRL/(xsLR + xsRL)**2     # dA/dx
dAdy = -2.0*xsLR/(xsLR + xsRL)**2     # dA/dy
dA = np.sqrt( (dAdx*dxsLR)**2 + (dAdy*dxsRL)**2)

# Let's also define units of R  (86.8 nb/ (ECM in GeV)**2)
QED = 86.8e3/(ECM**2)
print(QED)
dQED=0.0

Ratio = xsU/QED
dRatio = dxsU/QED

plotfn.PlotCustomize()                # Customize fonts and font sizes etc

plotfn.PlotDataWithGrid(1, ECM, xsLR, dxsLR, r'$e^{-}_{L} e^{+}_{R} \to \mu^{-} \mu^{+}$', 'Center-of-Mass Energy (GeV)', 'Cross-Section (pb)')
plotfn.PlotDataWithGrid(2, ECM, xsRL, dxsRL, r'$e^{-}_{R} e^{+}_{L} \to \mu^{-} \mu^{+}$', 'Center-of-Mass Energy (GeV)', 'Cross-Section (pb)')
plotfn.PlotDataWithGrid(3, ECM, xsU, dxsU, r'$e^{-} e^{+}\to \mu^{-} \mu^{+}$', 'Center-of-Mass Energy (GeV)', 'Cross-Section (pb)')
plotfn.PlotDataWithGrid(4, ECM, ALR, dA, r'$e^{-} e^{+}\to \mu^{-} \mu^{+}$', 'Center-of-Mass Energy (GeV)', 'Left-Right Asymmetry')
plotfn.PlotDataWithGrid(5, ECM, QED, dQED, r'$e^{-} e^{+}\to \mu^{-} \mu^{+}$', 'Center-of-Mass Energy (GeV)', 'QED-only Born Cross-Section (pb)')
plotfn.PlotDataWithGrid(6, ECM, Ratio, dRatio, r'$e^{-} e^{+}\to \mu^{-} \mu^{+}$', 'Center-of-Mass Energy (GeV)', 'Cross-Section Ratio to QED-only')
   
# Fit the ALR values vs sqrt(s)
ALR_model   = fitting.FitPoly(ECM, ALR, dA)
# Fit the ratio of unpolarized cross-section to QED point-like cross-section vs sqrt(s)
Ratio_model = fitting.FitBW(ECM, Ratio, dRatio)

plotfn.PlotModel2(4, ECM, ALR_model, r'$e^{-} e^{+}\to \mu^{-} \mu^{+}$', 'Center-of-Mass Energy (GeV)', 'Left-Right Asymmetry', 'magenta')
plotfn.PlotModel2(6, ECM, Ratio_model, r'$e^{-} e^{+}\to \mu^{-} \mu^{+}$', 'Center-of-Mass Energy (GeV)', 'Cross-Section Ratio to QED-only', 'magenta')

# We can also now reconstruct the models for xsLR, xsRL, and xsU, and superimpose those on the respective plots
xsU_model = Ratio_model*QED
xsLR_model = 2.0*xsU_model*(1.0 + ALR_model)
xsRL_model = 2.0*xsU_model*(1.0 - ALR_model)

plotfn.PlotModel2(1, ECM, xsLR_model, r'$e^{-}_{L} e^{+}_{R} \to \mu^{-} \mu^{+}$', 'Center-of-Mass Energy (GeV)', 'Cross-Section (pb)', 'magenta')
plotfn.PlotModel2(2, ECM, xsRL_model, r'$e^{-}_{R} e^{+}_{L} \to \mu^{-} \mu^{+}$', 'Center-of-Mass Energy (GeV)', 'Cross-Section (pb)', 'magenta')
plotfn.PlotModel2(3, ECM, xsU_model, r'$e^{-} e^{+}\to \mu^{-} \mu^{+}$', 'Center-of-Mass Energy (GeV)', 'Cross-Section (pb)', 'magenta')

plt.show()
