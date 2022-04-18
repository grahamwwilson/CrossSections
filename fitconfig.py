# fitconfig.py
#
# Specific fit configurations such as FitBW and FitPoly
# Note that this contains also the initial values for the parameters and 
# optionally whether some parameters are fixed.
#
import mymodels
from iminuit import Minuit
import MyLeastSquares
from fitting import FixParameters, GenericFit

def FitBW(ECM, Ratio, dRatio):

    print(" ")
    print(" ")
    print("---------------------------")
    print("Running fitconfig.FitBW fit")
    print("---------------------------")
    print(" ") 
    
    lsq = MyLeastSquares.LsqDriver(mymodels.rbw, ECM, Ratio, dRatio)
    m = Minuit(lsq, a0 = 3.0, a1 = 0.0, a2 = 0.0, bw=100.0, mz = 91.1876, gz = 2.4952 )
    FixParameters(m, ["mz", "gz"])
    GenericFit(m, lsq, len(ECM))
    
    y_model = mymodels.rbw(ECM, *m.values)
    
    return y_model
    
def FitPoly(x, y, dy):

    print(" ")
    print(" ")
    print("-----------------------------")
    print("Running fitconfig.FitPoly fit")
    print("-----------------------------")
    print(" ")
    
    lsq = MyLeastSquares.LsqDriver(mymodels.polynomial, x, y, dy)
    m = Minuit(lsq, a0 = 0.2, a1 = 0.0, a2 = 0.0, a3=0.0, a4=0.0, a5=0.0, a6=0.0, a7=0.0, a8=0.0)
    FixParameters(m, ["a2", "a4", "a5", "a6", "a7", "a8"])    
    GenericFit(m, lsq, len(x))
    
    y_model = mymodels.polynomial(x, *m.values)
    
    return y_model
