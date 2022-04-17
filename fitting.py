# fitting.py
# The generic code is in GenericFit 
# and the fit configurations are in the specific cases, like FitBW
import mymodels
from scipy import stats
from iminuit import Minuit
import MyLeastSquares

def GenericFit(m, lsq, ndata):
# 1. Global minuit settings. Can be defined after m is constructed.
    m.print_level = 2
    m.errordef = 1.0
# 2. Do the fit
    m.migrad(ncall=100000)

# 3. Check how many parameters are fitted (ie not fixed) in the fit
    nfitted = 0
    nfixed = 0
    for p in m.parameters:
        if m.fixed[p] == True:
            print('Parameter ',p,' is fixed ')
            nfixed +=1
        else:
            print('Parameter ',p,' is fitted ')
            nfitted +=1
            
# 4. Do error analysis
    if nfitted>0:  # Before proceeding with hesse and error analysis, make sure 
                   #that there are free parameters to fit
        m.hesse(ncall=100000)
   
# 5. Hypothesis test results based solely on chi-squared
    print(' ')
    print('A) Hypothesis test: IS THIS A REASONABLE FIT? ')
    print('chi-squared                          : ', m.fval)
    print('number of data points                : ', ndata)
    print('total number of fit parameters       : ', nfitted+nfixed, 
          ' of which ',nfitted,' are/is fitted and ',nfixed,'are/is fixed')
    print('number of ACTUAL fitted parameters   : ', nfitted)
    ndof = ndata - nfitted
    print('number of degrees of freedom (d.o.f.): ', ndof)
    print('chi-squared ratio is ',m.fval/ndof)

    pvaluepercent = 100.0*(1.0 - stats.chi2.cdf(m.fval, ndof ))
    print(' ')
    print('Observed chi-squared p-value of',pvaluepercent,'%')
    
# 6. Parameter Estimation results including covariance matrix and correlations
    print(' ')
    print('B) Parameter Estimation: Given a reasonable fit, what are the model parameters and uncertainties?')
    print('Is the covariance matrix accurate? ',m.accurate)
    print('Fitted parameters:')
    for p in m.parameters:
        if m.fixed[p] == False:
            print("{} = {} +- {} ".format(p,m.values[p], m.errors[p] ))
        else:
            print("{} = {}  (Fixed = {})".format(p,m.values[p], m.fixed[p] ))
    print(' ')

    if nfitted > 0:  # Likewise only do this if there are free parameters
        print('Correlation Coefficient Matrix')
        print(m.covariance.correlation())
        print(' ')
        
# 7. Use fitted parameter values to evaluate Run Test statistic
    rpval = lsq.runspvalue(*m.values)  # this is a method in MyLeastSquares.py
    print('Observed run test p-value (%) = ',rpval)
    return    

#FIXME DRY principle says we should have an interface to fix a list of parameters
# rather than repeating this.

def FitBW(ECM, Ratio, dRatio):
    print(" ")
    print(" ")
    print("---------------------------")
    print("Running fitting.FitBW fit")
    print("---------------------------")
    print(" ")    
    lsq = MyLeastSquares.LsqDriver(mymodels.rbw, ECM, Ratio, dRatio)
    m = Minuit(lsq, a0 = 3.0, a1 = 0.0, a2 = 0.0, bw=100.0, mz = 91.1876, gz = 2.4952 )
    m.fixed["gz"] = True
    m.fixed["mz"] = True
    ndata = len(ECM)
    GenericFit(m, lsq, ndata)
    y_model = mymodels.rbw(ECM, *m.values)
    return y_model
    
def FitPoly(x, y, dy):
    print(" ")
    print(" ")
    print("---------------------------")
    print("Running fitting.FitPoly fit")
    print("---------------------------")
    print(" ")
    lsq = MyLeastSquares.LsqDriver(mymodels.polynomial, x, y, dy)
    m = Minuit(lsq, a0 = 0.2, a1 = 0.0, a2 = 0.0, a3=0.0, a4=0.0, a5=0.0, a6=0.0, a7=0.0, a8=0.0)
    m.fixed["a2"] = True
    m.fixed["a4"] = True
    m.fixed["a5"] = True
    m.fixed["a6"] = True
    m.fixed["a7"] = True
    m.fixed["a8"] = True
    ndata = len(x)
    GenericFit(m, lsq, ndata)
    y_model = mymodels.polynomial(x, *m.values)
    return y_model    
