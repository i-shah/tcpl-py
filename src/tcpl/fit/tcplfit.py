import matplotlib.text as text
import pandas as pd
import numpy as np
import pylab as pl
import scipy as sp
import re
import copy

#from IPython.core.debugger import set_trace

import rpy2
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri,FloatVector,ListVector
from rpy2.rinterface import NULL
import scipy.optimize as optz
import seaborn as sns
from scipy import stats
import rpy2.robjects.packages as rpackages
import seaborn as sns
from statsmodels import robust
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib as mpl
from functools import reduce

rtcpl = rpackages.importr('tcplfit2')

def fix_type(x):
    t = type(x)
    if t in [rpy2.rinterface.NARealType,
             rpy2.rinterface.NAIntegerType
            ]:
        return None
    else:
        return x
             
def as_dict(vector):
    """Convert an RPy2 ListVector to a Python dict"""
    result = {}
    for i, name in enumerate(vector.names):
        if isinstance(vector[i], ListVector) :
            result[name] = as_dict(vector[i])
        elif type(vector[i])!=rpy2.rinterface.NULLType:
            if len(vector[i]) == 1:
                x = vector[i][0]
                #if not np.isnan(x):
                result[name]=x
            else:
                result[name] = list(vector[i])
    return result

                
def tcplFit(conc,resp,bmad=False):
    kwargs={'logc':FloatVector(conc),
            'resp':FloatVector(resp),
           }
    if bmad:
        kwargs['bmad']=1.0*bmad
    else:
        kwargs['force.fit']=True
        
    Y0 = rtcpl.tcplFit(**kwargs)
    Y1= pd.Series(as_dict(Y0))
    F0 = []
    for m in ['cnst','hill','gnls']:
        fit=Y1.pop(m)
        K = [i for i in Y1.index if i.startswith(m)]
        Y2= Y1[K]
        Y2.index=[re.sub('%s_?'%m,'',i) for i in K]
        Y2['model']=m
        if fit==fit:
            F0.append(Y2.to_dict())
        Y1= Y1.drop(K)
        
    Y1['bmad']=bmad
    R0 = {}
    R0['fits']=F0
    #set_trace()
    R0['cr_info']=Y1.to_dict()
    R0['cr_data']=dict(conc=list(conc),resp=list(resp))

    # FIgure out best fit
    F1 = pd.DataFrame(F0)
    R0['best_fit']=F1.sort_values('aic').iloc[0].to_dict()
    
    return R0

def tcplFit2(conc,resp,cutoff=False):
    kwargs={'conc':FloatVector(conc),
            'resp':FloatVector(resp),
            'bidirectional':True
           }
    if cutoff:
        kwargs['cutoff']=1.0*cutoff
    else:
        kwargs['force.fit']=True

    try:
      Y0 = as_dict(rtcpl.tcplfit2_core(**kwargs))
    except:
      return {}
    R0 = {}
    Fits = []
    #pd.DataFrame([{k:M[k] for k in ['aic','top','ac50'] if k in M} for M in [Y[m] for m in Y['modelnames']]])
    for m in Y0['modelnames']:
        M = Y0[m]
        F = dict(model=m)
        F.update({k:M[k] for k in ['aic','top','ac50'] if k in M})
        Y0[m]['model']=m
        PR = {}
        for p in M.pop('pars'):
          if p in M: 
            PR[p]=M.pop(p)
        F['params']=PR
        SD = {}
        if 'sds' in M:
          for sd in M.pop('sds'):
            if sd in M: 
              SD[sd]=M.pop(sd)
          F['param_sds']=SD
        if 'modl' in M:
          F['resp_fit'] = M.pop('modl') 
        Fits.append(F)
    Y0.pop('modelnames') 
    R0['cr_data']=dict(conc=list(conc),resp=list(resp))
    R0['fits']=Fits
    
    # FIgure out best fit
    F1 = pd.DataFrame(Fits)
    R0['best_fit']=F1.sort_values('aic').iloc[0].model
    
    return R0

# exp2 f(x) = a*(e^(x/b)- 1)
# exp3 f(x) = a*(e^((x/b)^p)- 1) 
# exp4 f(x) = tp*(1-2^(-x/ga))
# exp5 f(x) = tp*(1-2^(-(x/ga)^p))
# poly1 f(x) = a*x
# poly2 f(x) = a*(x/b + x^2/b^2)
# pow  f(x) = a*x^p

def poly1F(x=1,a=1,**kws):
  return a*x

def poly2F(x=1,a=1,b=1,**kws):
  return a*(x/b + (x**2)/(b**2))

def powF(x=1,a=1,p=1,**kws):
  return a*(x**p)

def hillF(x=1,tp=1,ga=1,p=1,**kws):
  """f(x) = tp/[(1 + (ga/x)^p )]"""
  return tp/(1 + (ga/x)**p)

def gnlsF(x=1,tp=1,ga=1,p=1,la=1,q=1,**kws):
  """gnls f(x) = tp/[(1 + (ga/x)^p )(1 + (x/la)^q )]"""
  return tp/((1 + (ga/x)**p )*(1 + (x/la)**q))

def exp2F(x=1,a=1,b=1,**kws):
  return a*(np.exp(x/b)-1)

def exp3F(x=1,a=1,b=1,p=1,**kws):
  return a*(np.exp((x/b)**p)-1)

def exp4F(x=1,tp=1,ga=1,**kws):
  return tp*(1.0-2**(-x/ga))

def exp5F(x=1,tp=1,ga=1,p=1,**kws):
  """f(x) = tp*(1-2^(-(x/ga)^p))"""
  return tp*(1.0-2**(-(x/ga)**p))


def calc_Resp(Tcpl):
    BF = Tcpl['best_fit']
    b0 = Tcpl['cr_info']['bmad']
    ch = 1 #BF['ch']
    
    model=None
    if BF['model']=='hill':
        kw = {o:v for o,v in iter(BF.items()) if o in ['tp','ga','gw'] }
        model=hillF
    elif BF['model']=='gnls':
        kw = {k:v for o,v in iter(BF.items()) if 
              o in ['tp','ga','gw','la','lw'] }
        model=gnlsF
    elif BF['model']=='cnst':
        kw = {}
        def model(x): 
            y=np.median(Tcpl['cr_data']['resp'])
            return [y]*len(x)
    
    # CR Data
    cr_c = np.array(Tcpl['cr_data']['conc'])
    cr_r = np.array(Tcpl['cr_data']['resp'])    
    c0 = cr_c.min()-0.8

    # Response
    ci,cf = cr_c.min()-1,cr_c.max()+1
    C = np.linspace(ci,cf,50)
    R = model(C,**kw)

    return cr_c,cr_r,model,kw

    
def tcplPlot(Tcpl,ax=None,show_data=True,fnsz=8,r_max=None,
             show_bmds=dict(E=[10,20,50],Z=[1,2,3]),show_legend=True,
             cols=dict(data='green',fit='seagreen',ac50='red',bmad='skyblue')
            ):
    BF = Tcpl['best_fit']
    b0 = Tcpl['cr_info']['bmad']
    ch = 1
    
    if not ax:
        ax = pl.subplot(111)
    
    # CR Data
    cr_c,cr_r,model,kw = calc_Resp(Tcpl)

    c0 = cr_c.min()-0.8

    if show_data:
        ax.scatter(cr_c,cr_r,marker='+',s=30,c=cols['data'],alpha=0.7)

    # Plot model results
    ci,cf = cr_c.min()-1,cr_c.max()+1
    C = np.linspace(ci,cf,50)
    R = model(C,**kw)
    ax.plot(C,R,color=cols['fit'],label='fit',alpha=0.7)

    if not r_max:
        r_max = np.max([max(R),cr_r.max()])

    # SHow the bmads
    ax.hlines(b0*np.array([-2,-2,-1,1,2,3]),ci,cf,
              lw=0.3,color=cols['bmad'],linestyle='-.')
    
    # Show AC50
    if 'ga' in BF:
       ax.vlines(BF['ga'],0,r_max,color='red',lw=1)
       ax.text(c0,r_max*0.9,r"$AC_{50}$=%5.2f$\mu$M" % 10**BF['ga'],
              fontdict=dict(size=fnsz,color=cols['ac50'],alpha=0.6))

    # Find the BMD10

    ax.set_xlim(ci,cf)
    
    #Conc=np.linspace(cr_c.min()-0.5,cr_c.max()+0.5,9)
    Conc=np.array([0.1,0.5,1,5,10,25,50,100,250])
    Clab=['0.1','0.5','1','5','10','25','50','100','250']
    #Clab=['%5.2f'%c for c in Conc]
    
    ax.set_xticks(np.log10(Conc), minor=False)
    xl=ax.set_xticklabels(Clab,rotation=90)
    #for tick in ax.get_xticklines(): 
    #    tick.set_visible(False)
    for tick in ax.get_xticklabels(): 
        tick.set_fontsize(8)

    if show_legend: ax.legend(fontsize='xx-small')
    
    ax.set_ylim(-0.2,r_max*1.1)    
    
    
def calc_BMDs(Fit,BMR=[1,2,3],C=None,
              ret='dict',add_info=False,
              dbg=False):
  """
  Calculate benchmark doses corresponding to bmrs
  """
  ci,cf = C.min(),C.max()

  BMD=[]
  P = copy.deepcopy(Fit['params'])
  
  for bmr in BMR:
    def bmdf(c): 
      P['x']=c
      return eval("{}F(**P)".format(Fit['model'])) - bmr
    
    try:
        soln = optz.fsolve(bmdf,[ci,cf])
        soln = soln[np.logical_and(soln>ci,soln<cf)]
        bmd0 = np.min(soln)
    except:
        bmd0 = None
        if dbg: print("Failed E %0.2f" % bmr)
    else:
        BMD.append(dict(bmr=bmr,bmd=bmd0))
  
  return BMD
  