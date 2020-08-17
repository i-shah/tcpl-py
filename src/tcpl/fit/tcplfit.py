import matplotlib.text as text
import pandas as pd
import numpy as np
import pylab as pl
import scipy as sp
import re

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

rtcpl = rpackages.importr('tcpl')

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
                if not np.isnan(x):
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

def hillF(x,tp,ga,gw):
    return tp/(1+10**((ga-x)*gw))

def gnlsF(x,tp,ga,gw,lw,la):
    gn = 1/(1+10**((ga-x)*gw))
    ls = 1/(1+10**((x-la)*lw))
    return tp*gn*ls

def calc_Resp(Tcpl):
    BF = Tcpl['best_fit']
    b0 = Tcpl['cr_info']['bmad']
    ch = 1 #BF['ch']
    
    model=None
    if BF['model']=='hill':
        kw = {o:v for o,v in iter(BF.items()) if o in ['tp','ga','gw'] }
        model=hillF
    elif BF['model']=='gnls':
        kw = {o:v for o,v in iter(BF.items()) if o in ['tp','ga','gw','la','lw'] }
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
    
    
