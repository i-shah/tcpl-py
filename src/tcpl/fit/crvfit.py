import pandas as pd
import numpy as np
import pylab as pl
import scipy as sp
import re

import rpy2
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri,FloatVector,ListVector
from rpy2.rinterface import NULL
import scipy.optimize as optz
import seaborn as sns
from scipy import stats
import rpy2.robjects.packages as rpackages
import rpy2.robjects as robjects
import seaborn as sns
from statsmodels import robust
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib as mpl
from functools import reduce


rtcpl = rpackages.importr('tcplfit2')

from .rtcplhit import tcplhit2_in_r

tcplhit2 = robjects.r(tcplhit2_in_r)

class CurveFit:
    
  def __init__(self,conc=[],resp=[],cutoff=None,r_fill=0,hit_call=True,bmr_magic=1.349):
    self.C = conc
    self.R = resp
    self.r0= cutoff
    self.bmr_magic=bmr_magic
    self.r_fill = r_fill
    self.hit_call=hit_call
    self.RFits= None
    self.Fits = None
    self.Hit  = None
    
  def __call__(self, conc=[],resp=[],cutoff=None,assay=None, 
               onesd=1,bmr_magic=None,summary=False,**kwargs):
    bmr_magic = bmr_magic if bmr_magic else self.bmr_magic
    self.fit(conc=conc,resp=resp,cutoff=cutoff)
    if self.hit_call: 
      self.hit(onesd=onesd,bmr_magic=bmr_magic,cutoff=cutoff)
    if assay: self.assay = assay
    if summary:
      return self.get_summary(BMR=[1])
    
    
  def get_summary(self,BMR=[-3,-2,-1,1,2,3]):
    BF  = self.get_best_fit()
    Summary = {k:v for k,v in BF.items() if k in['model','aic','top','ac50']}
    #BMD = self.calc_bmds(BMR,model=None)
    #if len(BMD)>0:
    #  Summary.update(BMD[0])
    if self.hit_call: Summary.update(self.Hit.to_dict())
    if self.assay: Summary.update(dict(assay=self.assay))
      
    return Summary
    
  def fit(self,conc=[],resp=[],cutoff=None):
    self.C = conc if len(conc)>0 else self.C
    self.R = resp if len(resp)>0 else self.R
    self.r0= cutoff if cutoff else self.r0
    
    kwargs={'conc':FloatVector(self.C),
            'resp':FloatVector(self.R),
            'bidirectional':True
           }
    if cutoff:
        kwargs['cutoff']=1.0*self.r0
    else:
        kwargs['force.fit']=True
        
    try:
      RFits = rtcpl.tcplfit2_core(**kwargs)
      Y0 = self.as_dict(RFits)
    except:
      return {'failed':1}

    Fits = []
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
    
    self.CR_data=dict(conc=list(conc),resp=list(resp))
    self.Fits  = Fits
    self.RFits = RFits 

    # FIgure out best fit
    F1 = pd.DataFrame(Fits)
    self.best=F1.sort_values('aic').iloc[0].model       
    
  def as_dict(self,vector):
    """Convert an RPy2 ListVector to a Python dict"""
    result = {}
    for i, name in enumerate(vector.names):
        if isinstance(vector[i], ListVector) :
            result[name] = self.as_dict(vector[i])
        elif type(vector[i])!=rpy2.rinterface.NULLType:
            if len(vector[i]) == 1:
                x = vector[i][0]
                #if not np.isnan(x):
                result[name]=x
            else:
                result[name] = list(vector[i])
    return result

  def get_fits(self):
    return self.Fits
  
  def get_model_fit(self,model):
    return next(i for i in self.Fits if i['model']==model)    
  
  def get_best_fit(self):
    return next(i for i in self.Fits if i['model']==self.best)

  def plot_fits(self,C=10**np.linspace(-1,2),interp=False,pli=pl,
                best_only=False,lab=None):
    pli.scatter(self.C,self.R)
    CRv = CRCurve()
    best= self.best
  
    for Fit in self.Fits:
      if not 'resp_fit' in Fit: continue
      if best_only and Fit['model']!=best: continue
        
      if interp:
        R = CRv.curve(Fit['model'],**Fit['params'])(C)
      else:
        R = Fit['resp_fit']
        C = self.C
      #label = Fit['model'] if Fit['model']!=best and show_best else best+' *'
      label = Fit['model'] if not lab else lab
      pli.plot(C,R,label=label)
      #pli.hlines([self.r0,-1*self.r0],self.C.min(),self.C.max(),
      #          linestyle='--',colors='grey')
      pli.axhline(self.r0,lw=0.5,color='grey')
      pli.axhline(-1*self.r0,lw=0.5,color='grey')
      
    if self.hit_call:
      pli.axhline(self.Hit.bmr,lw=0.8,ls='-',color='green',label='bmr')
      pli.axvspan(self.Hit.bmdl,self.Hit.bmdu,color='grey',alpha=0.2,label='bmd[ul]')
      pli.axvline(self.Hit.bmd,lw=0.8,ls='--',color='green',label='bmc')
      pli.axvline(self.Hit.ac50,lw=1,ls='-',color='red',label='ac50')
    pli.legend(fontsize='small')
    
  def calc_bmds(self,BMR=[1,2,3],
                model='hill',
                dbg=False):
    """
    Calculate benchmark doses corresponding to bmrs 
    """
    ci,cf = self.C.min(),self.C.max()
    if ci==0: ci += 1e-5
    BMD=[]
    F = self.get_best_fit() if not model else self.get_model_fit(model)
    P = F['params']
    CV = CRCurve()

    for bmr in BMR:
      BMRF = CV.bmrf(F['model'],bmr=bmr,**P)
      try:
          soln = optz.brentq(BMRF,ci,cf*10)
          bmd0 = np.min(soln)
      except:
          bmd0 = None
          if dbg: print("Failed E %0.2f" % bmr)
      else:
          BMD.append(dict(bmr=bmr,bmd=bmd0))

    return BMD
  
  def hit(self,onesd=1,bmr_magic=1.349,cutoff=None):
    
    kwargs={'params':self.RFits,
            'conc':FloatVector(self.C),
            'resp':FloatVector(self.R),
            'cutoff':cutoff or self.r0,
            'onesd':onesd,
            'bmr_magic':bmr_magic
            }
    RHit = tcplhit2(**kwargs)
    self.Hit = pd.Series(self.as_dict(RHit)).dropna()
    
  
class CRCurve:
  
  def __init__(self,name=""):
    self._name = name
    
  def calc_resp(self,x,params,**kws):
    pass
  
  def curve(self,model,**params):
     return lambda x: getattr(self,"{}F".format(model))(x,**params)
  
  def bmrf(self,model,bmr=1,**params):
    F = self.curve(model,**params)
    return lambda c: np.abs(F(c))-bmr
  
  def poly1F(self,x=1,a=1,**kws):
    return a*x

  def poly2F(self,x=1,a=1,b=1,**kws):
    return a*(x/b + (x**2)/(b**2))

  def powF(self,x=1,a=1,p=1,**kws):
    return a*(x**p)

  def hillF(self,x=1,tp=1,ga=1,p=1,**kws):
    """f(x) = tp/[(1 + (ga/x)^p )]"""
    return tp/(1 + (ga/x)**p)

  def gnlsF(self,x=1,tp=1,ga=1,p=1,la=1,q=1,**kws):
    """gnls f(x) = tp/[(1 + (ga/x)^p )(1 + (x/la)^q )]"""
    return tp/((1 + (ga/x)**p )*(1 + (x/la)**q))

  def exp2F(self,x=1,a=1,b=1,**kws):
    return a*(np.exp(x/b)-1)

  def exp3F(self,x=1,a=1,b=1,p=1,**kws):
    return a*(np.exp((x/b)**p)-1)

  def exp4F(self,x=1,tp=1,ga=1,**kws):
    return tp*(1.0-2**(-x/ga))

  def exp5F(self,x=1,tp=1,ga=1,p=1,**kws):
    """f(x) = tp*(1-2^(-(x/ga)^p))"""
    return tp*(1.0-2**(-(x/ga)**p))
  
  
  