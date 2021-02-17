import pandas as pd
import numpy as np
import pylab as pl
import scipy as sp
import sys
import rpy2 
import os 
import pdb

from joblib import Parallel, delayed, parallel_backend
from joblib import wrap_non_picklable_objects
import json

from functools import reduce
import random,time 
#from pathos.pools import ProcessPool

from tcpl.fit.crvfit import *

class MultiCurveFit:
    
    def __init__(self,*args,n_jobs=2,**kwargs):
        self._CF_kw = kwargs
        self._CF = CurveFit(**kwargs)
        self._Conc = None
        self._Resp = None
        self._Cutoff=None
        self._Onesd= None

        
    def __call__(self, Conc=[],Resp=[],Cutoff=None, 
                Onesd=1,bmr_magic=1,ret=False,
                n_jobs=10,
                **kwargs):
        """
        Parallel multiple curve-fitting
        
        Parameters:
        -----------
        
        
        Conc = List/pd.Series of concentrations 
        Resp = pd.DataFrame in which rows match concentrations in Conc and columns 
               are L2FC values for individual endpoints
        Cutoff=pd.Series in which index corresponds to the columns in Resp and values
               are the cutoffs for individual endpoints
        Onesd =pd.Series in which index corresponds to the columns in Resp and values
               are one standard deviation for individual endpoints
        n_jobs= number of parallel threads
        """

        self._Conc = Conc = list(Conc)
        self._Resp = Resp
        self._Cutoff=Cutoff
        self._Onesd= Onesd

        def myfunc(*args,**kwargs):
          CF = CurveFit(**self._CF_kw)
          H  = CF(*args,**kwargs)
          return json.dumps(H)

        Hits = Parallel(n_jobs=n_jobs)(delayed(myfunc)(Conc,list(Resp[assay]),
                                                   cutoff=float(Cutoff[assay]),
                                                   assay=assay,
                                                   onesd=float(Onesd[assay]),
                                                   summary=True)
                                       for assay in Resp.columns)       
            
        self._Hits = [json.loads(i) for i in Hits]
        if ret: return self._Hits
            
    def fits(self, Conc=[],Resp=[],Cutoff=None, 
                 Onesd=1,bmr_magic=1,ret=False,**kwargs):
        """
        Serial multiple curve-fitting
        Conc = List/pd.Series of concentrations 
        Resp = pd.DataFrame in which rows match concentrations in Conc and columns 
               are L2FC values for individual endpoints
        Cutoff=pd.Series in which index corresponds to the columns in Resp and values
               are the cutoffs for individual endpoints
        Onesd =pd.Series in which index corresponds to the columns in Resp and values
               are one standard deviation for individual endpoints 
        """

        self._Conc = Conc
        self._Resp = Resp
        self._Cutoff=Cutoff
        self._Onesd= Onesd

        CF = self._CF

        Hits = []

        for assay in Resp.columns:
            #pdb.set_trace()
            #print(assay)
            CF(list(Conc),list(Resp[assay]),
               cutoff=float(Cutoff[assay]),onesd=float(Onesd[assay]))
            F = CF.get_summary()
            F.update(assay=assay)
            Hits.append(F)
        
        self._Hits = Hits

        if ret: return Hits


#     def call_new(self, Conc=[],Resp=[],Cutoff=None, 
#                 Onesd=1,bmr_magic=1,ret=False,
#                 n_jobs=10,
#                 **kwargs):
#         """
#         Parallel multiple curve-fitting
        
        
#         Conc = List/pd.Series of concentrations 
#         Resp = pd.DataFrame in which rows match concentrations in Conc and columns 
#                are L2FC values for individual endpoints
#         Cutoff=pd.Series in which index corresponds to the columns in Resp and values
#                are the cutoffs for individual endpoints
#         Onesd =pd.Series in which index corresponds to the columns in Resp and values
#                are one standard deviation for individual endpoints
#         n_jobs= number of parallel threads
#         """

#         self._Conc = Conc = list(Conc)
#         self._Resp = Resp
#         self._Cutoff=Cutoff
#         self._Onesd= Onesd

#         CF = self._CF
#         PP = ProcessPool(n_jobs)
        
#         Hits = PP.map(lambda assay: CF(Conc,Resp[assay],cutoff=float(Cutoff[assay]),
#                                        assay=assay,onesd=float(Onesd[assay]),
#                                        summary=True),
#                                  Resp.columns)
        
            
#         self._Hits = Hits

#         if ret: return Hits
        