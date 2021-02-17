tcpl-py
==============================

Simple example of fitting concentration-response data from high-throughput screening (HTS) using the [R/tcpl](https://cran.r-project.org/web/packages/tcpl/index.html) package in Python3 using rpy2. This assumes you've successfully installed [rpy2](https://pypi.org/project/rpy2/) and other R requirements.

Read about the tcplFit approach in [The ToxCast Analysis Pipeline : An R Package for Processing and Modeling Chemical Screening Data](https://academic.oup.com/bioinformatics/article/33/4/618/2617576)


# Project Organization

    tcpl-py
    ├── LICENSE
    ├── notebooks
    │   └── dev
    │       └── 004-tcpl-v2.ipynb
    ├── README.md
    ├── requirements.txt
    ├── setup.py
    └── src
        └── tcpl
            └── fit
                ├── crvfit.py
                ├── multicf.py
                └── rtcplhit.py

# Installation

## environment 

Create tcpl-py/notebooks/.env file with the following:-

TOP=path_to_top/tcpl-py
LIB=path_to_top/tcpl-py/src
DAT=path_to_top/tcpl-py/data
FIG=path_to_top/tcpl-py/figs 


<p><small>Using <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/"> #cookiecutterdatascience</small></a></p>
