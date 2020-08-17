tcpl-py
==============================

Simple example of fitting concentration-response data from high-throughput screening (HTS) using the [R/tcpl](https://cran.r-project.org/web/packages/tcpl/index.html) package in Python3 using rpy2. This assumes you've successfully installed [rpy2](https://pypi.org/project/rpy2/) and other requirements.

Read about the tcplFit approach in [The ToxCast Analysis Pipeline : An R Package for Processing and Modeling Chemical Screening Data](https://academic.oup.com/bioinformatics/article/33/4/618/2617576)

The wrapper is in src/tcpl/fit/tcplfit.py and an example of usage is in notebooks/dev/002-tcpl-v1.ipynb.


Project Organization
------------

    ├── LICENSE
    ├── README.md          <- The top-level README for developers using this project.
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │   └── dev
    │         └──002-tcpl-v1.ipynb  
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment
    ├── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
    └── src                <- Source code for use in this project.
        └─ tcpl
           └─ fit
            └─ tcplfit.py    <- The R/tcpl wrapper


--------

<p><small>Using <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/"> #cookiecutterdatascience</small></a></p>
