# transferfactor
Python package for performing background estimation using 2D transfer factor (TF) method

## Quick start

On lxplus, do:
```
$ cd your/work/directory
$ git clone git@github.com:asogaard/transferfactor.git
$ cd transferfactor
$ git clone git@github.com:asogaard/rootplotting.git
$ source pythonenv.sh
$ python examples/closure.py -h
$ python examples/closure.py --mass 80 --show
```

From there on, try the different macros in [examples](examples), which all use a similar syntax. Change variables, binning, and ROOT TTree names in [transferfactor/config.py](transferfactor/config.py). Use your own files by setting the appropriate base path in [transferfactor/config.py](transferfactor/config.py), and file names in the macro in question.