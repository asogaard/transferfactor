# -*- coding: utf-8 -*-

""" Class for computing transfer factor for dijet + ISR analysis

@file:   calculator.py
@date:   1 May 2017
@author: Andreas Søgaard 
@email:  andreas.sogaard@cern.ch
"""

# Basic import(s)
import sys
import itertools

# Scientific import(s)
try:

    # ROOT
    import ROOT

    # Numpy
    import numpy as np
    from root_numpy import *

    # Matplotlib
    import matplotlib.pyplot as plt

    # Scikit-learn
    from sklearn.gaussian_process import GaussianProcess

except ImportError:
    print "WARNING: One or more scientific python packages were not found. If you're in lxplus, try running:"
    print "         $ source /cvmfs/sft.cern.ch/lcg/views/LCG_88/x86_64-slc6-gcc49-opt/setup.sh"
    sys.exit()
    pass

# Local import(s)
try:
    from rootplotting import ap
except ImportError:
    print "WARNING: Needs 'rootplotting' package to produce plots."
    sys.exit()
    pass
from utils import *


class calculator (object):
    """docstring for calculator"""

    def __init__(self, data, config, subtract=None, verbose=True):
        super(calculator, self).__init__()
        self._data     = data
        self._subtract = subtract
        self._config   = config
        self._verbose  = verbose

        self._mass   = None
        self._window = None

        self._clf = None
        self._toyclfs = list()
        self._partialbins = True
        self._emptybins   = False
        self._fitted     = False
        self._fullfitted = False
        self._shift = None
        return



    # Properties
    # --------------------------------------------------------------------------

    @property
    def mass (self):
        """ The 'mass' property """
        return self._mass

    @mass.setter
    def mass (self, value):
        # ...
        self._mass = value
        return
        
    @mass.deleter  
    def mass (self):
        del self._mass
        return


    @property
    def window (self):
        """ The 'window' property """
        return self._window

    @window.setter
    def window (self, value):
        # ...
        self._window = value
        return

    @window.deleter
    def window (self):
        del self._window
        return


    @property
    def partialbins (self):
        """ The 'partialbins' property """
        return self._partialbins

    @partialbins.setter
    def partialbins (self, value):
        # ...
        self._partialbins = value
        return

    @partialbins.deleter
    def partialbins (self):
        del self._partialbins
        return


    @property
    def emptybins (self):
        """ The 'emptybins' property """
        return self._emptybins

    @emptybins.setter
    def emptybins (self, value):
        # ...
        self._emptybins = value
        return

    @emptybins.deleter
    def emptybins (self):
        del self._emptybins
        return


    def theta (self):
        """ The classifier 'theta' accessor """
        return self._clf.theta_.flatten()




    # High-level method(s)
    # --------------------------------------------------------------------------

    def fit (self, theta=None):
        """ ... """

        # Check(s)
        if self._verbose: print "  Performing checks"

        # -- Mass window
        assert (self._mass is None) == (self._window is None), "Both mass and window width must either be specified, or not."
        assert (self._window is None) or (self._window < 1.), "Relative window width must be less than 100%."

        # -- self._config
        assert self._config, "Need configuration to fit."
        assert type(self._config) == dict, "Configuration must be of type dict."
        for field in ['params', 'pass', 'axes']:
            assert field in self._config, "Configuration must contain field '{:s}'.".format(field)
            pass
        assert len(self._config['params']) == len(self._config['axes']), "Number of dimensions in parametrisation ({:d}) and number of axes definitions ({:d}) must agree".format(len(self._config['params']), len(self._config['axes']))

        if self._fitted:
            print "  Overwriting previous fit"
            pass


        # Preparing transfer factor for fitting
        # ----------------------------------------------------------------------

        # Get pass/fail masks
        if self._verbose: print "  Getting pass- and fail masks for input data" 
        msk_data_pass = self._config['pass'](self._data)
        msk_data_fail = ~msk_data_pass

        if self._subtract is not None:
            msk_sub_pass = self._config['pass'](self._subtract)
            msk_sub_fail = ~msk_sub_pass
            pass

        # Get pass/fail histograms
        if self._verbose: print "  Generating pass- and fail histograms"
        # -- CR (fit)
        h_data_CR_pass     = get_histogram(self._data,     self._config['params'], self._config['axes'], mask=msk_data_pass &  in_CR(self._data,     self._mass, self._window))
        h_data_CR_fail     = get_histogram(self._data,     self._config['params'], self._config['axes'], mask=msk_data_fail &  in_CR(self._data,     self._mass, self._window))
        if self._subtract is not None:
            h_sub_CR_pass  = get_histogram(self._subtract, self._config['params'], self._config['axes'], mask=msk_sub_pass  &  in_CR(self._subtract, self._mass, self._window))
            h_sub_CR_fail  = get_histogram(self._subtract, self._config['params'], self._config['axes'], mask=msk_sub_fail  &  in_CR(self._subtract, self._mass, self._window))
            pass

        # -- SR (interpolation)
        h_data_SR_pass     = get_histogram(self._data,     self._config['params'], self._config['axes'], mask=msk_data_pass & ~in_CR(self._data,     self._mass, self._window))
        h_data_SR_fail     = get_histogram(self._data,     self._config['params'], self._config['axes'], mask=msk_data_fail & ~in_CR(self._data,     self._mass, self._window))
        if self._subtract is not None:
            h_sub_SR_pass  = get_histogram(self._subtract, self._config['params'], self._config['axes'], mask=msk_sub_pass  & ~in_CR(self._subtract, self._mass, self._window))
            h_sub_SR_fail  = get_histogram(self._subtract, self._config['params'], self._config['axes'], mask=msk_sub_fail  & ~in_CR(self._subtract, self._mass, self._window))
            pass

        # Subtract MC component from data (opt.)
        if self._subtract is not None:
            if self._verbose: print "  Subtracting component from pass- and fail histograms"
            if self._verbose: print "   -- (CR) Pass/fail component %.1f/%.1f from %.1f/%.1f" % (h_sub_CR_pass.Integral(), h_sub_CR_fail.Integral(), h_data_CR_pass.Integral(), h_data_CR_fail.Integral())
            h_data_CR_pass.Add(h_sub_CR_pass, -1)
            h_data_CR_fail.Add(h_sub_CR_fail, -1)
            if self._verbose: print "   -- (SR) Pass/fail component %.1f/%.1f from %.1f/%.1f" % (h_sub_SR_pass.Integral(), h_sub_SR_fail.Integral(), h_data_SR_pass.Integral(), h_data_SR_fail.Integral())
            h_data_SR_pass.Add(h_sub_SR_pass, -1)
            h_data_SR_fail.Add(h_sub_SR_fail, -1)
            pass

        # Compute ratio
        if self._verbose: print "  Getting TF ratio in SR and CR"
        self._TF_CR_mean, self._TF_CR_err = get_ratio_numpy(h_data_CR_pass, h_data_CR_fail)
        self._TF_SR_mean, self._TF_SR_err = get_ratio_numpy(h_data_SR_pass, h_data_SR_fail)

        # Remove partially filled bins (optional)
        if not self._partialbins:
            if self._verbose: print "  Removing partially filled bins"
            h_TF_CR = h_data_CR_pass.Clone('h_TF_CR')
            h_TF_CR.Divide(h_data_CR_fail)

            h_TF_CR = fixHist(h_TF_CR, self._mass, self._window)
            self._TF_CR_mean = hist2array(h_TF_CR).T
            self._TF_CR_err[self._TF_CR_mean == 0] = 0
            pass


        # Fitting transfer factor
        # ----------------------------------------------------------------------

        # Mesh
        if self._verbose: print "  Setting up grid for fitting"
        X1, X2 = np.meshgrid(*self._config['centres']) # *self._config['axes']
        X = np.vstack((X1.ravel(), X2.ravel())).T
        
        # Mean- and error arrays
        if self._verbose: print "  Setting up mean- and error arrays"
        y = self._TF_CR_mean.ravel()
        s = self._TF_CR_err .ravel()

        # Fit weights
        if self._verbose: print "  Computing fit weights"
        msk_fit = (y > 0) # non-zero content
        #w = np.ones_like(s) * eps
        #w[msk_fit] = np.power(s[msk_fit], -1)

        # Remove empty bins (opt.)
        if self._emptybins:
            X_fit = X[:,:]
            y_fit = y[:]
            #w_fit = w[:]
            s_fit = s[:]
        else:
            if self._verbose: print "  Removing empty bins"
            X_fit = X[msk_fit,:]
            y_fit = y[msk_fit]
            #w_fit = w[msk_fit]
            s_fit = s[msk_fit]
            pass

        # Perform fit
        if self._verbose: print "  Performing fit to TF profile"
        #s_fit = np.power(w_fit, -1)
        nugget = np.square(s_fit/(y_fit + eps)).ravel()
        if theta is None:
            # Using ML-optimised theta
            self._clf = GaussianProcess(theta0=[1E-01, 1E-01], thetaL=[1E-03, 1E-03], thetaU=[1E+01, 1E+01], nugget=nugget)
        else:
            # Using manually set theta
            self._clf = GaussianProcess(theta0=theta, nugget=nugget)
            pass
        self._clf.fit(X_fit, y_fit)

        print "  Best value(s) of theta found:", self._clf.theta_

        # ...

        self._fitted = True
        return


    def toysfit (self, N=1, theta=None):
        """ ... """

        # Check(s)
        if self._verbose: print "  Performing checks"

        # -- Mass window
        assert (self._mass is None) == (self._window is None), "Both mass and window width must either be specified, or not."
        assert (self._window is None) or (self._window < 1.), "Relative window width must be less than 100%."

        # -- self._config
        assert self._config, "Need configuration to fit."
        assert type(self._config) == dict, "Configuration must be of type dict."
        for field in ['params', 'pass', 'axes']:
            assert field in self._config, "Configuration must contain field '{:s}'.".format(field)
            pass
        assert len(self._config['params']) == len(self._config['axes']), "Number of dimensions in parametrisation ({:d}) and number of axes definitions ({:d}) must agree".format(len(self._config['params']), len(self._config['axes']))

        if self._fitted:
            print "  Overwriting previous fit"
            pass


        # Preparing transfer factor for fitting
        # ----------------------------------------------------------------------

        # Get pass/fail masks
        if self._verbose: print "  Getting pass- and fail masks for input data" 
        msk_data_pass = self._config['pass'](self._data)
        msk_data_fail = ~msk_data_pass

        # Get pass/fail histograms
        if self._verbose: print "  Generating pass- and fail histograms"
        # -- CR (fit)
        h_data_CR_pass     = get_histogram(self._data,     self._config['params'], self._config['axes'], mask=msk_data_pass &  in_CR(self._data,     self._mass, self._window))
        h_data_CR_fail     = get_histogram(self._data,     self._config['params'], self._config['axes'], mask=msk_data_fail &  in_CR(self._data,     self._mass, self._window))

        # -- SR (interpolation)
        h_data_SR_pass     = get_histogram(self._data,     self._config['params'], self._config['axes'], mask=msk_data_pass & ~in_CR(self._data,     self._mass, self._window))
        h_data_SR_fail     = get_histogram(self._data,     self._config['params'], self._config['axes'], mask=msk_data_fail & ~in_CR(self._data,     self._mass, self._window))

        # Set up TF grid
        if self._verbose: print "  Setting up grid for fitting"
        X1, X2 = np.meshgrid(*self._config['centres'])
        X = np.vstack((X1.ravel(), X2.ravel())).T


        # Fitting transfer factor
        # ----------------------------------------------------------------------
        thetas = list()
        for idx in range(N):
            print "== Fitting experiment %*d/%d:" % (int(np.ceil(np.log10(N + 1))), idx + 1, N),

            # Generate toy pass/fail histograms
            h_toys_CR_pass = h_data_CR_pass.Clone('h_toys_CR_pass')
            h_toys_CR_fail = h_data_CR_fail.Clone('h_toys_CR_fail')
            h_toys_SR_pass = h_data_SR_pass.Clone('h_toys_SR_pass')
            h_toys_SR_fail = h_data_SR_fail.Clone('h_toys_SR_fail')

            h_toys_CR_pass.Reset()
            h_toys_CR_fail.Reset()
            h_toys_SR_pass.Reset()
            h_toys_SR_fail.Reset()
            
            for binx, biny in itertools.product(1 + np.arange(h_toys_CR_pass.GetXaxis().GetNbins()),
                                                1 + np.arange(h_toys_CR_pass.GetYaxis().GetNbins())):

                # -- Get x- and y-coordinates
                x = h_toys_CR_pass.GetXaxis().GetBinCenter(binx)
                y = h_toys_CR_pass.GetYaxis().GetBinCenter(biny)

                # -- Fill all toys histograms
                content = np.random.poisson(h_data_CR_pass.GetBinContent(binx, biny))
                for _ in range(content):
                    h_toys_CR_pass.Fill(x, y)
                    pass

                content = np.random.poisson(h_data_CR_fail.GetBinContent(binx, biny))
                for _ in range(content):
                    h_toys_CR_fail.Fill(x, y)
                    pass

                content = np.random.poisson(h_data_SR_pass.GetBinContent(binx, biny))
                for _ in range(content):
                    h_toys_SR_pass.Fill(x, y)
                    pass

                content = np.random.poisson(h_data_SR_fail.GetBinContent(binx, biny))
                for _ in range(content):
                    h_toys_SR_fail.Fill(x, y)
                    pass
                pass

            # ...

            # Compute ratio
            if self._verbose: print "  Getting TF ratio in SR and CR"
            self._TF_CR_mean, self._TF_CR_err = get_ratio_numpy(h_toys_CR_pass, h_toys_CR_fail)
            self._TF_SR_mean, self._TF_SR_err = get_ratio_numpy(h_toys_SR_pass, h_toys_SR_fail)

            # Mean- and error arrays
            if self._verbose: print "  Setting up mean- and error arrays"
            y = self._TF_CR_mean.ravel()
            s = self._TF_CR_err .ravel()
                                   
            # Remove empty bins        
            if self._verbose: print "  Removing empty bins"
            msk_fit = (y > 0) # non-zero content
            X_fit = X[msk_fit,:]
            y_fit = y[msk_fit]
            s_fit = s[msk_fit]

            # Perform fit
            if self._verbose: print "  Performing fit to TF profile"
            nugget = np.square(s_fit/(y_fit + eps)).ravel()
            
            if theta is None:
                # Using ML-optimised theta
                if self._verbose: print "  Using ML-optimised theta"
                self._clf = GaussianProcess(theta0=[1E-01, 1E-01], thetaL=[1E-03, 1E-03], thetaU=[1E+01, 1E+01], nugget=nugget)
            else:
                # Using manually set theta
                if self._verbose: print "  Using manually set theta"
                self._clf = GaussianProcess(theta0=theta, nugget=nugget)
                pass
            
            self._clf.fit(X_fit, y_fit)
            
            print "  Best value(s) of theta found:", self.theta()
            
            # Store lengths scales for each toy experiment
            thetas.append(self.theta())
            self._toyclfs.append(self._clf)
            pass

        print ""
        print "== Mean toys thetas:", np.mean(thetas, axis=0)
        print "== +1 sigma range:  ", np.percentile(thetas,        84.135, axis=0)
        print "== -1 sigma range:  ", np.percentile(thetas, 100. - 84.135, axis=0)
        print "== +2 sigma range:  ", np.percentile(thetas,        97.725, axis=0)
        print "== -2 sigma range:  ", np.percentile(thetas, 100. - 97.725, axis=0)

        print ""
        self.fit()
        base_theta = self.theta()
        print "== Base theta:", base_theta

        print ""
        percx = np.sum(base_theta[0] > np.array([th[0] for th in thetas])) / float(len(thetas))
        percy = np.sum(base_theta[1] > np.array([th[1] for th in thetas])) / float(len(thetas))
        print "== Percentile of base fit wrt. toy experiments along x-axis: %.1f%%" % (percx * 100.)
        print "== Percentile of base fit wrt. toy experiments along y-axis: %.1f%%" % (percy * 100.)

        print ""
        pvaluex = 2. * min(percx, 1 - percx)
        pvaluey = 2. * min(percy, 1 - percy)
        print "== p-value (two-sided) for base fit wrt. toy experiments along x-axis: %.1f%%" % (pvaluex * 100.)
        print "== p-value (two-sided) for base fit wrt. toy experiments along y-axis: %.1f%%" % (pvaluey * 100.)

        return thetas


    def fullfit (self):
        """ Perform full, three-step fit as used in the di-jet + ISR search. """

        # Check(s)
        assert 'massbins' in self._config, "Need 'massbins' entry in config to perform full fit."

        # Pass/fail masks
        msk_data_pass = self._config['pass'](self._data)
        msk_data_fail = ~msk_data_pass
        
        # -- First 20% fit
        print "First 20% fit"
        self.window = 0.2
        self.fit()
        theta = self.theta()
        params = self._clf.get_params(deep=True)
        print "  -- Optimal theta for 20% fit:", theta

        # -- 30% validation fit
        print "30% validation fit"
        self.window = 0.3
        self.fit(theta=theta)
        w_nom  = self.weights(self._data[msk_data_fail])
        w_up   = self.weights(self._data[msk_data_fail], shift=+1)
        w_down = self.weights(self._data[msk_data_fail], shift=-1)
        

        # Compare data to estimate in VRs
        c = ap.canvas(batch=True)
        
        bins = self._config['massbins']

        # -- Main pad
        h_est_nom  = c.hist(self._data['m'][msk_data_fail], bins=bins, weights=self._data['weight'][msk_data_fail] * w_nom,  display=False)
        h_est_up   = c.hist(self._data['m'][msk_data_fail], bins=bins, weights=self._data['weight'][msk_data_fail] * w_up,   display=False)
        h_est_down = c.hist(self._data['m'][msk_data_fail], bins=bins, weights=self._data['weight'][msk_data_fail] * w_down, display=False)
        h_data     = c.plot(self._data['m'][msk_data_pass], bins=bins, weights=self._data['weight'][msk_data_pass],          display=False)

        # -- Data/estimate agreement
        bincentres = bins[:-1] + 0.5 * (bins[1] - bins[0])
        a_data = hist2array(h_data)

        a_est  = hist2array(h_est_nom)
        a_stat = np.array(map(h_est_nom.GetBinError, range(1, h_est_nom.GetXaxis().GetNbins() + 1))).reshape(a_est.shape)

        a_up   = hist2array(h_est_up)
        a_down = hist2array(h_est_down)
        a_syst = np.maximum(np.abs(a_up - a_est), np.abs(a_up - a_est))

        msk_VR = (np.abs(bincentres - self._mass) / self._mass < 0.3) & \
                 (np.abs(bincentres - self._mass) / self._mass > 0.2)
        a_dev = np.abs(a_data - a_est) / np.sqrt(np.square(a_stat) + np.square(a_syst) + eps)
        delta = np.mean(a_dev[msk_VR])

        self._shift = max(1, delta)
        print "  -- Data-estimate agreement in VRs:", delta, "-->", self._shift


        # -- Second 20% fit
        print "Second 20% fit"
        self.window = 0.2
        self.fit(theta=theta)
        #self._clf.set_params(**params)
        
        # ...
        self._fullfitted = True
        return


    def weights (self, data, shift=0, exact=False):
        """ Get transfer factor (TF) weights for data, optionally with multiple sigma shift. """

        # Check(s)
        assert self._fitted, "Must have called 'fit' before 'weights'."

        # Mesh (fine)
        X1, X2 = np.meshgrid(*self._config['axes_fine']) # *self._config['axes_fine'])
        mesh = np.vstack((X1.ravel(), X2.ravel())).T

        # Format input data
        if len(data.shape) == 2 and data.shape[1] == 2:
            # Array of shape (N,2) with assumed correct parameters
            X = data
        else:
            # Structured array
            x = data[self._config['params'][0]]
            y = data[self._config['params'][1]]
            X = np.column_stack((x, y))
            pass

        if exact:
            # Get _exact_ predictions for input data
            TF_pred, TF_err = asyncPredict(self._clf, X, quiet=True, eval_MSE=True)
        else:
            # Get indices to nearest point in fine mesh
            N = X.shape[0]

            xaxis = self._config['axes_fine'][0] # axes_fine
            yaxis = self._config['axes_fine'][1] # axes_fine
            idx1 = np.round(np.clip((X[:,0] - xaxis[0])/(xaxis[-1] - xaxis[0]),0,1) * (xaxis.size - 1)).astype(int)
            idx2 = np.round(np.clip((X[:,1] - yaxis[0])/(yaxis[-1] - yaxis[0]),0,1) * (yaxis.size - 1)).astype(int)

            # Get _approximate_ prediction for input data
            TF_pred, TF_err = asyncPredict(self._clf, mesh, quiet=True, eval_MSE=True, num_processes=1)
            TF_pred = TF_pred.reshape(X1.shape)[idx2, idx1].ravel()
            TF_err  = TF_err .reshape(X1.shape)[idx2, idx1].ravel()
            pass

        # Go from variance to r.m.s.
        TF_err = np.sqrt(TF_err)

        return TF_pred + shift * TF_err


    def toysweights (self, data, shift=1, exact=False):
        """ Get transfer factor (TF) weights for toy experiment data. """

        # Check(s)
        assert len(self._toyclfs) > 0, "Must have called 'toysfit' before 'toysweights'."

        # Store base classifier
        baseclf = self._clf

        # Compute weights for each toy classifier
        result = list()
        N = len(self._toyclfs)
        print ""
        for idx, clf in enumerate(self._toyclfs):
            print "== Getting weights from experiment %*d/%d." % (int(np.ceil(np.log10(N + 1))), idx + 1, N)
            self._clf = clf
            result.append( (self.weights(data, shift=0,      exact=exact), 
                            self.weights(data, shift=+shift, exact=exact), 
                            self.weights(data, shift=-shift, exact=exact)) )
            pass
        
        # Restore base classifier
        self._clf = baseclf

        return result
        


    def fullweights (self, data, exact=False):
        """ Get transfer factor (TF) weights for data, after the possible inflation of the systematics using validation regions. """

        # Check(s)
        assert self._fullfitted, "Must have called 'fullfit' before 'fullweights'."

        return self.weights(data, shift=0,            exact=exact), \
               self.weights(data, shift=+self._shift, exact=exact), \
               self.weights(data, shift=-self._shift, exact=exact)
        

    def plot (self, show=True, save=False, MC=True, prefix=''):
        """ ... """

        # Style
        plt.style.use('ggplot')
        plt.close('all')

        # Check(s)
        assert self._fitted, "Must have called 'fit' before 'plot'."

        if (not show) and (not save):
            print "Niether showing nor saving plot, so why bother making it."
            return

        # Try to create target directory, if necessary
        dirs = prefix.split('/')
        if len(dirs) > 1:
            check_make_dir('/'.join(dirs[:-1]))
            pass


        # (1) TF profile and residuals
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

        X1, X2 = np.meshgrid(*self._config['centres']) # (*self._config['axes'])

        # Get TF weights on mesh
        X = np.column_stack((X1.ravel(), X2.ravel()))
        TF_pred = asyncPredict(self._clf, X, quiet=True, num_processes=1).reshape(self._TF_CR_mean.shape)
        
        X1, X2 = np.meshgrid(*self._config['axes']) # (*self._config['axes'])

        # Compute pulls
        TF_CR_pulls = np.zeros_like(self._TF_CR_mean)
        msk = (self._TF_CR_err > 0)
        TF_CR_pulls[msk] = (self._TF_CR_mean[msk] - TF_pred[msk])/self._TF_CR_err[msk]

        TF_SR_pulls = np.zeros_like(self._TF_SR_mean)
        if self._window is not None:
            msk = (self._TF_SR_err > 0)
            TF_SR_pulls[msk] = (self._TF_SR_mean[msk] - TF_pred[msk])/self._TF_SR_err[msk]
            pass

        # Plot
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex = True, sharey = True, figsize = (10,8))
        im1 = ax1.pcolormesh(X1, X2, self._TF_CR_mean,  vmin=0,  vmax=1, cmap='magma')
        _   = ax2.pcolormesh(X1, X2, TF_pred,     vmin=0,  vmax=1, cmap='magma')
        im2 = ax3.pcolormesh(X1, X2, TF_CR_pulls, vmin=-2, vmax=2, cmap='RdBu')
        _   = ax4.pcolormesh(X1, X2, TF_SR_pulls, vmin=-2, vmax=2, cmap='RdBu')

        ax1.set_title('TF profile', size=16) # r'Data - W/Z (MC) profile', size=16)
        ax2.set_title('GP fit', size=16)
        ax3.set_title('Fit region', size=16)
        if self._mass is not None:
            ax4.set_title(r'Interp. region: $m \in %.0f \pm %.0f$ GeV (%.0f%%)' % (self._mass, self._mass * self._window, 100. * self._window), size=16)
        else:
            ax4.set_title(r'Interp. region', size=16)
            pass

        for ax in [ax1, ax2, ax3, ax4]:
            ax.set_xlim(self._config['axes'][0][0], self._config['axes'][0][-1])
            ax.set_ylim(self._config['axes'][1][0], self._config['axes'][1][-1])
            pass

        for ax in [ax3, ax4]:
            ax.set_xlabel(r'$\rho^{\mathrm{DDT}}$', fontsize=14)
            pass
        for ax in [ax1, ax3]:
            ax.set_ylabel(r'$\log(p_{\perp,J})$ [$\log(\mathrm{GeV})$]', fontsize=14)
            pass

        fig.subplots_adjust(left=0.1, right=0.85)
        cbar_ax1 = fig.add_axes([0.88, 0.53, 0.04, 0.36])
        cbar_ax2 = fig.add_axes([0.88, 0.10, 0.04, 0.36])
        fig.colorbar(im1, cax=cbar_ax1).set_label(label=r'$N_{\mathrm{pass}}/N_{\mathrm{fail}}$', size=14)
        fig.colorbar(im2, cax=cbar_ax2).set_label(label='Residual pulls', size=14)
        if save: plt.savefig('./' + prefix + ('profiles_%dGeV_pm%d.pdf' % (self._mass, self._window * 100.) if self._mass else 'profiles_full.pdf'))
        if show: plt.show()


        # (2) 1D pulls
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        
        bins = np.linspace(-10, 5, 30 + 1, endpoint=True)

        c2 = ap.canvas(batch=not show)
       
        h_CR_pulls = c2.hist(TF_CR_pulls[(self._TF_CR_err > 0)], bins=bins, normalise=True, display=False)
        h_SR_pulls = c2.hist(TF_SR_pulls[(self._TF_SR_err > 0)], bins=bins, normalise=True, display=False)

        # Fitting
        f_CR_pulls = ROOT.TF1('f_CR_pulls', 'gaus', -3, 3)
        f_SR_pulls = ROOT.TF1('f_SR_pulls', 'gaus', -3, 3)
        f_CR_pulls.SetLineColor(ROOT.kRed)
        f_SR_pulls.SetLineColor(ROOT.kBlue)
        f_CR_pulls.SetLineStyle(2)
        f_SR_pulls.SetLineStyle(2)

        h_CR_pulls.Fit('f_CR_pulls', 'QR0+')
        h_SR_pulls.Fit('f_SR_pulls', 'QR0+')

        # Drawing
        h_CR_pulls = c2.hist(h_CR_pulls, linecolor=ROOT.kRed)
        h_SR_pulls = c2.hist(h_SR_pulls, linecolor=ROOT.kBlue, display=bool(self._mass))

        f_CR_pulls.DrawCopy('LC SAME')
        if self._mass is not None:
            f_SR_pulls.DrawCopy('LC SAME')
            pass
        f_CR_pulls.SetLineColor(ROOT.kGray + 2)

        c2.padding(0.5)

        c2.xlabel("Transfer factor fit residual pull")
        c2.ylabel("Number of bins (a.u.)")

        c2.text(["#sqrt{s} = 13 TeV,  L = %.1f fb^{-1}" % (36.1)] + 
                 (["Sherpa incl. #gamma MC"] if MC else []) +
                 ["Trimmed anti-k_{t}^{R=1.0} jets",
                 "ISR #gamma selection"] +
                 ["No mass window" if not self._mass else \
                     ("Window: %2d GeV #pm %d%%" % (self._mass, self._window * 100.) if self._window == 0.3 else
                      "Window: %2d GeV #pm %d%%" % (self._mass, self._window * 100.)),
                 ], qualifier='Simulation Internal')

        # Legend
        ymax = 0.60
        xmin = 0.20
        step = 0.05
        legend = ROOT.TLegend(xmin, ymax - (4 if not self._mass else 7) * step, 0.5, ymax)
        legend.AddEntry(h_CR_pulls, "Fit region:",         'L')
        legend.AddEntry(None,       "  Mean: #scale[0.5]{ }%.2f #pm %.2f" % (f_CR_pulls.GetParameter(1), f_CR_pulls.GetParError(1)), '')
        legend.AddEntry(None,       "  Width: %.2f #pm %.2f" %              (f_CR_pulls.GetParameter(2), f_CR_pulls.GetParError(2)), '')
        if self._mass:
            legend.AddEntry(h_SR_pulls, "Interp. region:", 'L')
            legend.AddEntry(None,       "  Mean: #scale[0.5]{ }%.2f #pm %.2f" % (f_SR_pulls.GetParameter(1), f_SR_pulls.GetParError(1)), '')
            legend.AddEntry(None,       "  Width: %.2f #pm %.2f" %              (f_SR_pulls.GetParameter(2), f_SR_pulls.GetParError(2)), '')
            pass
        legend.AddEntry(f_CR_pulls, "Central fit", 'L')
        legend.Draw()

        # Save/show
        if save: c2.save('./' + prefix + ('residuals_distributions_%dGeV_pm%d.pdf' % (self._mass, self._window * 100) if self._mass else 'residuals_distributions.pdf'))
        if show: c2.show()

        return
    
    pass
