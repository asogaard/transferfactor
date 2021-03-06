# -*- coding: utf-8 -*-

""" Common configuration for transfer factor (TF) background estimation class.

@file:   config.py
@author: Andreas Søgaard
@date:   1 May 2017
@email:  andreas.sogaard@cern.ch
"""

import numpy as np

config = {

    # Base path for input files
    #'base_path' : '/eos/atlas/user/a/asogaard/Analysis/2016/BoostedJetISR/outputObjdef/2017-06-22/',
    #'base_path' : '/afs/cern.ch/user/a/asogaard/Analysis/2016/BoostedJetISR/AnalysisTools/outputObjdef/',
    #'base_path' : '/eos/atlas/user/a/asogaard/Analysis/2016/BoostedJetISR/outputObjdef/2017-07-23/',
    'base_path' : '/eos/atlas/user/a/asogaard/Analysis/2016/BoostedJetISR/outputObjdef/2017-08-06/',

    # Path to cross-sections file
    #'xsec_file' : '/afs/cern.ch/user/a/asogaard/public/Analysis/2016/BoostedJetISR/sampleInfo.csv',
    'xsec_file' : '/afs/cern.ch/user/a/asogaard/Analysis/2016/BoostedJetISR/AnalysisTools/share/sampleInfo.csv',

    # Name of tree from which to read input data
    'tree'      : 'BoostedJet+ISRgamma/Nominal/EventSelection/Pass/NumLargeRadiusJets/Postcut',

    # Name of tree from which to read input data _afte_ tau21DDT cut
    'finaltree' : 'BoostedJet+ISRgamma/Nominal/EventSelection/Pass/Jet_tau21DDT/Postcut',

    # Common prefix for branches in main tree; to be removed
    'prefix'    : 'Jet_',

    # Name of tree from which to read DSID and weight information
    'outputtree': 'BoostedJet+ISRgamma/Nominal/outputTree',

    # Luminosity
    'lumi'      : 36.1,

    # Parametrisation
    'params'    : ('rhoDDT', 'logpt'),

    # Pass criterion
    'pass'      : lambda data: data['tau21DDT'] < 0.5,

    # Axis definitions
    'axes'      : [
        np.linspace(1.5,         5.0,          70 + 1, endpoint=True), # x-axis: rhoDDT
        np.linspace(np.log(200), np.log(1000), 10 + 1, endpoint=True), # y-axis: log(pT)
        ],

    'axes_fine' : [
        np.linspace(1.5,         5.0,          10 * 70 + 1, endpoint=True), # x-axis: rhoDDT
        np.linspace(np.log(200), np.log(1000), 10 * 10 + 1, endpoint=True), # y-axis: log(pT)
        ],

    # Jet mas spectrum bin edges
    #'massbins' : np.linspace(50, 300, 50 + 1, endpoint=True),
    'massbins': np.array([50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 126, 132, 138, 145, 152, 159, 167, 175, 183, 191, 200, 210, 220, 232, 245, 259, 274, 300], dtype=float),

    }

# Bin centres
config['centres'] = [
    config['axes'][0][:-1] + 0.5 * np.diff(config['axes'][0]),
    config['axes'][1][:-1] + 0.5 * np.diff(config['axes'][1])
]

config['centres_fine'] = [
    config['axes_fine'][0][:-1] + 0.5 * np.diff(config['axes_fine'][0]),
    config['axes_fine'][1][:-1] + 0.5 * np.diff(config['axes_fine'][1])
]
