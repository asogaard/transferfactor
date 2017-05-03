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
    'base_path' : '/eos/atlas/user/a/asogaard/Analysis/2016/BoostedJetISR/outputObjdef/2017-04-27/',

    # Path to cross-sections file
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
        ]

    
    }
