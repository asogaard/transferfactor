# -*- coding: utf-8 -*-

""" Common configuration for transfer factor (TF) background estimation class.

@file:   config.py
@author: Andreas Søgaard
@date:   1 May 2017
@email:  andreas.sogaard@cern.ch
"""

import numpy as np

config = {

    # Parametrisation
    'params'   : ('rhoDDT', 'logpt'),

    # Pass criterion
    'pass'     : lambda data: data['tau21DDT'] < 0.5,

    # Axis definitions
    'axes'     : [
        np.linspace(1.5,         5.0,          70 + 1, endpoint=True), # x-axis: rhoDDT
        np.linspace(np.log(200), np.log(1000), 10 + 1, endpoint=True), # y-axis: log(pT)
        ],

    'axes_fine': [
        np.linspace(1.5,         5.0,          10 * 70 + 1, endpoint=True), # x-axis: rhoDDT
        np.linspace(np.log(200), np.log(1000), 10 * 10 + 1, endpoint=True), # y-axis: log(pT)
        ]

    
    }
