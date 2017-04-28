#!/usr/bin/env python
import numpy as np
import io
import math
import matplotlib.pyplot as plt
import lebedev
from spherical_harmonics import sh

def metric_sh(
    N,
    L,
    ):

    leb = lebedev.Lebedev.build(N)
    w = leb.w
    theta = leb.theta
    phi = leb.phi

    Y = sh.sh(theta,phi,L)

    M = np.zeros((len(Y),)*2)
    
    lm_keys = list(sorted(Y.keys()))

    for lm1ind, lm1 in enumerate(lm_keys):
        for lm2ind, lm2 in enumerate(lm_keys):
            M[lm1ind,lm2ind] = np.sum(w * Y[lm1] * Y[lm2])

    
    print np.max(np.abs(M - np.eye(len(lm_keys))))
    print leb.order

metric_sh(
    N=170,
    L=10,
    )



    




