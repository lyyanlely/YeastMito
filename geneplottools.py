from __future__ import division
import numpy as np
import pysam
from collections import Counter
import os
import counter_stats as cs
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
import numpy.ma as ma


def add_features(fig, features, ignore_features): 
    # add genetic features to a gene plot
    for feature in features:
        name = ''
        if feature.type == 'gene':
            color = 'C1'
            try:
                name = feature.qualifiers['gene'][0]
            except KeyError:
                continue
        elif feature.type == 'rep_origin':
            color = 'C2'
            try:
                name = feature.qualifiers['note'][0]
            except KeyError:
                pass
        else:
            continue
        if name in ignore_features:
            continue
        start = feature.location.start.position
        end = feature.location.end.position
        fig.gca().axvspan(start, end, alpha=0.25, color=color)
        gene_box = mpatches.Patch(color='C1', alpha=0.25)
        ori_box = mpatches.Patch(color='C2', alpha=0.25)
        return fig, gene_box, ori_box