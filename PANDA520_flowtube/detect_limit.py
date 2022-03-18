# -*- coding: utf-8 -*-
"""
Created on Sat Nov 27 11:31:17 2021

@author: jiali
"""
import statistics


def get_detect_llimit(Sa_stage):
    # sa=merg_stage(num,A)
    sigma = statistics.stdev(Sa_stage) * 3 + Sa_stage.mean()
    sigma = "{:.2e}".format(sigma)
    return sigma
