# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 18:52:58 2022

@author: jiali
"""


# convert the species name

def get_formula(plot_spec):
    formula = []
    for s in plot_spec:
        add_str = '_'
        for i in range(len(s)):
            if s[i].isdigit():
                add_str = add_str + s[i]
                s = s.replace(s[i], add_str)
                s = '$\mathdefault{' + s + '}$'

        formula.append(s)
    return (formula)
