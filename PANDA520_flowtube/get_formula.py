# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 18:52:58 2022

@author: jiali
"""


# convert the species name

def get_formula(plot_spec):
    formula = []
    for s in plot_spec:
        j=-1
        for i in range(len(s)):
            add_str = '_'
            j=j+1
            if s[j].isdigit():
                add_str = add_str + s[j]
                s = s.replace(s[j], add_str)
                j=j+1

        formula.append(s)
    return (formula)
