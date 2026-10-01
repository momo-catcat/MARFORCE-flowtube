def itProd(modelparams):
    if modelparams.OHsource == 'Continuous':
        itProd = modelparams.Itx
    else:
        itProd = 0
    return itProd