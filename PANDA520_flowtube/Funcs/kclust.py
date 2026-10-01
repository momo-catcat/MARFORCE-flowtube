def kClust(modelparams):
    if hasattr(modelparams, 'kClust_value'):
        return modelparams.kClust_value
    else:
        return 4e-14  # default value


def kDimer(modelparams):
    """Rate constant for SA + SA -> dimer."""
    if hasattr(modelparams, 'kDimer_value'):
        return modelparams.kDimer_value
    else:
        return 4e-14


def kTrimer(modelparams):
    """Rate constant for dimer + SA -> trimer and higher clustering."""
    if hasattr(modelparams, 'kTrimer_value'):
        return modelparams.kTrimer_value
    else:
        return 5e-13
