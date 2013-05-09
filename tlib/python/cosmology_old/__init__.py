"""The CosmoPy Theory Module

In this module you will find all the relevant classes for theoretical
computations.  The list of things to be implemented here includes
cosmological distance calculator, power spectrum and correlation
function derivation and much more...

Current revision:
    ID:         $Id: __init__.py 56 2008-01-08 00:53:36Z neyrinck $
    Date:       $Date: 2008-01-07 14:53:36 -1000 (Mon, 07 Jan 2008) $
    Revision:   $Revision: 56 $
"""

REVISION = '$Revision: 56 $'

def import_all():
    """By default the nested modules are not imported automatically.
    Call this function if you would like to import them all.
    This may be useful for autocompletion in interactive mode."""

    import bessint
    import cex
    import example
    import fkernel
    import frw
    import halo
    import hankel
    import info
    import massfn
    import param
    import proj
    import pt
    import utils
