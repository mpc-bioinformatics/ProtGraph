
def _get_mass_dict(factor=1000000000, type=int):
    ''' Return a Dictionary containing the masses of each aminoacid 

        We explicitly convert them by a factor of 1 000 000 000 (default) into integers

        The values are taken from: https://proteomicsresource.washington.edu/protocols06/masses.php
    '''
    return dict(   #  In format: AA = (MONO_MASS, AVG_MASS)
        G=(type(57.021463735*factor),  type(57.05132*factor)  ),
        A=(type(71.037113805*factor),  type(71.0779*factor)   ),
        S=(type(87.032028435*factor),  type(87.0773*factor)   ),
        P=(type(97.052763875*factor),  type(97.11518*factor)  ),
        V=(type(99.068413945*factor),  type(99.13106*factor)  ),
        T=(type(101.047678505*factor), type(101.10388*factor) ),
        C=(type(103.009184505*factor), type(103.1429*factor)  ),
        L=(type(113.084064015*factor), type(113.15764*factor) ),
        I=(type(113.084064015*factor), type(113.15764*factor) ),
        N=(type(114.042927470*factor), type(114.10264*factor) ),
        D=(type(115.026943065*factor), type(115.0874*factor)  ),
        Q=(type(128.058577540*factor), type(128.12922*factor) ),
        K=(type(128.094963050*factor), type(128.17228*factor) ),
        E=(type(129.042593135*factor), type(129.11398*factor) ),
        M=(type(131.040484645*factor), type(131.19606*factor) ),
        H=(type(137.058911875*factor), type(137.13928*factor) ),
        F=(type(147.068413945*factor), type(147.17386*factor) ),
        U=(type(150.953633405*factor), type(150.3079*factor)  ),
        R=(type(156.101111050*factor), type(156.18568*factor) ),
        Y=(type(163.063328575*factor), type(163.17326*factor) ),
        W=(type(186.079312980*factor), type(186.2099*factor)  ),
        O=(type(237.147726925*factor), type(237.29816*factor) ),

        J=(type(113.084064015*factor), type(113.1594*factor)  ),
        X=(type(0.0*factor),           type(0.0*factor)       ), # Unknown Amino Acid
        Z=(type(128.55059*factor),     type(128.6231*factor)  ),
        B=(type(114.53495*factor),     type(114.5962*factor)  ),

        __start__=(type(0), type(0)),
        __end__  =(type(0), type(0))
    )
