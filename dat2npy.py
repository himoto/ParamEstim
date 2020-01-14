import os
import re
import numpy as np


def main():
    """Convert fitparam/n/*.dat -> out/n/*.npy
    For BioMASS (https://github.com/okadalabipr/biomass)
    
    Usage
    -----
    $ python dat2npy.py
    $ mv dat2npy/out/ path_to_biomass/
    
    """
    n_file = 0
    fitparam_files = os.listdir('./fitparam')
    for file in fitparam_files:
        if re.match(r'\d', file):
            n_file += 1
    for i in range(n_file):
        os.makedirs('./dat2npy/out/%d' % (i+1), exist_ok=True)
        ith_fitparam_files = os.listdir('./fitparam/%d' % (i+1))
        for dat_file in ith_fitparam_files:
            if 'fit' in dat_file:
                """
                - fit_param%d.dat -> fit_param%d.npy
                - best_fitness.dat -> best_fitness.npy
                """
                data = np.loadtxt(
                    './fitparam/%d/%s' % (i+1, dat_file), dtype='float'
                )
            else:
                """
                - count_num.dat -> count_num.npy
                - generation.dat -> generation.npy
                """
                data = np.loadtxt(
                    'fitparam/%d/%s' % (i+1, dat_file), dtype='int'
                )
            np.save(
                './dat2npy/out/%d/' % (i+1) +
                dat_file.replace('.dat', '.npy'), data
            )
        

if __name__ == '__main__':
    main()
