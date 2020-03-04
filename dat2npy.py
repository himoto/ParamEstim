import os
import shutil
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
    n_file = []
    fitparam_files = os.listdir('./fitparam')
    for file in fitparam_files:
        if re.match(r'\d', file):
            n_file.append(int(file))
    for _, nth_paramset in enumerate(n_file):
        os.makedirs(
            './dat2npy/out/{:d}'.format(
                nth_paramset
            ), exist_ok=True
        )
        ith_fitparam_files = os.listdir(
            './fitparam/{:d}'.format(
                nth_paramset
            )
        )
        for dat_file in ith_fitparam_files:
            if 'fit' in dat_file:
                """
                - fit_param%d.dat -> fit_param%d.npy
                - best_fitness.dat -> best_fitness.npy
                """
                data = np.loadtxt(
                    './fitparam/{:d}/{}'.format(
                        nth_paramset, dat_file
                    ), dtype='float'
                )
            else:
                """
                - count_num.dat -> count_num.npy
                - generation.dat -> generation.npy
                """
                data = np.loadtxt(
                    'fitparam/{:d}/{}'.format(
                        nth_paramset, dat_file
                    ), dtype='int'
                )
            np.save(
                './dat2npy/out/{:d}/'.format(
                    nth_paramset
                ) + dat_file.replace(
                    '.dat', '.npy'
                ), data
            )
        if os.path.isfile('./logs/{:d}.log'.format(nth_paramset)):
            shutil.copyfile(
                './logs/{:d}.log'.format(
                    nth_paramset
                ), './dat2npy/out/{:d}/out.log'.format(
                    nth_paramset
                )
            )

if __name__ == '__main__':
    main()
