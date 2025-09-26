"""
This is a program to take the output of w0RE - i.e. jackknifes (bootstrap not supported but w0RE_boot does supply enough information)
and then use pyerrors to determine the autocorrelation
Usage

python processPyErrors PathToFortJacks.csv ensembleName w0PhysMean w0PhysError w0PhysName

"""

import pyerrors as pe
import pandas as pd  # for read_csv
import sys
import os

gamAutoCorrList = [0.0, 2.0]
gamAutoCorrLab = ['None', 'Std']


def main(args):
    if len(args) != 5:
        sys.exit(f'bad number of arguments: {args}')

    # Read the csv
    df = pd.read_csv(args[0])
    # Get the w0_ij data as jackknifes
    # 0 is the means
    w0ij_Jack = df['w0ij(xig)'].values
    xig_Jack = df['xig'].values
    # Make it pyerrors
    w0ij_PE = pe.import_jackknife(w0ij_Jack, name=args[1])
    xig_PE = pe.import_jackknife(xig_Jack, name=args[1])
    # Get the w0Phys value as pyerrors
    w0Phys = pe.cov_Obs(float(args[2]), float(args[3])**2.0, name=args[4])
    # And now calculate spacings
    a_s = w0Phys / w0ij_PE
    a_t = a_s / xig_PE

    # Get the 'input/output' dir
    outDir = '/'.join(args[0].split('/')[:-1])
    print(outDir)
    descriptionBase = f'This file contains the pyerrors observables for xi_g, a_s and a_t for the data from {args[0]} and the w0Phys value (a covObs) {w0Phys} used. Ordering is xi_g, a_s, a_t, w0Phys'
    ObsList = [xig_PE, a_s, a_t, w0Phys]
    pe.input.json.dump_to_json(ObsList,
                               os.path.join(outDir, 'jacks_Pyerrors.json.gz'),
                               description=descriptionBase,
                               gz=True)

    for lab, gamAutoCorr in zip(gamAutoCorrLab, gamAutoCorrList):
        print('Amount of autocorrelation assumed is ')
        print(lab, gamAutoCorr)

        for var, varName in zip([xig_PE, a_s, a_t], ['xig', 'a_s', 'a_t']):
            var.gm(S=gamAutoCorr)
            print('variable is ', varName)
            print(var.details())
            print('final number', var)
            print('')


if __name__ == '__main__':
    main(sys.argv[1:])
