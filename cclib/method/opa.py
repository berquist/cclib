# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Calculation of overlap population analysis based on cclib data."""

from cclib.method.population import Population

import numpy


def func(x):
    if x == 1:
        return 1
    else:
        return x + func(x - 1)


class OPA(Population):
    """Overlap population analysis."""

    def __init__(self, *args):
        super().__init__(logname="OPA", *args)

    def __str__(self):
        """Return a string representation of the object."""
        return f"OPA of {self.data}"

    def __repr__(self):
        """Return a representation of the object."""
        return f'OPA("{self.data}")'

    def calculate(self, indices=None, fupdate=0.05):
        """Perform an overlap population analysis given the results of a parser"""
        if not indices:
            # Build list of groups of orbitals in each atom for atomresults.
            if hasattr(self.data, "aonames"):
                names = self.data.aonames
            elif hasattr(self.data, "foonames"):
                names = self.data.fonames

            atoms = []
            indices = []

            name = names[0].split("_")[0]
            atoms.append(name)
            indices.append([0])

            for i in range(1, len(names)):
                name = names[i].split("_")[0]
                try:
                    index = atoms.index(name)
                except ValueError:  # not found in atom list
                    atoms.append(name)
                    indices.append([i])
                else:
                    indices[index].append(i)

        # Determine number of steps, and whether process involves beta orbitals.
        nfrag = len(indices)  # nfrag
        nstep = func(nfrag - 1)
        unrestricted = len(self.data.mocoeffs) == 2
        alpha = len(self.data.mocoeffs[0])

        self.logger.info("Creating attribute results: array[4]")
        results = [numpy.zeros([nfrag, nfrag, alpha], "d")]
        if unrestricted:
            beta = len(self.data.mocoeffs[1])
            results.append(numpy.zeros([nfrag, nfrag, beta], "d"))
            nstep *= 2

        if hasattr(self.data, "aooverlaps"):
            overlap = self.data.aooverlaps
        elif hasattr(self.data, "fooverlaps"):
            overlap = self.data.fooverlaps

        # intialize progress if available
        if self.progress:
            self.progress.initialize(nstep)

        step = 0

        for spin in range(len(self.data.mocoeffs)):
            two = numpy.array([2.0] * len(self.data.mocoeffs[spin]), "d")

            # OP_{AB,i} = \sum_{a in A} \sum_{b in B} 2 c_{ai} c_{bi} S_{ab}

            for A in range(len(indices) - 1):
                for B in range(A + 1, len(indices)):
                    if self.progress:  # usually only a handful of updates, so remove random part
                        self.progress.update(step, "Overlap Population Analysis")

                    for a in indices[A]:
                        ca = self.data.mocoeffs[spin][:, a]

                        for b in indices[B]:
                            cb = self.data.mocoeffs[spin][:, b]
                            temp = ca * cb * two * overlap[a, b]
                            results[spin][A, B] = numpy.add(results[spin][A, B], temp)
                            results[spin][B, A] = numpy.add(results[spin][B, A], temp)

                    step += 1

        temparray2 = numpy.swapaxes(results[0], 1, 2)
        self.results = [numpy.swapaxes(temparray2, 0, 1)]
        if unrestricted:
            temparray2 = numpy.swapaxes(results[1], 1, 2)
            self.results.append(numpy.swapaxes(temparray2, 0, 1))

        if self.progress:
            self.progress.update(nstep, "Done")

        return True
