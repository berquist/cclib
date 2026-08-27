# Copyright (c) 2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
"""Unit tests for Q-Chem parser edge cases."""

import re

from cclib.parser.qchemparser import QChem


def test_orbital_energies_without_virtual_orbitals_use_last_orbital_as_homo():
    """A one-basis-function closed-shell job has no ``Virtual`` marker."""

    parser = QChem.__new__(QChem)
    parser.re_dashes_and_spaces = re.compile(r"^[\s-]+$")
    lines = iter(
        [
            "header before MOs\n",
            "MOs\n",
            "-- Occupied --\n",
            "-0.917\n",
            "--------------------------------------------------------------\n",
        ]
    )

    energies, symbols, homo = parser.parse_orbital_energies_and_symmetries(lines)

    assert energies == [-0.917]
    assert symbols == []
    assert homo == 0
