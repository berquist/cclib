#
# Copyright (c) 2020-2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for the QCSchema writer."""

import json
import os
import unittest
from pathlib import Path

import cclib
from cclib.parser.data import ccData
from cclib.parser.utils import find_package


_found_qcelemental = find_package("qcelemental")
_found_qcschema = find_package("qcschema")


__filedir__ = os.path.dirname(__file__)
__filepath__ = Path(os.path.realpath(__filedir__))
__datadir__ = __filepath__.joinpath("..", "..").resolve()


def validate_output(outputpath):
    data = cclib.io.ccread(str(outputpath))
    writer = cclib.io.qcschemawriter.QCSchemaWriter(data)
    cclib.io.qcschemawriter.validate_qcschema_output(writer.as_dict())


def hf_water_ccdata():
    """cclib's internal representation of Q-Chem 6.2 HF_water.out."""
    return ccData(
        {
            "atomcoords": [
                [
                    [0.0, 0.0, 0.11464268],
                    [-0.7538139, 0.0, -0.45857073],
                    [0.7538139, 0.0, -0.45857073],
                ]
            ],
            "atomnos": [8, 1, 1],
            "charge": 0,
            "homos": [4],
            "metadata": {
                "package": "QChem",
                "package_version": "6.2.2",
                "methods": ["DFT"],
                "functional": "hf",
                "basis_set": "3-21G",
                "success": True,
            },
            "moments": [[0.0, 0.0, 0.0], [0.0, 0.0, -2.4181]],
            "mult": 1,
            "natom": 3,
            "nbasis": 13,
            "nmo": 13,
            "moenergies": [
                [
                    -555.79253965,
                    -36.29998766,
                    -18.85748984,
                    -14.61251377,
                    -13.06146482,
                    7.26543981,
                    9.95936693,
                    32.7625076,
                    35.81018273,
                    48.51789954,
                    50.66759896,
                    54.91257503,
                    84.84509859,
                ]
            ],
            "scfenergies": [-2056.7703591],
            "scfvalues": [
                [
                    [0.244],
                    [0.0432],
                    [0.0204],
                    [0.00213],
                    [0.000437],
                    [3.75e-05],
                    [2.73e-06],
                    [6.39e-07],
                    [1.06e-07],
                    [1.13e-08],
                    [9.66e-10],
                ]
            ],
        }
    )


class QCSchemaWriterTest(unittest.TestCase):
    @unittest.skipUnless(_found_qcelemental, "qcelemental is not installed")
    def test_qcelemental_atomic_result_hf_water_fixture(self):
        """HF-water ccData fixture generates a QCElemental QCSchema v1 result."""
        import qcelemental as qcel

        result = qcel.models.AtomicResult(
            **cclib.io.qcschemawriter.QCSchemaWriter(hf_water_ccdata()).as_dict(validate=False)
        )

        self.assertEqual(result.model.method, "hf")
        self.assertAlmostEqual(result.molecule.geometry[0][2], 0.21664326737869366)
        self.assertAlmostEqual(result.properties.return_energy, -75.58491988999288)
        self.assertEqual(result.properties.scf_dipole_moment.shape, (3,))

    def test_flat_extras_export_all_ccdata_attributes_in_documented_units(self):
        """Extras preserve each ccData attribute with a sibling unit entry."""
        data = hf_water_ccdata()
        data.setattributes(
            {
                "gbasis": [],
                "etrotats": [0.0],
                "geotargets": [0.001],
                "geovalues": [0.001],
                "hessian": [[1.0]],
                "polarizabilities": [[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]],
                "scanparm": [(0, 0.0)],
                "scftargets": [[1.0e-6]],
                "scfvalues": [[[0.244]]],
                "transprop": {"transition": ([1.0], [0.1])},
            }
        )
        result = cclib.io.qcschemawriter.QCSchemaWriter(data).as_dict(validate=False)
        extras = result["extras"]

        for attribute in data._attrlist:
            if hasattr(data, attribute) and attribute != "metadata":
                self.assertIn(attribute, extras)
        self.assertNotIn("metadata", extras)
        self.assertNotIn("metadata_unit", extras)
        self.assertAlmostEqual(extras["atomcoords"][0][0][2], 0.21664326737869366)
        self.assertAlmostEqual(extras["scfenergies"][0], -75.58491988999288)
        self.assertAlmostEqual(extras["moenergies"][0][0], -20.42500000013781)
        self.assertNotIn("atomcoords_unit", extras)
        self.assertNotIn("scfenergies_unit", extras)
        self.assertNotIn("moenergies_unit", extras)
        self.assertNotIn("moments_unit", extras)
        self.assertNotIn("homos_unit", extras)
        for attribute in (
            "etrotats",
            "gbasis",
            "geotargets",
            "geovalues",
            "hessian",
            "polarizabilities",
            "scanparm",
            "scftargets",
            "scfvalues",
            "transprop",
        ):
            self.assertNotIn(f"{attribute}_unit", extras)
        json.dumps(extras, allow_nan=False)

    @unittest.skipUnless(_found_qcelemental, "qcelemental is not installed")
    def test_qcelemental_atomic_result_mp2_energy(self):
        """MP2 output uses the correlated total energy as its result."""
        import qcelemental as qcel

        fpath = __datadir__ / "data" / "QChem" / "basicQChem5.1" / "water_mp2.out"
        data = cclib.io.ccread(str(fpath))
        result = qcel.models.AtomicResult(
            **cclib.io.qcschemawriter.QCSchemaWriter(data).as_dict(validate=False)
        )

        self.assertEqual(result.model.method, "mp2")
        self.assertAlmostEqual(result.properties.return_energy, -75.00228214)

    def test_qcschema_numeric_units(self):
        """Emitted dimensional fields use QCSchema atomic units."""
        fpath = __datadir__ / "data" / "QChem" / "basicQChem5.1" / "dvb_dispersion_bp86_d3zero.out"
        data = cclib.io.ccread(str(fpath))
        result = cclib.io.qcschemawriter.QCSchemaWriter(data).as_dict(validate=False)

        self.assertAlmostEqual(result["properties"]["scf_total_energy"], -382.3281329557)
        self.assertAlmostEqual(
            result["properties"]["scf_dispersion_correction_energy"], -0.01471992326976388
        )
        self.assertAlmostEqual(result["wavefunction"]["scf_eigenvalues_a"][0], -9.743)

    def test_as_dict_is_json_serializable(self):
        """Parsed NumPy values are converted to types the json module writes."""
        fpath = __datadir__ / "data" / "QChem" / "basicQChem5.1" / "dvb_sp.out"
        data = cclib.io.ccread(str(fpath))
        result = cclib.io.qcschemawriter.QCSchemaWriter(data).as_dict(validate=False)

        json.dumps(result, allow_nan=False)
        self.assertIsInstance(result["properties"]["calcinfo_nalpha"], int)
        self.assertIsInstance(result["molecule"]["molecular_charge"], int)

    @unittest.skipUnless(_found_qcschema, "qcschema is not installed")
    def test_validate_output_b3lyp_energy(self):
        fpath = __datadir__ / "data" / "QChem" / "basicQChem5.1" / "dvb_sp.out"
        validate_output(fpath)

    @unittest.skipUnless(_found_qcschema, "qcschema is not installed")
    def test_validate_output_b3lyp_gradient(self):
        # This is not quite right, because a geometry optimization is a
        # procedure (repeated use of the gradient driver) and not a driver
        # itself. In the future, when requesting QCSChema output, a geometry
        # optimization should split into multiple force outputs (might depend
        # on #657).
        fpath = __datadir__ / "data" / "QChem" / "basicQChem5.1" / "dvb_gopt.out"
        validate_output(fpath)

    @unittest.skipUnless(_found_qcschema, "qcschema is not installed")
    def test_validate_output_b3lyp_hessian(self):
        fpath = __datadir__ / "data" / "QChem" / "basicQChem5.1" / "dvb_ir.out"
        validate_output(fpath)
