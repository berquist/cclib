#
# Copyright (c) 2018-2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""A writer for MolSSI quantum chemical JSON (QCSchema) files."""

import datetime as dt
import json
import math
from collections.abc import Mapping

from cclib.io.cjsonwriter import CJSON as CJSONWriter
from cclib.io.cjsonwriter import JSONIndentEncoder, NumpyAwareJSONEncoder
from cclib.parser.utils import convertor, find_package

import numpy as np


_found_qcschema = find_package("qcschema")
if _found_qcschema:
    import qcschema


_AU_CONVERSIONS = {
    "atomcoords": ("Angstrom", "bohr"),
    "scancoords": ("Angstrom", "bohr"),
    "vibdisps": ("Angstrom", "bohr"),
    "ccenergies": ("eV", "hartree"),
    "dispersionenergies": ("eV", "hartree"),
    "moenergies": ("eV", "hartree"),
    "mpenergies": ("eV", "hartree"),
    "scanenergies": ("eV", "hartree"),
    "scfenergies": ("eV", "hartree"),
    "etenergies": ("wavenumber", "hartree"),
    "vibanharms": ("wavenumber", "hartree"),
    "vibfreqs": ("wavenumber", "hartree"),
    "time": ("fs", "time_au"),
}
_AU_UNITS = {
    "atomcharges": "e",
    "charge": "e",
    "enthalpy": "hartree",
    "freeenergy": "hartree",
    "zpve": "hartree",
    "grads": "hartree/bohr",
    "etdips": "ebohr",
    "etveldips": "ebohr",
    "etmagdips": "ebohr",
    "moments": "a.u.",
}
_DIMENSIONLESS = {
    "aonames",
    "aooverlaps",
    "atombasis",
    "atomnos",
    "atomspins",
    "coreelectrons",
    "etoscs",
    "etsecs",
    "etsyms",
    "fonames",
    "fooverlaps",
    "fragnames",
    "frags",
    "homos",
    "mocoeffs",
    "mosyms",
    "mult",
    "natom",
    "nbasis",
    "nmo",
    "nocoeffs",
    "nooccnos",
    "nsocoeffs",
    "nsooccnos",
    "optdone",
    "optstatus",
    "scannames",
    "vibsyms",
}
_NO_UNIT_LABEL = {
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
}
_NATIVE_UNITS = {
    "atommasses": "Da",
    "entropy": "hartree/(particle*K)",
    "nmrcouplingtensors": "Hz",
    "nmrtensors": "ppm",
    "pressure": "atm",
    "rotconsts": "GHz",
    "temperature": "K",
    "vibfconsts": "mDyne/Angstrom",
    "vibirs": "km/mol",
    "vibramans": "Angstrom^4/Da",
    "vibrmasses": "Da",
}


def _convert_tree(value, fromunits, tounits):
    if isinstance(value, np.ndarray):
        return _convert_tree(value.tolist(), fromunits, tounits)
    if isinstance(value, np.generic):
        return _convert_tree(value.item(), fromunits, tounits)
    if isinstance(value, (list, tuple)):
        return [_convert_tree(item, fromunits, tounits) for item in value]
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return convertor(value, fromunits, tounits)
    return value


def _json_key(key):
    if isinstance(key, str):
        return key
    return "__cclib_key__:" + json.dumps(_json_safe(key), sort_keys=True, separators=(",", ":"))


def _json_safe(value):
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, dt.timedelta):
        return {"__cclib_timedelta_seconds__": value.total_seconds()}
    if isinstance(value, (dt.datetime, dt.date, dt.time)):
        return value.isoformat()
    if isinstance(value, float):
        return value if math.isfinite(value) else {"__cclib_float__": str(value)}
    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, Mapping):
        return {_json_key(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    raise TypeError(f"Cannot serialize cclib extras value {type(value)!r}")


def validate_qcschema_output(qcschema_dict):
    """Validate a QCSchema output dict against the MolSSI JSON schema.

    Fields that QCElemental defines but the pinned MolSSI schema version
    doesn't know about yet are dropped first, since the schema forbids
    additional properties.  Currently that is only wavefunction "restricted".
    """
    if not _found_qcschema:
        raise ImportError("The qcschema package is required to validate QCSchema output")
    wavefunction = {
        key: value
        for key, value in qcschema_dict.get("wavefunction", {}).items()
        if key != "restricted"
    }
    qcschema.validate({**qcschema_dict, "wavefunction": wavefunction}, schema_type="output")


def _unit_for(attribute):
    if attribute in _AU_CONVERSIONS:
        return _AU_CONVERSIONS[attribute][1]
    if attribute in _AU_UNITS:
        return _AU_UNITS[attribute]
    if attribute in _DIMENSIONLESS:
        return "dimensionless"
    if attribute in _NATIVE_UNITS:
        return _NATIVE_UNITS[attribute]
    return "unknown"


class QCSchemaWriter(CJSONWriter):
    """A writer for QCSchema files."""

    def __init__(self, ccdata, *args, **kwargs):
        super().__init__(ccdata, *args, **kwargs)

    def as_dict(self, validate=True):
        """Build QCSchema output using JSON-serializable native Python values.

        Args:
            validate: Validate the output against the MolSSI JSON schema.
        """
        metadata = self.ccdata.metadata

        qcschema_dict = {
            "schema_name": "qcschema_output",
            "schema_version": 1,
            "molecule": {
                "geometry": convertor(self.ccdata.atomcoords[-1], "Angstrom", "bohr")
                .flatten()
                .tolist(),
                "molecular_charge": self.ccdata.charge,
                "molecular_multiplicity": self.ccdata.mult,
                "schema_name": "qcschema_molecule",
                "schema_version": 2,
                "symbols": [self.pt.element[atomno] for atomno in self.ccdata.atomnos],
                "validated": True,
            },
            "provenance": {
                "creator": metadata["package"],
                "version": metadata["package_version"],
                "routine": f"{self.__module__}.{self.__class__.__qualname__}",
            },
            "success": metadata["success"],
            "error": None,
            "stdout": None,
            "stderr": None,
        }

        qcschema_dict["extras"] = {}
        for attribute in self.ccdata._attrlist:
            if attribute != "metadata" and hasattr(self.ccdata, attribute):
                value = getattr(self.ccdata, attribute)
                if attribute in _AU_CONVERSIONS:
                    value = _convert_tree(value, *_AU_CONVERSIONS[attribute])
                qcschema_dict["extras"][attribute] = value
                if (
                    attribute
                    not in set(_AU_CONVERSIONS) | set(_AU_UNITS) | _DIMENSIONLESS | _NO_UNIT_LABEL
                ):
                    qcschema_dict["extras"][f"{attribute}_unit"] = _unit_for(attribute)

        # TODO This should be derived from a parsed job type.  It's also not
        # quite right when considering that many jobs (for example, geometry
        # optimizations) are actually procedures and not drivers.
        if hasattr(self.ccdata, "hessian"):
            driver = "hessian"
        elif hasattr(self.ccdata, "grads"):
            driver = "gradient"
        else:
            driver = "energy"
        # TODO what is "properties" for?
        qcschema_dict["driver"] = driver

        # FIXME
        # if "keywords" in metadata:
        #     qcschema_dict["keywords"] = set(metadata["keywords"])
        # else:
        qcschema_dict["keywords"] = dict()

        # TODO methods should be an enum
        if metadata["methods"]:
            method = metadata["methods"][-1]
            if method == "DFT":
                method = metadata["functional"]
        else:
            # FIXME
            method = "HF"

        basis_set_name = metadata["basis_set"].lower()

        qcschema_dict["model"] = {"method": method.lower(), "basis": basis_set_name}

        scf_total_energy = convertor(self.ccdata.scfenergies[-1], "eV", "hartree")
        mp2_correlation_energy = None
        mp2_total_energy = None
        ccsd_correlation_energy = None
        ccsd_total_energy = None

        if method == "HF":
            return_energy = scf_total_energy
        elif metadata["methods"][-1] == "DFT":
            return_energy = scf_total_energy
        elif method == "MP2":
            mp2_total_energy = convertor(self.ccdata.mpenergies[-1][-1], "eV", "hartree")
            mp2_correlation_energy = convertor(
                self.ccdata.mpenergies[-1][-1] - self.ccdata.scfenergies[-1], "eV", "hartree"
            )
            return_energy = mp2_total_energy
        elif method == "CCSD":
            if hasattr(self.ccdata, "mpenergies"):
                mp2_total_energy = convertor(self.ccdata.mpenergies[-1][-1], "eV", "hartree")
                mp2_correlation_energy = convertor(
                    self.ccdata.mpenergies[-1][-1] - self.ccdata.scfenergies[-1], "eV", "hartree"
                )
            ccsd_total_energy = convertor(self.ccdata.ccenergies[-1], "eV", "hartree")
            ccsd_correlation_energy = convertor(
                self.ccdata.ccenergies[-1] - self.ccdata.scfenergies[-1], "eV", "hartree"
            )
            return_energy = ccsd_total_energy
        else:
            raise RuntimeError(f"Don't know what to do with method {method}")

        scf_dipole_moment = None
        if hasattr(self.ccdata, "moments"):
            # FIXME We currently have no easy way to tell if the moments from
            # a correlated calculation are true correlated (relaxed) density
            # moments, or just SCF density moments.
            dipole_moment = convertor(self.ccdata.moments[1], "Debye", "ebohr")
            scf_dipole_moment = dipole_moment.tolist()

        qcschema_dict["properties"] = {
            "calcinfo_nalpha": self.ccdata.homos[0] + 1,
            "calcinfo_natom": self.ccdata.natom,
            "calcinfo_nbasis": self.ccdata.nbasis,
            "calcinfo_nbeta": self.ccdata.homos[-1] + 1,
            "calcinfo_nmo": self.ccdata.nmo,
            "return_energy": return_energy,
            "scf_iterations": self.ccdata.scfvalues[-1].shape[0],
            "scf_total_energy": scf_total_energy,
            # TODO These properties aren't parsed yet.
            # ccsd_dipole_moment
            # ccsd_doubles_energy
            # ccsd_iterations
            # ccsd_opposite_spin_correlation_energy
            # ccsd_prt_pr_correlation_energy
            # ccsd_prt_pr_dipole_moment
            # ccsd_prt_pr_total_energy
            # ccsd_same_spin_correlation_energy
            # ccsd_singles_energy
            # ccsdt_correlation_energy
            # ccsdt_dipole_moment
            # ccsdt_iterations
            # ccsdt_total_energy
            # ccsdtq_correlation_energy
            # ccsdtq_dipole_moment
            # ccsdtq_iterations
            # ccsdtq_total_energy
            # mp2_dipole_moment
            # mp2_doubles_energy
            # mp2_opposite_spin_correlation_energy
            # mp2_same_spin_correlation_energy
            # mp2_singles_energy
            # nuclear_repulsion_energy
            # scf_one_electron_energy
            # scf_two_electron_energy
            # scf_vv10_energy
            # scf_xc_energy
        }
        if hasattr(self.ccdata, "dispersionenergies"):
            qcschema_dict["properties"]["scf_dispersion_correction_energy"] = convertor(
                self.ccdata.dispersionenergies[-1], "eV", "hartree"
            )
        if scf_dipole_moment is not None:
            qcschema_dict["properties"]["scf_dipole_moment"] = scf_dipole_moment
        if mp2_correlation_energy is not None:
            qcschema_dict["properties"].update(
                {
                    "mp2_correlation_energy": mp2_correlation_energy,
                    "mp2_total_energy": mp2_total_energy,
                }
            )
        if ccsd_correlation_energy is not None:
            qcschema_dict["properties"].update(
                {
                    "ccsd_correlation_energy": ccsd_correlation_energy,
                    "ccsd_total_energy": ccsd_total_energy,
                }
            )

        # QCElemental's WavefunctionProperties requires "restricted", but the
        # MolSSI JSON schemas (both v2 and dev) don't have it yet, so
        # validate_qcschema_output has to drop it.
        qcschema_dict["wavefunction"] = {
            "basis": {"name": basis_set_name, "center_data": {}, "atom_map": []},
            "restricted": len(self.ccdata.homos) == 1,
        }

        has_beta = len(self.ccdata.homos) == 2

        if hasattr(self.ccdata, "gbasis"):
            pass

        # FIXME Again, since we currently do not know if the (natural) orbital
        # coefficients and eigenvalues come from diagonalizing a (correlated)
        # density matrix, we assume that they are from SCF.
        if hasattr(self.ccdata, "moenergies"):
            qcschema_dict["wavefunction"]["scf_eigenvalues_a"] = convertor(
                self.ccdata.moenergies[0], "eV", "hartree"
            ).tolist()
            if has_beta:
                qcschema_dict["wavefunction"]["scf_eigenvalues_b"] = convertor(
                    self.ccdata.moenergies[1], "eV", "hartree"
                ).tolist()

        if hasattr(self.ccdata, "mocoeffs"):
            mocoeffs_a = self.ccdata.mocoeffs[0]
            # Sometimes we don't parse all MOs and fill missing ones with NaN,
            # which isn't allowed by QCSchema.
            if not np.isnan(mocoeffs_a).any():
                qcschema_dict["wavefunction"]["scf_orbitals_a"] = mocoeffs_a.transpose().tolist()
            if has_beta:
                mocoeffs_b = self.ccdata.mocoeffs[1]
                if not np.isnan(mocoeffs_b).any():
                    qcschema_dict["wavefunction"]["scf_orbitals_b"] = (
                        mocoeffs_b.transpose().tolist()
                    )

        if driver == "energy":
            return_result = return_energy
        elif driver == "gradient":
            return_result = self.ccdata.grads[-1].flatten().tolist()
        elif driver == "hessian":
            return_result = self.ccdata.hessian.flatten().tolist()
        else:
            raise RuntimeError(f"Don't know driver {driver}")
        qcschema_dict["return_result"] = return_result

        # Normalize once at the serialization boundary so direct callers do not
        # need a NumPy-aware JSON encoder.
        qcschema_dict = _json_safe(qcschema_dict)

        if validate:
            validate_qcschema_output(qcschema_dict)

        return qcschema_dict

    def generate_repr(self):
        """Generate the QCSchema representation of the logfile data."""
        qcschema_dict = self.as_dict()
        if self.terse:
            return json.dumps(qcschema_dict, cls=NumpyAwareJSONEncoder)
        else:
            return json.dumps(qcschema_dict, cls=JSONIndentEncoder, sort_keys=True, indent=4)
