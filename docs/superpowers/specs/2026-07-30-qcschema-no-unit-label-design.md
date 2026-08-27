# QCSchema extras: omit labels for undocumented units

## Goal

Do not write an `*_unit` entry in QCSchema `extras` when the corresponding cclib attribute has no physical unit specified in the parsed-data table at <https://cclib.github.io/data.html>.

## Design

- Add a dedicated `_NO_UNIT_LABEL` set in `cclib/io/qcschemawriter.py`.
- Populate it from table entries whose unit column is empty: `gbasis`, `geotargets`, `geovalues`, `hessian`, `metadata`, `polarizabilities`, `scanenergies`, `scanparm`, `scftargets`, `scfvalues`, and `transprop`.
- Exclude this set alongside converted atomic-unit attributes, already-atomic-unit attributes, and dimensionless attributes when deciding whether to write a sibling `*_unit` key.
- Preserve `_DIMENSIONLESS` as a distinct classification: dimensionless is an explicit semantic unit, whereas `_NO_UNIT_LABEL` denotes that cclib does not document a unit.
- Retain existing value serialization unchanged.

## Testing

Extend QCSchema-writer tests with `ccData` attributes representing undocumented units and assert that their serialized `extras` values have no sibling `*_unit` entries. The test must fail under the current behavior and pass once the set is used in the emission condition.

## Scope

This change only controls extras unit-label emission. It does not add conversions or assign inferred units.
