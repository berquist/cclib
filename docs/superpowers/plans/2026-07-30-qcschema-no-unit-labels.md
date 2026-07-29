# QCSchema Undocumented Unit Labels Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Omit QCSchema extras `*_unit` labels for cclib attributes whose parsed-data table entry has no specified physical unit.

**Architecture:** `QCSchemaWriter` will use a `_NO_UNIT_LABEL` set to distinguish attributes with undocumented units from attributes whose units are known. The extras emission condition will suppress labels for that set while retaining existing value conversion and serialization.

**Tech Stack:** Python, NumPy, pytest, cclib QCSchema writer.

## Global Constraints

- Derive unitless attributes from <https://cclib.github.io/data.html>.
- Do not introduce conversions or infer units for attributes without documented units.
- Keep `_DIMENSIONLESS` distinct from `_NO_UNIT_LABEL`.
- Keep existing atomic-unit conversion behavior unchanged.

---

### Task 1: Suppress labels for documented-unitless extras

**Files:**
- Modify: `cclib/io/qcschemawriter.py:25-90`
- Modify: `test/io/testqcschemawriter.py:35-96`

**Interfaces:**
- Consumes: `ccData._attrlist` and each attribute in `QCSchemaWriter.as_dict()`.
- Produces: `extras` values without sibling `<attribute>_unit` keys for documented-unitless attributes.

- [ ] **Step 1: Write the failing test**

Add attributes without a documented unit to `hf_water_ccdata`, then assert their serialized extras have no unit labels:

```python
"geotargets": [0.001],
"hessian": [1.0],
"polarizabilities": [[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]],
"scanparm": [(0, 0.0)],
"scftargets": [[1.0e-6]],
"scfvalues": [[[0.244]]],
```

```python
for attribute in ("geotargets", "hessian", "polarizabilities", "scanparm", "scftargets", "scfvalues"):
    self.assertNotIn(f"{attribute}_unit", extras)
```

- [ ] **Step 2: Run test to verify it fails**

Run:

```bash
uv run --no-project --with pytest --with numpy --with scipy --with periodictable --with pyyaml python -m pytest test/io/testqcschemawriter.py::QCSchemaWriterTest::test_flat_extras_export_all_ccdata_attributes_in_documented_units -q
```

Expected: FAIL because `geotargets_unit`, `hessian_unit`, `polarizabilities_unit`, `scanparm_unit`, `scftargets_unit`, and `scfvalues_unit` are currently emitted as `mixed/unknown` or `unknown`.

- [ ] **Step 3: Write minimal implementation**

Define the documented-unitless set beside the other unit classifications:

```python
_NO_UNIT_LABEL = {
    "gbasis", "geotargets", "geovalues", "grads", "hessian", "metadata",
    "polarizabilities", "scanenergies", "scanparm", "scftargets", "scfvalues", "transprop",
}
```

Extend the extras label condition so it suppresses labels when `attribute in _NO_UNIT_LABEL`:

```python
if attribute not in set(_AU_CONVERSIONS) | set(_AU_UNITS) | _DIMENSIONLESS | _NO_UNIT_LABEL:
    qcschema_dict["extras"][f"{attribute}_unit"] = _unit_for(attribute)
```

- [ ] **Step 4: Run focused and module tests to verify they pass**

Run:

```bash
uv run --no-project --with pytest --with numpy --with scipy --with periodictable --with pyyaml python -m pytest test/io/testqcschemawriter.py -q
git diff --check
```

Expected: all executable tests pass, optional-schema tests may skip, and `git diff --check` has no output.

- [ ] **Step 5: Commit**

```bash
git add cclib/io/qcschemawriter.py test/io/testqcschemawriter.py
git commit -m "io: omit undocumented QCSchema unit labels"
```
