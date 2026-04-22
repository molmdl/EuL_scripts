# Plan: Filter Positive minimizedAffinity Values

## Objective

Add filtering to exclude models with `minimizedAffinity > 0` before any selection.

## Changes Required

### 1. `select_best_across_sdfiles()` (lines 520-553)

Add filtering at the start of the model loop (after line 520):

```python
for model in models:
    # Filter: skip models with positive minimizedAffinity
    if "minimizedAffinity" in model["scores"]:
        if model["scores"]["minimizedAffinity"] > 0:
            continue
    
    # ... rest of existing logic
```

### 2. `_select_model()` (lines 1231-1236)

Add filtering before candidate list creation:

```python
def _select_model(models, metric=DEFAULT_METRIC):
    # Filter out positive minimizedAffinity first
    filtered = []
    for m in models:
        if "minimizedAffinity" in m["scores"] and m["scores"]["minimizedAffinity"] > 0:
            continue
        filtered.append(m)
    
    candidates = [m for m in filtered if metric in m["scores"]]
    # ... rest unchanged
```

### 3. `list_models_across_sdfiles()` (lines 600-613)

Optionally mark filtered models or skip them (user preference needed).

For listing, two options:
- **Option A**: Skip filtered models entirely
- **Option B**: Show filtered models with a marker (e.g., "FILTERED")

**Recommendation**: Option A (skip) - cleaner output

### 4. Update Help Text

Add note in docstring about filtering behavior.

## Implementation Order

1. Add helper function `_is_model_valid()` for reusable filter logic
2. Update `select_best_across_sdfiles()`
3. Update `_select_model()`
4. Update `list_models_across_sdfiles()`
5. Update module docstring

## Validation

After changes:
```bash
# Test with existing data (all negative, should work)
python dock2com_2.py -i lig_g.itp -s hsa*-phe_sssD.sdf --list-models

# Verify selection still works
python dock2com_2.py -i lig_g.itp -s hsa0-phe_sssD.sdf -r topol.top -t phe_sssD.mol2
```

## Edge Cases

1. **All models filtered**: Raise error "No valid models after filtering"
2. **Missing minimizedAffinity**: Keep model (don't filter if metric not present)
3. **NaN values**: Already handled (lines 537-538)
