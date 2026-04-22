# Research: Filter Positive minimizedAffinity Values

## Current Behavior

In `dock2com_2.py`, model selection happens in two functions:

### 1. `select_best_across_sdfiles()` (lines 434-561)
- Iterates through all SDF files
- For each model, checks if metric exists in scores
- Compares scores to find best (lowest for minimizedAffinity, highest for CNNscore)
- **No filtering** - considers all models with the metric

### 2. `_select_model()` (lines 1231-1236)
- Used by `sdf_pose_to_gro()` when `model_num=None`
- Filters models that have the metric, then sorts
- **No filtering** by value

### 3. `list_models_across_sdfiles()` (lines 564-613)
- Lists all models with scores
- **No filtering** by value

## minimizedAffinity Values

From `hsa0-phe_sssD.sdf`:
- All values are negative (e.g., -7.95796, -7.59528)
- Negative = favorable binding (typical for docking scores)
- Positive = unfavorable binding (should be filtered)

## User Requirement

Filter out entries where `minimizedAffinity > 0` **before** any selection logic.

## Implementation Locations

Need to add filtering in:
1. `select_best_across_sdfiles()` - main selection function
2. `_select_model()` - single-SDF selection helper
3. `list_models_across_sdfiles()` - listing function (optional, for consistency)

## Key Constants

- `LOWER_IS_BETTER = {"minimizedAffinity", "sasa"}` (line 43)
- `DEFAULT_METRIC = "minimizedAffinity"` (line 44)

## Filter Logic

```python
# Before selection, filter models:
if "minimizedAffinity" in model["scores"]:
    affinity = model["scores"]["minimizedAffinity"]
    if affinity > 0:
        continue  # Skip this model
```

This should be applied early in the loop, before any metric-specific logic.
