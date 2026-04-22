# Verification Report: Plan B (Single-Pass SASA Metric)

**Assessor:** gsd-plan-checker  
**Date:** March 26, 2026  
**Focus:** Plan B - Single-Pass SASA Metric Integration

---

## 1. Assessment of Plan B

### 1.1 Technical Feasibility

**Status: FEASIBLE with minor issues**

Plan B is technically feasible because:
- The core SASA algorithm from `quick_sasa.py` is well-tested and fast (~0.7ms per pose)
- Integration into `select_best_across_sdfiles()` maintains the existing code structure
- Coordinate conversion (Å↔nm) is handled consistently across the codebase
- Element extraction from SDF is already implemented in `parse_sdf()`

### 1.2 Preservation of Original Functionality

**Status: PRESERVED**

Plan B maintains full backward compatibility:
- Default `--metric` remains `minimizedAffinity`
- SASA calculation only triggers when explicitly requested via `--metric sasa`
- All existing output files (best.gro, com.gro, rec.itp, sys.top) are unchanged
- The original 6-step workflow remains intact

### 1.3 Dependency Identification

**Dependencies correctly identified:**

| Dependency | Source | Status in Plan |
|------------|--------|----------------|
| numpy | Existing (dock2com_1.py line 24) | ✓ Already present |
| scipy.spatial.cKDTree | quick_sasa.py | ✓ Plan includes import |
| SASA algorithm | quick_sasa.py | ✓ Plan adapts algorithm |

---

## 2. Identified Issues and Gaps

### Issue 1: Missing Element Inference Function (BLOCKER)

**Description:** The plan includes `infer_element_from_atom_name()` in the code skeleton (plan.md lines 945-994) but does NOT include adding this function to dock2com_1.py in the task list.

**Impact:** Receptor GRO parsing cannot determine element types without this function, causing SASA calculation to fail.

**Plan Task Reference:** Task 2 (GRO parser) mentions parsing but element inference is only in the skeleton, not as a deliverable.

**Recommendation:** Add `infer_element_from_atom_name()` as a required function in dock2com_1.py.

---

### Issue 2: Incomplete Task 3 Pseudocode (WARNING)

**Description:** The pseudocode for `select_best_across_sdfiles()` modification (plan.md lines 513-562) has a logic error:

```python
if metric == "sasa" or (metric not in model["scores"]):
```

This condition calculates SASA for ALL metrics not in scores (e.g., if a pose lacks CNNscore, SASA would be calculated). This is inefficient and may not be the intended behavior.

**Recommended fix:**
```python
if metric == "sasa":
    # Always calculate SASA when that's the selection metric
    if receptor_atoms:
        sasa = calculate_pose_sasa(...)
        model["scores"]["sasa"] = sasa
```

---

### Issue 3: Receptor Loading Inside Loop (WARNING)

**Description:** Plan Task 3 loads receptor atoms inside the per-model loop:

```python
for model in models:
    if metric == "sasa":
        rec_gro, _ = derive_receptor_gro_from_sdf(sdf_path, ...)
        if os.path.exists(rec_gro):
            receptor_atoms = parse_gro_for_sasa(rec_gro)  # Called per model!
```

**Impact:** If an SDF file has 32 poses, the receptor GRO would be parsed 32 times unnecessarily.

**Recommendation:** Move receptor loading outside the model loop (per SDF file).

---

### Issue 4: Missing VDW Radii Elements (MINOR)

**Description:** The VDW radii table in the plan skeleton (plan.md lines 814-820) is incomplete compared to `quick_sasa.py`. Missing elements include: He, Ne, Ar, Kr, Xe, Ga, As, Se, etc.

**Impact:** Atoms with uncommon elements may default to Carbon radius (1.70 Å), reducing SASA accuracy for unusual ligands.

**Recommendation:** Copy the full VDW radii table from `quick_sasa.py` (lines 178-189).

---

### Issue 5: No Graceful Degradation for Missing Receptor (INFO)

**Description:** Plan B lacks explicit handling when receptor GRO file is missing. If `--metric sasa` is specified but the receptor file doesn't exist, the workflow will fail.

**Comparison with Plan A:** Plan A mentions graceful degradation with warning and penalty SASA value.

**Recommendation:** Add fallback handling:
```python
if not os.path.exists(rec_gro):
    print(f"WARNING: Cannot calculate SASA without receptor GRO: {rec_gro}")
    print("  Falling back to original selection metric")
    # Fall back to default metric
```

---

### Issue 6: CLI Help Text Incomplete (INFO)

**Description:** Plan Task 4 mentions updating `--metric` help to include "sasa" but doesn't specify adding "sasa" to the metavar choices for argparse autocomplete.

**Current:** `metavar="NAME"` allows any value  
**Recommended:** Add "sasa" to choices list for better UX

---

### Issue 7: Performance Consideration - All Poses Calculate SASA (WARNING)

**Description:** Unlike Plan A (two-stage) which only calculates SASA for top-K candidates, Plan B calculates SASA for ALL poses when `--metric sasa` is used.

**Impact:**
- If N = 320 poses across 10 SDF files
- Plan B: 320 × 0.7ms = ~224ms (acceptable but higher than Plan A)
- Plan A (top-k=10): 10 × 0.7ms = ~7ms for final stage

**Note:** This is by design for Plan B (simplicity over optimization). It's documented correctly in the plan.

---

## 3. Verification Summary

### Dimension Checklist

| Dimension | Status | Notes |
|-----------|--------|-------|
| Requirement Coverage | ✓ COVERED | SASA as metric is the requirement |
| Task Completeness | ⚠ PARTIAL | Task 3 pseudocode has logic issue |
| Dependency Correctness | ✓ VALID | All dependencies identified |
| Key Links Planned | ✓ COMPLETE | GRO→SASA→selection wired |
| Scope Sanity | ✓ WITHIN | Plan B is simpler than Plan A |
| Backward Compatibility | ✓ PRESERVED | Original behavior unchanged |

### Issues Count by Severity

| Severity | Count | Issues |
|----------|-------|--------|
| Blocker | 1 | Missing element inference function |
| Warning | 3 | Incomplete pseudocode, receptor in loop, performance |
| Info | 3 | Missing VDR radii, graceful degradation, CLI help |

---

## 4. Recommendations for Improvement

### Priority 1 (Must Fix Before Implementation)

1. **Add `infer_element_from_atom_name()` function** to dock2com_1.py
   - Required for Task 2 to work
   - Copy logic from plan.md skeleton (lines 945-994)

2. **Fix Task 3 logic** - Change SASA calculation condition:
   ```python
   # Current (inefficient):
   if metric == "sasa" or (metric not in model["scores"]):
   
   # Recommended:
   if metric == "sasa":
   ```

### Priority 2 (Should Fix)

3. **Move receptor loading outside model loop:**
   ```python
   # Load once per SDF file, not per model
   rec_gro, prefix = derive_receptor_gro_from_sdf(sdf_path, ...)
   receptor_atoms = None
   if os.path.exists(rec_gro):
       receptor_atoms = parse_gro_for_sasa(rec_gro)
   
   for model in models:
       # Now use receptor_atoms without reloading
   ```

4. **Add graceful degradation** for missing receptor GRO

### Priority 3 (Nice to Have)

5. Use full VDW radii table from quick_sasa.py

6. Add "sasa" to argparse choices for `--metric`

---

## 5. Final Verdict

### VERDICT: NEEDS_REVISION

**Reason:** 1 blocker issue prevents implementation.

**Summary:**
- Plan B is technically sound and maintains backward compatibility
- The core SASA algorithm integration is feasible
- However, the missing `infer_element_from_atom_name()` function is a blocker
- The pseudocode logic issue should also be corrected

**Required Actions Before Execution:**
1. Add `infer_element_from_atom_name()` function to dock2com_1.py (NEW FUNCTION)
2. Fix the SASA calculation condition in Task 3 pseudocode
3. Move receptor loading outside model loop (efficiency)

**After fixes:** Plan B is ready for implementation.

---

## 6. Alternative Consideration

If the issues cannot be resolved in the current plan, consider implementing **Plan A (Two-Stage)** instead, which includes:
- Better performance (only top-K candidates get SASA calculated)
- More explicit graceful degradation for missing receptors
- Clearer separation of filtering vs. rescoring stages

However, Plan B remains simpler and fully achievable once the identified issues are addressed.

---

*End of Verification Report*