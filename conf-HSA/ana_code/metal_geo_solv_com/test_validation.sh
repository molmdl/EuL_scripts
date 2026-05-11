#!/bin/bash
echo "=== VALIDATION TESTS FOR metal_geo_analysis_v3.py ==="
echo ""

# Test L enantiomers
for geom in sap tsap; do
    echo "Testing phe_rrrL_${geom}..."
    python metal_geo_analysis_v3.py --system phe_rrrL_${geom} --gro test/solv/phe_rrrL_${geom}.gro --outdir test_results --no-shape --validate 2>&1 | grep -A3 "Geometry (twist)"
    echo ""
done

# Test D enantiomers
for geom in sap tsap; do
    echo "Testing phe_rrrD_${geom}..."
    python metal_geo_analysis_v3.py --system phe_rrrD_${geom} --gro test/solv/phe_rrrD_${geom}.gro --outdir test_results --no-shape --validate 2>&1 | grep -A3 "Geometry (twist)"
    echo ""
done

# Test Me ligands
echo "Testing me_sssL_sap..."
python metal_geo_analysis_v3.py --system me_sssL_sap --gro test/solv/me_sssL_sap.gro --outdir test_results --no-shape 2>&1 | grep -A3 "Geometry (twist)"
echo ""

echo "=== VALIDATION COMPLETE ==="
