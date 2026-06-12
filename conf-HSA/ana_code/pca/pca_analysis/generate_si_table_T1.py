#!/usr/bin/env python3
"""
generate_si_table_T1.py

Generate a publication-quality LaTeX table from si_table_T1_correlations.csv,
including FDR-BH corrected p-values.

Input:  pca_analysis/si_table_T1_correlations.csv
Output: pca_analysis/si_table_T1_correlations.tex

The script:
1. Reads the existing correlation CSV
2. Computes FDR-BH corrected p-values for all p-value columns
3. Generates a booktabs-style LaTeX table with all columns including FDR
4. Uses appropriate formatting: 2 dp for correlations, 3 dp for p-values,
   scientific notation for very small values
"""

import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def format_p_value(p, threshold=0.001):
    """Format p-value: use <0.001 for very small, otherwise 3 decimal places."""
    if p < threshold:
        return "$<$0.001"
    else:
        return f"{p:.3f}"


def format_r(r):
    """Format correlation coefficient to 2 decimal places."""
    return f"{r:.2f}"


def main():
    csv_path = os.path.join(SCRIPT_DIR, "si_table_T1_correlations.csv")
    tex_path = os.path.join(SCRIPT_DIR, "si_table_T1_correlations.tex")
    
    if not os.path.exists(csv_path):
        print(f"ERROR: {csv_path} not found")
        return
    
    df = pd.read_csv(csv_path)
    print(f"Read {len(df)} rows from {csv_path}")
    print(f"Columns: {list(df.columns)}")
    
    # Extract p-values for FDR correction
    # We have: p-value (Pearson), Spearman p, Partial p
    p_pearson = df['p-value'].values
    p_spearman = df['Spearman p'].values
    p_partial = df['Partial p'].values
    
    # FDR-BH correction for each set of p-values (6 tests per set)
    n_tests = len(df)
    
    _, p_pearson_fdr, _, _ = multipletests(p_pearson, method='fdr_bh')
    _, p_spearman_fdr, _, _ = multipletests(p_spearman, method='fdr_bh')
    _, p_partial_fdr, _, _ = multipletests(p_partial, method='fdr_bh')
    
    # Add FDR columns to dataframe
    df['Pearson p (FDR)'] = p_pearson_fdr
    df['Spearman p (FDR)'] = p_spearman_fdr
    df['Partial p (FDR)'] = p_partial_fdr
    
    # Build LaTeX table
    tex_lines = []
    tex_lines.append(r"\begin{table}[htbp]")
    tex_lines.append(r"\centering")
    tex_lines.append(r"\caption{PCA--Binding Correlation Statistics with FDR-BH Correction}")
    tex_lines.append(r"\label{tab:correlations}")
    tex_lines.append(r"\begin{tabular}{llrrrrrrrrrrr}")
    tex_lines.append(r"\toprule")
    
    # Header
    header = (
        r"Metric & PC & Pearson $r$ & $p$ & 95\% CI lower & 95\% CI upper "
        r"& Spearman $\rho$ & $p$ (Spearman) & Partial $r$ & $p$ (Partial) "
        r"& Pearson $p_{\text{FDR}}$ & Spearman $p_{\text{FDR}}$ & Partial $p_{\text{FDR}}$ \\"
    )
    tex_lines.append(header)
    tex_lines.append(r"\midrule")
    
    # Data rows
    for _, row in df.iterrows():
        metric = row['Metric']
        pc = row['PC']
        
        pearson_r = format_r(row['Pearson r'])
        pearson_p = format_p_value(row['p-value'])
        ci_lower = format_r(row[r'95\% CI lower'])
        ci_upper = format_r(row[r'95\% CI upper'])
        spearman_r = format_r(row['Spearman ρ'])
        spearman_p = format_p_value(row['Spearman p'])
        partial_r = format_r(row['Partial r'])
        partial_p = format_p_value(row['Partial p'])
        pearson_fdr = format_p_value(row['Pearson p (FDR)'])
        spearman_fdr = format_p_value(row['Spearman p (FDR)'])
        partial_fdr = format_p_value(row['Partial p (FDR)'])
        
        # Escape delta for LaTeX
        if 'ΔG' in metric or '\u0394G' in metric:
            metric_tex = r'$\Delta G$'
        else:
            metric_tex = metric
        
        line = (
            f"{metric_tex} & {pc} & {pearson_r} & {pearson_p} & {ci_lower} & {ci_upper} "
            f"& {spearman_r} & {spearman_p} & {partial_r} & {partial_p} "
            f"& {pearson_fdr} & {spearman_fdr} & {partial_fdr} \\\\"
        )
        tex_lines.append(line)
    
    tex_lines.append(r"\bottomrule")
    tex_lines.append(r"\end{tabular}")
    tex_lines.append("")
    tex_lines.append(r"\raggedright")
    tex_lines.append(r"\footnotesize")
    tex_lines.append(r"FDR-BH correction applied separately to each set of 6 p-values (Pearson, Spearman, Partial). "
                     r"$p_{\text{FDR}} < 0.05$ indicates significance after controlling the false discovery rate. "
                     r"$\Delta G$ = TOTAL\_corrected binding free energy; Contacts = mean total contacts. "
                     r"95\% CI computed via Fisher z-transform ($n = 16$)."
                     )
    tex_lines.append(r"\end{table}")
    
    # Write
    with open(tex_path, 'w') as f:
        f.write('\n'.join(tex_lines) + '\n')
    
    print(f"Wrote LaTeX table to {tex_path}")
    
    # Also update the CSV with FDR columns
    df_updated = df.copy()
    # Rename FDR columns for CSV compatibility
    df_updated = df_updated.rename(columns={
        'Pearson p (FDR)': 'Pearson p (FDR)',
        'Spearman p (FDR)': 'Spearman p (FDR)',
        'Partial p (FDR)': 'Partial p (FDR)',
    })
    csv_updated_path = os.path.join(SCRIPT_DIR, "si_table_T1_correlations_with_fdr.csv")
    df_updated.to_csv(csv_updated_path, index=False)
    print(f"Wrote updated CSV with FDR columns to {csv_updated_path}")
    
    # Print summary
    print("\n--- FDR Summary ---")
    for col, fdr_col in [('p-value', 'Pearson p (FDR)'), 
                         ('Spearman p', 'Spearman p (FDR)'),
                         ('Partial p', 'Partial p (FDR)')]:
        significant = (df[fdr_col] < 0.05).sum()
        print(f"  {col}: {significant}/{len(df)} significant after FDR correction")


if __name__ == "__main__":
    main()
