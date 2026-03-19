with open('fes_ana.py', 'r') as f:
    content = f.read()

# Replace the calculate_kd function with calculate_log_k
old_func = '''def calculate_kd(dg_kj):
    """
    Calculate dissociation constant Kd from binding free energy.

    Parameters:
    dg_kj (float): Binding free energy in kJ/mol.

    Returns:
    tuple: (Kd_M, Kd_uM, Kd_nM) - Kd in Molar, micromolar, and nanomolar.

    Notes:
    Kd = exp(ΔG / RT) where ΔG is the standard binding free energy.
    For a negative ΔG (favorable binding), Kd < 1 M.
    """
    try:
        kd_molar = np.exp(dg_kj / KT)
        kd_um = kd_molar * 1e6
        kd_nm = kd_molar * 1e9
        return kd_molar, kd_um, kd_nm
    except Exception as e:
        logging.error(f"Error calculating Kd: {e}")
        raise'''

new_func = '''def calculate_log_k(dg_kj):
    """
    Calculate log10 of binding constant K from binding free energy.

    Parameters:
    dg_kj (float): Binding free energy in kJ/mol.

    Returns:
    float: log10(K) where K is the binding (association) constant in M^-1.

    Notes:
    K = exp(-ΔG / RT) where ΔG is the standard binding free energy.
    log10(K) = -ΔG / (RT * ln(10))
    For a negative ΔG (favorable binding), log10(K) > 0.
    """
    try:
        log_k = -dg_kj / (KT * np.log(10))
        return log_k
    except Exception as e:
        logging.error(f"Error calculating log K: {e}")
        raise'''

content = content.replace(old_func, new_func)

# Update compute_dg
content = content.replace("kd_m, kd_um, kd_nm = calculate_kd(dg_kj)", "log_k = calculate_log_k(dg_kj)")
content = content.replace("return dg_kj, dg_kcal, kd_m, kd_um, kd_nm", "return dg_kj, dg_kcal, log_k")

# Update calculate_mode
content = content.replace("dg_kj, dg_kcal, kd_m, kd_um, kd_nm = compute_dg(", "dg_kj, dg_kcal, log_k = compute_dg(")
content = content.replace('print(f"Kd: {kd_um:.2e} uM ({kd_nm:.2e} nM)")', 'print(f"log K: {log_k:.2f}")')
content = content.replace('f.write(f"{dg_kj:.2f} {dg_kcal:.2f} {kd_m:.6e} {kd_um:.2e} {kd_nm:.2e}\\n")', 'f.write(f"{dg_kj:.2f} {dg_kcal:.2f} {log_k:.2f}\\n")')

# Update convergence_mode
content = content.replace("dg_kj, _, _, _, _ = compute_dg(", "dg_kj, _, _ = compute_dg(")

with open('fes_ana.py', 'w') as f:
    f.write(content)

print("Patch applied successfully")
