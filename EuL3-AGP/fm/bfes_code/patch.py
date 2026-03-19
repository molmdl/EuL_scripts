import re

with open('fes_ana.py', 'r') as f:
    content = f.read()

# Add calculate_kd function after constants
constants_end = content.find("np.seterr")
kd_function = '''
def calculate_kd(dg_kj):
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
        raise

'''
content = content[:constants_end] + kd_function + content[constants_end:]

# Update compute_dg to return Kd as well
old_return = "        return dg_kj, dg_kcal\n    except Exception as e:"
new_return = "        kd_m, kd_um, kd_nm = calculate_kd(dg_kj)\n        return dg_kj, dg_kcal, kd_m, kd_um, kd_nm\n    except Exception as e:"
content = content.replace(old_return, new_return)

# Update calculate_mode output
old_calc = '''        dg_kj, dg_kcal = compute_dg(
            fes_info, args.xmin, args.xmax,
            args.ymin, args.ymax, args.wref, args.rcyl
        )
        print(f"DeltaG: {dg_kj:.2f} kJ/mol ({dg_kcal:.2f} kcal/mol)")
        out_path = os.path.join(args.output_dir, 'deltaG_single.txt')
        with open(out_path, 'w') as f:
            f.write(f"{dg_kj:.2f} {dg_kcal:.2f}\\n")'''
new_calc = '''        dg_kj, dg_kcal, kd_m, kd_um, kd_nm = compute_dg(
            fes_info, args.xmin, args.xmax,
            args.ymin, args.ymax, args.wref, args.rcyl
        )
        print(f"DeltaG: {dg_kj:.2f} kJ/mol ({dg_kcal:.2f} kcal/mol)")
        print(f"Kd: {kd_um:.2e} uM ({kd_nm:.2e} nM)")
        out_path = os.path.join(args.output_dir, 'deltaG_single.txt')
        with open(out_path, 'w') as f:
            f.write(f"{dg_kj:.2f} {dg_kcal:.2f} {kd_m:.6e} {kd_um:.2e} {kd_nm:.2e}\\n")'''
content = content.replace(old_calc, new_calc)

# Update convergence_mode - update the compute_dg call
old_conv_dg = "            dg_kj, _ = compute_dg(fes_info, args.xmin, args.xmax, args.ymin, args.ymax, args.wref, args.rcyl)"
new_conv_dg = "            dg_kj, _, _, _, _ = compute_dg(fes_info, args.xmin, args.xmax, args.ymin, args.ymax, args.wref, args.rcyl)"
content = content.replace(old_conv_dg, new_conv_dg)

with open('fes_ana.py', 'w') as f:
    f.write(content)

print("Patch applied successfully")
