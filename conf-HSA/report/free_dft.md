# Condition
 Temperature:     298.150 K
 Pressure:          1.000 atm

# Distribution 
  SAP    Relative G=    0.000 kJ/mol     Boltzmann weight= 100.000 %
 TSAP    Relative G=   30.789 kJ/mol     Boltzmann weight=   0.000 %

# Concepts (translated from Shermo blog)

## Equation 1: Exponential Relationship
$$e^{-E_i / RT} = e^{-(\Delta E_i + E_{\text{Ref}}) / RT} = C e^{-\Delta E_i / RT}$$ 

## Equation 2: Boltzmann Distribution / Relative Probability
$$p_i = \frac{e^{-\Delta E_i / RT}}{\sum_{j} e^{-\Delta E_j / RT}} = \frac{Q_{i(\text{Relat})}}{Q_{(\text{Relat})}}$$ 

Where:

* $p$ is the population ratio.
* $i$ is the conformation number.
* $n_i$ is the number of molecules in the $i$-th conformation.
* $E$ refers to the energy of the conformation.
* $T$ is the temperature (K).
* $R$ is the ideal gas constant.
* $Q$ is the partition function.

In this formula, $E_{Ref}$ represents the lowest energy value (reference) among all conformations, and $\Delta E$ is the relative value. The constant $C$ corresponding to the reference value cancels out in the numerator and denominator when calculating $p$. The subscript Relat on $Q$ stands for "Relative."
