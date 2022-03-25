## BAR Readme

### BAR Protocol
* 1) For each consecutive tuple of lambdas, compute the free energy, enthalpy, entropy and error estimates via BAR
* 2) Add the total energies from each consecutive tuple of lambdas.

### BAR Output file
* In BARResults.csv, in UTF-8 format, read README_TABLES to open in excel
* Table explanation below
```
Ele-Lambda             Electrostatic Lambda
Vdw-Lambda             Vdw Lambda
Rest-Lambda            Restraints (all types of restraints) Lambda
ΔGᶜᵒᵐᵖ                 Complexation free energy (default via BAR iteration)
ΔGᶜᵒᵐᵖᵉʳʳ              Complexation free energy error
ΔGᶜᵒᵐᵖᵇᵃʳᶦᵗᵉʳ          Complexation free energy via BAR iteraiton
ΔGᶜᵒᵐᵖᵉʳʳᵇᵃʳᶦᵗᵉʳ       Complexation free energy error via BAR iteraiton
ΔGᶜᵒᵐᵖᵇᵒᵒᵗˢᵗʳᵃᵖ        Complexation free energy via bootstrapping
ΔGᶜᵒᵐᵖᵉʳʳᵇᵒᵒᵗˢᵗʳᵃᵖ     Complexation free energy error via bootstrapping
ΔGᶜᵒᵐᵖᶠʷᵈ              Complexation free energy via forward BAR iteraiton
ΔGᶜᵒᵐᵖᵉʳʳᶠʷᵈ           Complexation free energy error via forward BAR iteraiton
ΔGᶜᵒᵐᵖᵇʷᵈ              Complexation free energy via backward BAR iteraiton
ΔGᶜᵒᵐᵖᵉʳʳᵇʷᵈ           Complexation free energy error via backward BAR iteraiton
ΔHᶜᵒᵐᵖ                 Complexation enthalpy
ΔHᶜᵒᵐᵖᵉʳʳ              Complexation enthalpy error
ΔSᶜᵒᵐᵖ                 Complexation entropy
ΔSᶜᵒᵐᵖᵉʳʳ              Complexation entropy error
``` 
