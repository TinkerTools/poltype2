## Tables Readme

### Tables Protocol
* Store all program inputs and outputs (free energies, entropy, enthalpy) into dictionaries and output into CSV (UTF-8) format (viewable in excel).

* Gibbs Free Energy Change Table 
* solv=Solvation, comp=complexation, err=error, corr=corrected, uncorr=uncorrected, bind=binding, ana=analytical
* ΔG units of kcal/mol, ΔH units of kcal/mol, ΔS units of kcal/mol/K 
``` 
ΔGˢᵒˡᵛ
ΔGˢᵒˡᵛᵉʳʳ
ΔGᶜᵒᵐᵖᶜᵒʳʳ
ΔGᶜᵒᵐᵖᵘⁿᶜᵒʳʳ
ΔGᶜᵒᵐᵖᶜᵒʳʳᵉʳʳ
ΔGᵃⁿᵃᶜᵒᵐᵖᶜᵒʳʳ
ΔGᵇᶦⁿᵈᶜᵒʳʳ
ΔGᵇᶦⁿᵈᶜᵒʳʳᵉʳʳ
```
* Enthalpy, Entropy, Gibbs Energy Change Table
* Modename	
```
Writeout Freq (ps)	
Dyntamic Time Step (fs)	
Total Time (ns)	
Total Dyn Extended Time	Total Dynamic Steps	
Total Extended Dynamic Steps	
Total ARC File Space Needed
ΔGˢᵒˡᵛ
ΔGˢᵒˡᵛᵉʳʳ	
ΔHˢᵒˡᵛ	
ΔHˢᵒˡᵛᵉʳʳ	
ΔSˢᵒˡᵛ	
ΔSˢᵒˡᵛᵉʳʳ	
ΔGᶜᵒᵐᵖᶜᵒʳʳ	
ΔGᶜᵒᵐᵖᵘⁿᶜᵒʳʳ	
ΔGᶜᵒᵐᵖᶜᵒʳʳᵉʳʳ	
ΔHᶜᵒᵐᵖ	
ΔHᶜᵒᵐᵖᵉʳʳ
```
* Simulation Info Table
* Prod = Production (not equilibritation dynamics)
```
Total Atom Number	
Average Box Size	
Prod MD Ensemb	
Prod MD Time	
Prod MD Steps	
Prod MD Arc File Space	
Dynamic Writeout Frequency (ps)	
Dynamic Time Step (fs)	
Equil Time NPT	
Equil Time NVT	
Physio Counterions	
Neut Counterions
```

### Excel Instructions
* Open Excel
* Go to data tab
* On the left, click get data
* Click get data from file
* Open CSV file
* On the new window, select unicode (UTF-8) from the drop down. In delimiter, make sure it says comma.
* Click okay.
