## Torsion Fitting Readme

The program **torsionfit.py** extracts the energies from the output files of the calculations performed on the MM and QM optimized structures, along the torsion scan for each torsion parameters. Pre-fit MM energy, QM energy and angles from all rotatable bonds are collected and average energy of multiple instances of the same torsion, which can be found, repeated, in different part of molecule (same symmetry types).

The fitting procedure is performed, and repeated **at most** twice if necessary, for all rotatable bonds, using Boltzmann fitting when trying for the second time. 
Consecutive coupled torsion can be parameterized, first by parameterizing single torsions are fitted (1D torsion) then by using a spline fit the QM-MM2, where MM2 is parameters from 1D torsion fitting.

Additionally, torsion rotations are categorized by symmetry, based on the atoms which are forming that torsion, for storing fitting results (*e.g.*: C-C-C-C or H-C-C-H torsions, *etc.*).

To proceed with the fitting, first QM and MM energies are **normalized**, and the QM value is then subtracted from the MM one. If a second try for the fitting procedure is required, the Boltzmann weight defined as

```math
 e^{-E_{QM}/8}.
```

From all the QM - MM energies differences a min/max criterion is used to set the boundaries to evaluate the fitting quality. Those boundaries are identified by the amplitude between the max and min. The parameters for fitting are then initialized, for each fold of the torsion potential energy term. In addition to the min/max boundaries, also boundaries for the phase angle (0-2$\pi$) for fitting phase angle and also an additional criterion is defined for non-aromatic ring puckering torsion parameters.

The cost function to be minimized in the fitting procedure is defined from subtracting the energy of torsion parameters only from the QM-MM energies. The least squares procedure of scipy is used for this optimization. The Boltzmann weights are applied to the energy difference , when it fails at first attempt of fitting.

The fitting parameters' quality is then checked: all parameters are compared to the min/max boundary, and when the actual fitted parameter is greater than that, the fitting is repeated and boundaries modified (essentially remove one of the parameters from fitting by keeping it fixed near a lower specific value). When combinations of torsions are considered, a specific criterion is applied to first and third fold torsion terms (try to prevent parameter cancellation). Graphs are generated to also check the results of the fitting, comparing pre-fitting MM, QM and post-fitting MM energies, as functions of the torsion angle values.

If everything goes smoothly, the fitted parameters are written to the .key file with all the parameters derived. Moreover, from the QM optimized structures, tinker minimization are performed again, using the new parameters, analyzing the newly optimized MM structures.

If either the QM energy, the WBO bond index, the pre-fit MM or post-fit MM energies is missing for a specific point along the torsion scan, that point is removed from the scan.



<ol>
<li>    Copy .key_4 from parent poltype directory into qm-torsion folder </li> 
<li>    Change directory to qm-torsion </li> 
<li>    If deriving parameters for tortor, decompose all tuples of two rotatable bonds into single rotatable bonds for fitting 1D torsion prior to fitting tortor </li> 
<li>    Extract MM1 (prefit MM energy), QM energy and angles from all rotatable bonds tuples (can fit to many same time) </li> 
<ol>
<li>    If QM, or MM energy missing, remove points from both QM and MM arrays </li> 
<li>    Average energy of multiple instances of class keys (same torsion different part of molecule) </li> 
</ol>
<li>    Try to fit all rotatable bonds </li> 
<ol>
<li>    Don’t do more than twice, use Boltzmann fitting second time if first time does not pass fitting criteria </li> 
<li>    For each torsion set (can be multiple rotatable bonds) </li> 
<li>    First skip arrays with more than 1 rotatable bond (fit 1D torsion first) </li> 
<li>    For each torsion in torset </li> 
<ol>
<li>    Grab symmetry type numbers for the main torsion (first heaviest torsion (a,d in a-b-c-d attempt to find non-hydrogen)) </li> 
<li>    Grab all torsion indices a-b-c-d and associate them with the main torsion symmetry class key, for storing fitting results </li> 
</ol>
<li>    Normalize QM, and MM energy arrays </li> 
<li>    Subtract QM from MM </li>  
<li>    Define Boltzmann weights as e^(-QM/8) </li> 
<li>    Only used if first round fitting fails </li> 
<li>    Take the max and min of QM-MM array, call this max_amp (max amplitude) </li> 
<li>    If max_amp > 20 kcal/mol, then set max_amp=20 </li> 
<li>    This  will be used to define bounds for torsion fitting </li> 
<li>    Initialize parameters for fitting </li> 
<ol>
<li>    If keyword for firstfoldphase is on (fit phase in first fold torsion), then add to parameter dictionary </li> 
<li>    For each fold (cosine term), add parameter to be fit </li> 
</ol>
<li>    Initialize boundaries for fitting	</li> 
<ol>
<li>    0-2pi for fitting phase angle </li> 
<li>    -max_amp – max_amp if not refining parameters </li> 
<li>    parameter-.3*parameter , parameter+.3*parameter for refining (non-aromatic ring puckering) torsion parmameters </li> 
</ol>
<li>    Cost function </li> 
<ol>
<li>    Subtract energy of torsion parameters only from (QM-MM array). </li> 
<li>    Minimize this with optimize least_squares in scipy </li> 
<li>    If failed fitting previously, apply Boltzmann weights to the difference </li> 
</ol>
<li>    Minimize the cost function </li> 
<li>    Check fitting parameter quality </li> 
<ol>
<li>    For all parameters, if parameter is greater than max_amp </li> 
<ol>
<li>    Then modify boundaries for parameter to remain between -1.5 and 1.5 and redo fitting </li> 
</ol>
<li>    For all combinations of two parameters of same cosine term for all torsion around rotatable bond </li> 
<ol>
<li>    Compute the torsion energy for all angles, only for the torsions given and only for given cosine term </li> 
<li>    Normalize the energy array </li> 
<li>    Take the max energy and the max parameter value from combination of two parameters, compute the ratio (max energy/max parameter) </li> 
<ol>
<li>    If this is <.1, set new boundaries of -1.5, 1.5 for this parameter and redo fitting. </li> 
</ol>
<li>    If cosine term is first or third, if all parameters in combination are > 15 kcal/mol </li> 
<ol>
<li>    Then pick the first parameter and fix at new boundaries of -1.5, 1.5 and redo fitting </li> 
</ol>
</ol>
<li>    Generate plots of QM vs angle </li> 
</ol>
</ol>
<li>    Write out fitting results to keyfile </li> 
<ol>
<li>    Including phase fitting results if keyword for firstfoldfitting turned on </li> 
</ol>
<li>    Remove filenames with post occurring twice in it (for if redoing fitting with Boltzmann) </li> 
<li>    Starting from QM opt structures, redo tinker minimize and then analyze with newly fit parameters </li> 
<li>    Evaluate parameters post fitting </li> 
<ol>
<li>    For every torsion set </li> 
<ol>
<li>    Skip this step if only fitting 1D currently and torsion set has more than one rotatable bond </li> 
<li>    Grab prefit and postfit MM energies </li> 
<li>    Grab QM energies and WBO bond indices </li> 
<li>    If QM or MM1 (prefit) or MM2 (postfit) or WBO values are missing, then remove point from all arrays </li> 
<li>    Normalize energies </li> 
<li>    If using boltzman weights, normalize weights </li> 
<li>    Compute QM-MM array, compute minimum RMSE and relative RMSE (divide by QM) </li> 
<li>    Generate plots, MM1, MM2, QM, MM1+fit </li> 
<li>    If RMSE and relative RMSE fails criteria 1.8 and .2 respectively, then redo fitting </li> 
<ol>
<li>    If fails again after Boltzmann then raise error </li> 
</ol>
</ol>
</ol>
<li>    If was fitting 1D torsion and have tortor to fit </li> 
<ol>
<li>    Grab QM, MM2 energy arrays </li> 
<li>    Rearrange flattened arrays into matrix </li> 
<li>    Subtract QM-MM2 </li> 
<li>    Generate class key for tortor (symmetry type numbers) </li> 
<li>    Add copies of first row and first column so that angles repeat themselves on other sides of matrix (tinker manual says needs this) </li> 
<li>    Plot heat maps for QM, MM2, QM-MM2 </li> 
<li>    Add matrix of QM-MM2 into tinker key file </li> 
</ol>
<li>    Copy final key file to .key_5 in parent directory </li> 
</ol>
