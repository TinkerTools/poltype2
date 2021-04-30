## Symmetry Readme
<ol>
<li>	Initialize an array of symmetry classes based on openbabels graph invariants method via a set of graph invariant indexes using the graph theoretical distance, number of connected heavy atoms, aromatic boolean, ring boolean, atomic number, and summation of bond orders connected to the atom. </li>
<li>	Count and relabel classes in order such as 1,2,3… </li>
<li>	Loop over symmetry classes and look at the symmetry classes of atoms bonded to each atom. If two atoms were originally of the same symmetry class but bonded to atoms of different symmetry classes, then assign different symmetry classes for them. Count and relabel classes. </li>
<li>	Repeat above until the number of symmetry classes does not change. </li>
<li>	Terminal atoms (valence 1) of same element are assigned the same type. </li>
<li>	Adjust symmetry classes to start at new number (default 401). </li>
</ol>
