## Python program that computes Morse-Bott spectral sequence for symplectic homology
Consider the following [link](https://drive.google.com/file/d/1DW6EujfaF360wOdVaIEZrbmbc4VfO03i/view?usp=sharing) for a picture of how the $E^1$-page of the spectral sequence looks like that comes out as the output of this code.

### Input
Given by coefficient matrix: for example, the Brieskorn polynomial $x^2+y^2+z^3$ is given to the standard input as
```
Enter Polynomial

Number of Variables: 3

Enter number of monomials in polynomial: 3
Enter coefficients for each monomial: 

Enter coefficient:2
------------------------------

Enter coefficient:0
------------------------------

Enter coefficient:0
------------------------------
Next Monomial

Enter coefficient:0
------------------------------

Enter coefficient:2
------------------------------

Enter coefficient:0
------------------------------
Next Monomial

Enter coefficient:0
------------------------------

Enter coefficient:0
------------------------------

Enter coefficient:3
------------------------------
```

### Output Example
```
Reduced coefficient matrix is:  [[2, 0, 0], [0, 2, 0], [0, 0, 3]]
Polynomial is weighted homogeneous of weights [3, 3, 2, 6]
The Z-homology of the Milnor Fiber is: [1, 0, 2, 0, 0]
The Q-homology of the link of singularity is: [1, 0, 0, 1]
The middle dimension Z-torsion of the link of singularity is: [3]

Spectral Sequence for positive homology has terms:

[[[1, 0, 2]], [[1, 1, 2]], [[2, 1, 2], [3, 0, 1]], [[2, 2, 2]], [], [[3, 3, 1]]]

Mean Euler Characteristic for S^1 Equivariant SH is:

[[[1, 0, 2]], [], [[2, 1, 2], [3, 0, 1]], [], [[3, 2, 1]]]

-3/2
```
