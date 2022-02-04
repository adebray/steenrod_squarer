This repo consists of a program for calculating Steenrod squares in mod 2 cohomology: you specify the Steenrod
squares on generators and it handles multiplication and the Cartan formula for you.

**Warning**: this code works properly when the cohomology ring is free, and when there are simple relations, but
for more complicated rings there can be errors. There is definitely a bug for quaternion groups.

### Usage

`./steenrod_squarer.py config/file.json`

Here `config/file.json` is a configuration file: you specify the generators of the cohomology ring and their
Steenrod squares (along with the Steenrod squares of the Thom class, if you are working with a Thom spectrum). You
can also specify relations in the ring. There are many examples in the `config/` directory of the repo.

The program `steenrod_squarer.py` contains several useful utilities. Someday I'll organize these into a library.

* `str_to_poly`: internally, elements of a cohomology ring are represented as members of a class called
	  `polynomial`. The function `str_to_poly` takes a string representing a monomial, such as `ab`, and converts
	  it into a member of the `polynomial` class. For powers, just repeat a variable, such as `bb` for b^2. If you
	  want non-monomials such as ab + b^2, you can use the `+` operator: `str_to_poly('ab') + str_to_poly('bb')`.
	  
* `identify_a1_summands`: takes as input a list of `polynomial`s and determines which ones generate free
	  A(1)-module summands.
	  
* `compute_span_of` takes as input a list of polynomials and repeatedly applies Sq^1 and Sq^2 until eventually
	  reaching 0, and outputs a list of everything that can be reached in this way.

Please let me know if you find any bugs, or if anything seems confusing or incorrect, etc.
