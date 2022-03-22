#!/usr/bin/env python3

# Arun Debray
# 25 February 2020
# (updated 19 March 2021)

# Little utility script to compute Steenrod squares in cohomology rings. I don't aim for this
# to be very sophisticated -- basically an implementation of the Cartan formula for Sq^1 and Sq^2;
# then, you define your cohomology ring generators (and maybe relations... maybe those won't be
# implemented) and how Sq^1 and Sq^2 act on them, and compute.

# probably the biggest TODO is to make this more useable by reading the data from a config file and
# exporting library functions (e.g. generates_an_A(1))
# I think the config file is working but I would like to test more.
# exporting library functions is still WIP

import collections
import functools
import itertools
import json
import sys

# Probably the most correct way to implement this would be to make it data of a class, and to pass that
# around. For now, though, it's easiest to make this a global. Eventually it will be a dict (see the examples
# in this directory for what it should look like)
coh_info = None

# ("x", 2) -> x^2, but also using some Unicode superscripts
def prettyprint_exponents(var: str, power: int):
	if power == 0: return ''
	elif power == 1: return var
	# indices 0 and 1 should not happen, so are represented with '!' to indicate to the user (i.e. me)
	# that the developer (also me) made a mistake somewhere
	elif power < 10: return '%s%s' % (var, "!!²³⁴⁵⁶⁷⁸⁹"[power])
	else: return '%s^{%d}' % (var, power)

# slightly surprise this doesn't exist yet
# takes in two collections.Counter objects, which are to be interpreted as "multisets"
# (i.e. unordered collections in which objects can repeat)
# returns whether the first one contains the second one
def contains(bigger: collections.Counter, smaller: collections.Counter) -> bool:
	return all(bigger[key] >= smaller[key] for key in smaller)

# utility method: produces collections.Counter of the elements of an iterable
# TODO: isn't this just collctions.Counter(xs) ?
#def iter_to_counter(xs: list) -> collections.Counter:
#	to_return = collections.Counter()
#	for x in xs:
#		to_return[x] += 1
#	return to_return

# trying to fix the contains_relation code to include non-monomial relations
# implemented as follows:
# relations is a dict from strings ("input data") to polynomials ("output data")
# if your monomial contains the string, factor it out and then replace it

# input a monomial *that you know contains the input data of the relation*
# factor out the input data
# old is a monomial
# output_data, and what is returned, are polynomials
def apply_relation(old, input_data: str, output_data):
	# first, factor out the "input data"
	factored_out = collections.Counter(old.terms) - collections.Counter(input_data)
	# then, convert this into a string
	factored_out_as_str = ''.join(factored_out.elements())
	# then, turn it into a polynomial
	# and multiply it by the output data
	return str_to_poly(factored_out_as_str) * output_data


# NOTE: the zero monomial is represented by a special case, where terms == None
# this results in annoying special cases but I don't know what the alternative is
# The 1 monomial is represented by an empty list
class monomial:
	terms = list()

	def __init__(self, terms: list):
		if terms is not None:
			terms.sort()
		self.terms = terms
	
	# if we've specified products that are equal to zero, this tests whether they're in this monomial
	# if so, this monomial is zero
	# TODO: update so relns is a dict
	#def contains_relation(self):
	#	unordered_terms = collections.Counter(self.terms)
	#	return any(contains(unordered_terms, collections.Counter(rel)) for rel in relations)

	# if I contain a relation, apply it and return the result (which is a polynomial)
	# to properly account for relations, this recurses.
	# NOTE: if you have an infinite loop, check that it isn't a subtle bug in the recursion here
	# here is why I think this doesn't stack overflow, though: when you apply a relation once, if
	# you specified the relations without being deliberately dumb, you do not need to apply it again
	# so there should be at most n layers, where there are n relations.
	# TODO: at present this is not working in some important cases, chiefly of the form where you want
	# to apply vw -> stuff and have v^2 w^2.
	def update_relations(self):
		unordered_terms = collections.Counter(self.terms)
		for key, val in coh_info["relations"].items():
			if contains(unordered_terms, collections.Counter(key)):
				# NOTE: when I switched to config file, `val` became a list of strings. Before that, it was a
				# polynomial. Also: I don't know if this handles relations of the form x = 0 well (as opposed
				# to x = y)
				to_return = apply_relation(self, key, list_to_poly(val))
				to_return.update_relations()
				return to_return
		# otherwise, no relations
		return polynomial({self})
	
	def zero():
		return monomial(None)
	
	def is_0(self):
		return self.terms is None #or self.contains_relation()
	def is_1(self):
		return self.terms is not None and len(self.terms) == 0
	
	def __str__(self):
		if self.is_0(): return '0'

		# not needed since we sort on init
		#self.terms.sort()
		c = collections.Counter(self.terms)
		return ''.join(prettyprint_exponents(key, value) for key, value in c.items())
	
	def __eq__(self, other):
		if self.is_0():
			return other.is_0()
		else:
			return collections.Counter(self.terms) == collections.Counter(other.terms)

	# so that we can use this in counters
	def __hash__(self):
		# (NOTE: this might behave poorly if there are a^0 terms in a monomial
		return hash(str(self))
	
	def __add__(self, other):
		if self.is_0(): return other
		elif other.is_0(): return self
		else:
			return polynomial([self, other])

	def __mul__(self, other):
		if self.is_0() or other.is_0():
			return zero()
		else:
			return monomial(self.terms + other.terms)
	
# implemented as a list of monomials.
# NOTE: here, polynomial.zero() is a polynomial with an empty list of terms
# and 1 is represented by a polynomial single monomial term, which is 1
class polynomial:
	terms = list()

	# removes any instances of monomial.zero() from the list of terms
	# this and _trim_mod_2 are preprocessing methods
	def _remove_0s(self):
		self.terms = list(filter(lambda m: not m.is_0(), self.terms))

	# remove anything appearing twice
	# TODO: test this method
	def _trim_mod_2(self):
		c = collections.Counter(self.terms)
		self.terms = list(key for key in c if c[key] % 2 != 0)

	def __init__(self, list_of_monomials):
		self.terms = list_of_monomials
		self._remove_0s()
		self._trim_mod_2()
	
	def zero():
		return polynomial(list())
	
	def __str__(self):
		if len(self.terms) == 0:
			return '0'
		else:
			return ' + '.join(map(str, self.terms))
	
	def __hash__(self):
		return hash(str(self))

	def __eq__(self, other):
		return collections.Counter(self.terms) == collections.Counter(other.terms)

	def __add__(self, other):
		# TODO: self and othr should not be None here
		return polynomial(self.terms + other.terms)
	
	def __mul__(self, other):
		return polynomial(m*n for m in self.terms for n in other.terms)

	# iterate through all monomials and all relations. if a monomial contains a relation, update it using
	# that relation. This returns a list of polynomials; add them together.
	def update_relations(self):
		#print(self)
		new_self = sum(map(lambda x: x.update_relations(), self.terms), polynomial.zero())
		self.terms = new_self.terms
		self._remove_0s()
		self._trim_mod_2()

# take a list of lists, return a single list of all elements
#def flatten(ls: list):
#	return list(item for sublist in ls for item in sublist)

# Steenrod squares are implemented at three levels:
#	1. Sq^i of a polynomial (sum of Sq^i of its terms)
#	2. Sq^i of a monomial (use Cartan formula)
#	3. Sq^i of a single variable, which I'll specify explicitly.
# NOTE: sum() assumes the initial value is 0 unless you specify an initial value,
# which is why that polynomial.zero() is there
def sq1_poly(p: polynomial) -> polynomial:
	to_return = sum((sq1_mono(term) for term in p.terms), polynomial.zero())
	to_return.update_relations()
	return to_return
def sq2_poly(p: polynomial) -> polynomial:
	to_return = sum((sq2_mono(term) for term in p.terms), polynomial.zero())
	to_return.update_relations()
	return to_return

def sq3_poly(p: polynomial) -> polynomial:
	to_return = sum((sq3_mono(term) for term in p.terms), polynomial.zero())
	to_return.update_relations()
	return to_return
def sq4_poly(p: polynomial) -> polynomial:
	to_return = sum((sq4_mono(term) for term in p.terms), polynomial.zero())
	to_return.update_relations()
	return to_return


def sq1_mono(m: monomial) -> polynomial:
	if m.is_0():
		print("\033[34mError: 0 passed to Sq^1!\033[0m", file=sys.stderr)
		exit(-2)

	if m.is_1(): return polynomial.zero()
	else:
		first = m.terms[0] # just a string
		rest = monomial(m.terms[1:]) # a monomial
		# both here and in Sq^2, we have to convert things from strings/monomials to
		# monomials/polynomials so that + and * always operate on two things of the same
		# type. Explicit is better than implicit, I guess
		return (sq1_var(first) * polynomial([rest]) +
				polynomial([monomial([first])]) * sq1_mono(rest))

# again, this is just the Cartan formula
def sq2_mono(m: monomial) -> polynomial:
	if m.is_0():
		print("\033[34mError: 0 passed to Sq^2!\033[0m", file=sys.stderr)
		exit(-2)

	if m.is_1(): return polynomial.zero()
	else:
		first = m.terms[0] # just a string
		rest = monomial(m.terms[1:]) # a monomial
		return (sq2_var(first) * polynomial([rest]) +
			    sq1_var(first) * sq1_mono(rest) +
			    polynomial([monomial([first])]) * sq2_mono(rest))

def sq3_mono(m: monomial) -> polynomial:
	if m.is_0():
		print("\033[34mError: 0 passed to Sq^3!\033[0m", file=sys.stderr)
		exit(-2)
	
	return sq1_poly(sq2_mono(m))

def sq4_mono(m: monomial) -> polynomial:
	if m.is_0():
		print("\033[34mError: 0 passed to Sq^4!\033[0m", file=sys.stderr)
		exit(-2)

	if m.is_1(): return polynomial.zero()
	else:
		first = m.terms[0] # just a string
		rest = monomial(m.terms[1:]) # a monomial
		return (sq4_var(first) * polynomial([rest]) +
			    sq3_var(first) * sq1_mono(rest) +
			    sq2_var(first) * sq2_mono(rest) +
			    sq1_var(first) * sq3_mono(rest) +
			    polynomial([monomial([first])]) * sq4_mono(rest))


# computes Sq^2Sq^2Sq^2 of a polynomial, returns whether that is zero or not.
# I intend to use this to quickly identify A(1) summands in a larger/messier A(1)-module.
def generates_an_A1(p: polynomial) -> (bool, polynomial):
	to_return = sq2_poly(sq2_poly(sq2_poly(p)))
	return to_return != polynomial.zero(), to_return

# given a list of polynomials, checks which (if any) generate an A(1). might miss some A(1)s,
# e.g. because you didn't specify enough inputs. For example, I plan on just testing monomials
# for now.
def identify_A1_summands(ps: list):
	for p in ps:
		gens_yes, result = generates_an_A1(p)
		if gens_yes:
			print('%s generates an A(1) summand ending in %s' % (str(p), str(result)))
	
def generates_a_Joker(p: polynomial) -> (bool, polynomial):
	first_check = sq2_poly(sq2_poly(p))
	second_check = sq1_poly(sq2_poly(p))
	to_return = first_check != polynomial.zero() and second_check == polynomial.zero()
	return to_return, first_check

def identify_Joker_summands(ps: list):
	for p in ps:
		gens_yes, result = generates_a_Joker(p)
		if gens_yes:
			print('%s generates a Joker ending in %s' % (str(p), str(result)))


# quick wrapper to convert some characters into a polynomial with a single term
def str_to_poly(s: str) -> polynomial:
	return polynomial([monomial(list(s))])


# another useful utility: given a starting set of summands, determine the submodule generated by
# those summands. Input is a set of polynomials
def compute_span_of(ps: set, sq_four = False):
	while len(ps) != 0:
		p = ps.pop()
		sq1_of_p = sq1_poly(p)
		sq2_of_p = sq2_poly(p)
		print('Sq1(%s) = %s' % (str(p), str(sq1_of_p)))
		print('Sq2(%s) = %s' % (str(p), str(sq2_of_p)))
		if sq1_of_p != polynomial.zero():
			ps.add(sq1_of_p)
		if sq2_of_p != polynomial.zero():
			ps.add(sq2_of_p)
		if sq_four:
			sq4_of_p = sq4_poly(p)
			print('Sq4(%s) = %s' % (str(p), str(sq4_of_p)))
			if sq4_of_p != polynomial.zero():
				ps.add(sq4_of_p)


# it would be nice to specify the ring structure in general, but that seems a bit complicated
# so instead, I'll specify a few monomials which are 0. This works for some cohomology rings
# (e.g. BS_4) but not everyone (e.g. BD_{2n}, n = 0 mod 4
# specify as a dict from strings  to polynomials
# it helps to have a relation that is guaranteed to only need finitely many applications
#relations = {
#	'uuu': str_to_poly('vw') + str_to_poly('vv') + str_to_poly('ww')
#	'xy': str_to_poly('yy')
#}
#{'abe': str_to_poly('ag') + str_to_poly('bf') + str_to_poly('ce')}
#{'uuu': str_to_poly('vw') + str_to_poly('vv') + str_to_poly('ww')}

#sq1_dict = {
#	'a': 'b',
#	'b': '0',
#	'U': '0'
#}
#
#sq2_dict = {
#	'a': 'aa',
#	'b': 'ab',
#	'U': 'Ua'
#}

#sq4_dict = {
#	'a': '0',
#	'b': '0',
#	'c': 'd',
#	'd': '0',
#	'e': 'ee',
#	'f': 'h',
#	'g': 'i',
#	'h': '0',
#	'i': '0'
#}

# convert xs (a list of strings) to a polynomial,
# e.g. the way Sqi is specified in the config files
def list_to_poly(xs: list):
	if xs == "0":
		return polynomial.zero()
	else:
		return functools.reduce(lambda x, y: x+y, map(str_to_poly, xs), polynomial.zero())

# helper method for Steenrod squares key = 'Sq1' or 'Sq2' or 'Sq4'
def square_var(ch: str, key: str):
	if ch in coh_info["classes"]:
		ans = coh_info["classes"][ch][key]
		if ans == '0':
			return polynomial.zero()
		else:
			return list_to_poly(ans)
	else:
		print("\033[33mError: %s of variable '%s' not implemented\033[0m" % (key, ch), file=sys.stderr)
		exit(-1)


# at this point, just specify them manually.
# note: this is a lot of nested dict accesses
# if I find this slows the program down, I should extract this info from the config file once,
# at the beginning, and then reference the smaller discts later
def sq1_var(ch: str):
	return square_var(ch, 'Sq1')
#	if ch in coh_info["classes"]:
#		ans = coh_info["classes"][ch]['Sq1']
#		if ans == '0':
#			return polynomial.zero()
#		else:
#			return list_to_poly(ans)
#	else:
#		print("\033[33mError: Sq^1 of variable '%s' not implemented\033[0m" % ch, file=sys.stderr)
#		exit(-1)

def sq2_var(ch: str):
	return square_var(ch, 'Sq2')	
#	if ch in coh_info["classes"]:
#		ans = coh_info["classes"][ch]['Sq2']
#		if ans == '0':
#			return polynomial.zero()
#		else:
#			return list_to_poly(ans)
#	if ch in sq2_dict:
#		if sq2_dict[ch] == '0':
#			return polynomial.zero()
#		else:
#			return str_to_poly(sq2_dict[ch])
#	elif ch == 'U':
#		return str_to_poly('Uw') + str_to_poly('Uxx')
#	else:
#		print("\033[33mError: Sq^2 of variable '%s' not implemented\033[0m" % ch, file=sys.stderr)
#		exit(-1)

def sq3_var(ch: str):
	return sq1_poly(sq2_var(ch))

def sq4_var(ch: str):
	return square_var(ch, 'Sq4')
#	if ch in coh_info["classes"]:
#		ans = coh_info["classes"][ch]['Sq4']
#		if ans == '0':
#			return polynomial.zero()
#		else:
#			return list_to_poly(ans)
#	else:
#		print("\033[33mError: Sq^4 of variable '%s' not implemented\033[0m" % ch, file=sys.stderr)
#		exit(-1)

# a few try/catch clauses to make sure things that I am assuming are present in the config file are actually there
# this is very much not comprehensive.
def validate(cohomology_info: dict):
	if "relations" not in cohomology_info:
			print("\033[33mError: 'relations' not specified in config file. Even if the cohomology ring is free, please supply an empty dict\033[0m", file=sys.stderr)
	# more to come

def the_deg(x: str) -> int:
	if x == 'a': return 1
	if x == 'b' or x == 'c': return 2
	if x == 'd' or x == 'e': return 3
	else:
		print('This wasn\'t supposed to happen, x = %s' % x)
	
def sum_deg(x: str) -> int:
	return sum(map(the_deg, x))

# true if total degree 9 or less
# quick n dirty helper fn
def not_too_big(s: str) -> bool:
	return sum_deg(s) <= 9

def main():
	if len(sys.argv) < 2:
		print("\033[33mError: please specify cohomology via config file.\033[0m")
		exit(-1)
	
	filename = sys.argv[1]
	global coh_info
	try:
		with open(filename, 'r', encoding='utf-8') as infile:
			coh_info = json.load(infile)
	except IOError:
		print("\033[33mError: config file %s cannot be loaded or is invalid JSON.\033[0m" % filename, file=sys.stderr)
		# TODO: one day I'll validate the input
		exit(-1)
	
	# ok, we've now loaded the cohomology information, let's validate it
	validate(coh_info)

	too_many_tuples = list(list(itertools.combinations_with_replacement('abcde', k)) for k in range(12))
	flattened_tuples = list(p for q in too_many_tuples for p in q)
	
	how_many = collections.Counter(sum_deg(x) for x in flattened_tuples)
#	print(how_many)

	all_tuples = sorted(list(filter(not_too_big, flattened_tuples)), key=sum_deg)

	to_test = ['U' + ''.join(mon) for mon in all_tuples]

#	to_test = ['U' + ''.join(mon) for k in range(10) for mon in list.sort(list(filter(not_too_big, itertools.combinations_with_replacement('abcde', k))), key=lambda x: sum(map(the_deg, ''.join(x))))]


	#for i in range(1, 10):
	#	print('Sq^1(x^%d) = %s' % (i, str(sq1_mono(monomial(['x']*i)))))
	#	print('Sq^2(x^%d) = %s' % (i, str(sq2_mono(monomial(['x']*i)))))
#	to_test = ['U' + 'a'*(n-i) + 'b'*i for n in range(0, 6) for i in range(n+1)]
#	to_test = ["U",
#			   "Ux", "Uy",
#			   "Uw", "Uxx", "Uyy",
#			   "Uxxx", "Uyyy", "Uwx", "Uwy",
#			   "Uxxxx", "Uyyyy", "Uwxx", "Uwyy", "Uww",
#			   "Uxxxxx", "Uyyyyy", "Uwxxx", "Uwyyy", "Uwwx", "Uwwy",
#			   "Uxxxxxx", "Uyyyyyy", "Uwxxxx", "Uwyyyy", "Uwwxx", "Uwwyy", "Uwww",
#			   "Uxxxxxxx", "Uyyyyyyy", "Uwxxxxx", "Uwyyyyy", "Uwwxxx", "Uwwyyy", "Uwwwx", "Uwwwy"]

	#to_test = ['U', 'Ux', 'Uy', 'Uz', 'Uw', 'Uxx', 'Uyy', 'Uxz', 'Uyz', 'Uzz',
	#	'Uxxx', 'Uyyy', 'Uzzz', 'Uxxz', 'Uxzz', 'Uyyz', 'Uyzz', 'Uwx', 'Uwy', 'Uwz']
	#to_test = ['Uaaabb', 'Uabbbb']

#	identify_A1_summands(map(str_to_poly, to_test))
	#print(generates_an_A1(str_to_poly('Uaaaaabbbb') + str_to_poly('Uabbbbbbbb')
	#))
	#identify_Joker_summands(map(str_to_poly, to_test))
	#print(generates_an_A1(str_to_poly('Uaaaaaaabb') + str_to_poly('Uaaaaabbbb'))[0])
	
	compute_span_of({str_to_poly('Uaaa')
	})


	#for term in to_test:
	#	m = monomial(list(term))
	#	print('Sq^1(%s) = %s' % (str(m), str(sq1_mono(m))))
	#	print('Sq^2(%s) = %s' % (str(m), str(sq2_mono(m))))
	
if __name__ == '__main__':
	main()



# For D_{2n}, n = 0 mod 4

