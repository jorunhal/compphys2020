from numpy import array, zeros, linspace, exp
import matplotlib.pyplot as plt
from time import perf_counter_ns

def solve_spec_tridag(g):
	# This function takes a vector g and solves the system of equations Av = b
	# for v, where A is a tridiagonal matrix with (-1, 2, -1) along the dia-
	# gonal

	n = len(g)

	# Compute values of r and s

	r = linspace(1.0, n, n)
	r = r/(1.0 + r)
	s = zeros(n - 1)
	s[0] = g[0]/2.0
	for i in range(n - 2):
		s[i + 1] = r[i + 1]*(s[i] + g[i + 1])

	# Backwards substitution to solve for the components of v

	v = zeros(n)
	v[n - 1] = r[n - 1]*(g[n - 1] + s[n - 2])
	for j in range(1, n):
		v[n - j - 1] = r[n - j - 1]*v[n - j] + s[n - j - 1]
	return v

def f(x):
	# Function for computing the source term values
	return 100.0*exp(-10.0*x)

def F(x):
	# Function to computing the values of the analytic solution
	return 1.0 - (1.0 - exp(-10.0))*x - exp(-10.0*x)

# Performing the algorithm for varying numbers of grid points, keeping track
# of time usage

p_ = array([ p for p in range(1, 7) ])	# Powers of 10
times_spec = zeros(len(p_))		# Times for specific tridiagonal method

for i in range(len(p_)):		# Loop through powers of 10
	n = 10**p_[i]
	x = linspace(0.0, 1.0, n + 2)
	h = 1.0/(n + 1)
	g = h**2*f(x[1:-1])		# Source term

	# Solve with specialized algorithm

	start = perf_counter_ns()/10**9
	v2 = solve_spec_tridag(g)
	end = perf_counter_ns()/10**9
	times_spec[i] = end - start

# Print CPU times

print("log10(h)\tTime (s)")
for i in range(len(p_)):
	print(str(p_[i]) + "\t\t" + str(times_spec[i]))

# log10(h)        Time (s)
# 1               0.00012579999999995373
# 2               0.0003671000000000646
# 3               0.0032730999999999177
# 4               0.03334119999999996
# 5               0.3499977999999999
# 6               3.5954535