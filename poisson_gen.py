from numpy import array, zeros, linspace, exp, log10, abs, max
import matplotlib.pyplot as plt
from time import perf_counter_ns

def solve_tridiag(a, b, c, g):
	# This function takes the diagonal and off-diagonal elements of a tri-
	# diagonal matrix A as vectors, and a vector f, and solves the set of 
	# linear equations
	# A*v = g for v.

	n = len(g)
	if (len(a) != n - 1) or (len(b) != n) or (len(c) != n - 1):
		return zeros(k)

	# Compute values r and s

	r, s = zeros((2, n - 1))
	r[0] = -c[0]/b[0]
	s[0] = g[0]/b[0]
	for i in range(n - 2):
		r[i + 1] = -c[i + 1]/(a[i]*r[i] + b[i + 1])
		s[i + 1] = (g[i + 1] - a[i]*s[i])/(a[i]*r[i] + b[i + 1])

	# Backwards substitution to compute components of and return of array v

	v = zeros(n)
	v[n - 1] = (g[n - 1] - a[n - 2]*s[n - 2])/(a[n - 2]*r[n - 2] + b[n - 1])
	for j in range(1, n):
		v[n - j - 1] = r[n - j - 1]*v[n - j] + s[n - j - 1]
	return v

def f(x):
	# Function for computing the source term values
	return 100.0*exp(-10.0*x)

def F(x):
	# Function to computing the values of the analytic solution
	return 1.0 - (1.0 - exp(-10.0))*x - exp(-10.0*x)


# Perform algorithms in order, keeping track of time usage and computing 
# relative errors

p_ = array([ p for p in range(1, 7) ])	# Powers of 10
times_gen = zeros(len(p_))		# Times for general tridiagonal method

errs = zeros(len(p_))			# Logarithm base-10 of relative errors

for i in range(len(p_)):		# Loop through powers of 10
	n = 10**p_[i]
	x = linspace(0.0, 1.0, n + 2)
	h = 1.0/(n + 1)
	g = h**2*f(x[1:-1])

	# Analytic solution

	u = F(x[1:-1])

	# Define diagonal elements and solve with general method

	a = zeros(n - 1) - 1.0
	b = zeros(n) + 2.0
	c = zeros(n - 1) - 1.0

	start = perf_counter_ns()/10**9
	v = solve_tridiag(a, b, c, g)
	end = perf_counter_ns()/10**9
	times_gen[i] = end - start

	# Plot results alongside analytic solution (for p <= 3)

	if p_[i] <= 3:
		plt.plot(x[1:-1], u)
		plt.plot(x[1:-1], v)
		plt.axis([0.0, 1.0, 0.0, 0.7])
		plt.title("Analytic vs. numerical solutions ($n = " + str(n) + "$)")
		plt.legend(("Analytic", "Numerical"))
		plt.xlabel("$x$")
		plt.ylabel("$u(x)$, $v_i$")
		plt.show()

	# Compute and plot relative error of general tridiagonal algorithm, 
	# should be the same for specific tridiagonal algorithm
	
	errs[i] = max(log10(abs((v - u)/u)))

# Plot relative error for numerical method

plt.plot(p_, errs)
plt.title("Logarithm base-10 of relative error")
plt.xlabel("$\\log_{10}(h)$")
plt.ylabel("$\\mathrm{max}\\:\\varepsilon_i$")
plt.show()

# Print CPU times

print("log10(h)\tTime (s)")
for i in range(len(p_)):
	print(str(p_[i]) + "\t\t" + str(times_gen[i]))

# log10(h)        Time (s)
# 1               9.48000000000615e-05
# 2               0.0016837000000000657
# 3               0.0179866999999998
# 4               0.1490172000000003
# 5               1.0724179999999999
# 6               10.035969200000002