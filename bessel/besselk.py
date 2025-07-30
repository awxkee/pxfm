from mpmath import mp
import sys

x = float(sys.argv[1])
mp.prec = 100
print(mp.besselk(0, x))