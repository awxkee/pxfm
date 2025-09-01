from mpmath import mp
import sys

x = float(sys.argv[1])
mp.prec = 120
print(mp.erfinv(x))