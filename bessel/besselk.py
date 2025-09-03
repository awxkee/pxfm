from mpmath import mp
import sys

x = float(sys.argv[1])
mp.prec = 110
# print(mp.besseli(2, x))
print(mp.exp(mp.mpf(x)*mp.mpf(x)) * mp.erfc(mp.mpf(x)))
# print(mp.erfinv(x))