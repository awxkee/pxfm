from mpmath import mp
import sys

x = float(sys.argv[1])
mp.prec = 100
# print(mp.besseli(2, x))
print(mp.erfinv(mp.mpf('1') - mp.mpf(x)))
# print(mp.erfinv(x))