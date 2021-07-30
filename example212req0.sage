"""example212req0.sage: SAGE script  for the example of twisted  digit system in  Section 7.2 of the paper Rational matrix digit systems by J. Jankauskas and J. M. Thuswaldner """

__author__ = "Jonas Jankauskas"
__copyright__ = "Copyright 2021, Jonas Jankauskas"
__license__ = "LGPLv3"
__version__ = "1.0"
__maintainer__ = "Jonas Jankauskas"
__email__ = "jonas.jankauskas@gmail.com"
__status__ = "Production"


reset()

attach('digsys.sage')

R.<x> = QQ['x']

f = x^2+1/2*x+1

A = companion_matrix(f)

rds = RotationDS(A, 'reduced', req0=True)

newhull = {0: [(-2, 0), (2, 1), (2, -1)], 1: [(0, 0), (2, 1), (2, -1)]}
rds.hull = {r: [vector(ZZ, v) for v in newhull[r]] for r in rds.res.keys()}
rds._build_digits_()

N = matrix(ZZ, 2, 2, [[0, 1],[0, 0]])	

tds = TwistedDS(rds, rds, N)

v = vector(ZZ, [1, 2, -3, 4])

print("\n Expand vector ", v, "=", tds.expand(v)[1])
