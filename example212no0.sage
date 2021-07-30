"""example212no0.sage: SAGE script  for the example of rotational digit system in  Section 7.2 of the paper Rational matrix digit systems by J. Jankauskas and J. M. Thuswaldner """

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

rds = RotationDS(A, 'reduced')

newhull = {0: [(-2, 0), (2, 1), (2, -1)], 1: [(0, 0), (2, 1), (2, -1)]}
rds.hull = {r: [vector(ZZ, v) for v in newhull[r]] for r in rds.res.keys()}
rds._build_digits_()
