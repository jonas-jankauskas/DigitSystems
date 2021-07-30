"""example.sage: SAGE script  for the examples of  Section 7.1 the paper Rational matrix digit systems by J. Jankauskas and J. M. Thuswaldner """

__author__ = "Jonas Jankauskas"
__copyright__ = "Copyright 2021, Jonas Jankauskas"
__license__ = "LGPLv3"
__version__ = "1.0"
__maintainer__ = "Jonas Jankauskas"
__email__ = "jonas.jankauskas@gmail.com"
__status__ = "Production"

reset()

load('digsys.sage')

M = matrix(QQ, 2, 2, [[3/5, -4/5],[4/5, 3/5]])

print('\nM=', M, '\n')

#Automatically build the digit system with a given matrix M with default hulls and smaller remainder vectors 
rds = RotationDS(M, rem='reduced')

#Report results
print('\nResidues:\n', rds.res, '\n')
print('\nOld hulls:\n', rds.hull, '\n')
#{0: (0, 0), 1: (1, 0), 2: (2, 0), 3: (-2, 0), 4: (-1, 0)}

#Rebuild the DS using a smaller manually chosen hulls
newhull = {0: [(-2, -1), (-1, 2), (3, -1)], 1: [(0, 0), (2, 1), (3, -1)], 2: [(0, 0), (2, 1), (3, -1)], 3: [(0, 0), (-2, -1), (-3, 1)], 4: [(0, 0), (-2, -1), (-3, 1)]}
rds.hull = {r: [vector(ZZ, v) for v in newhull[r]] for r in rds.res.keys()}

rds._build_digits_()

#Report the results
print('\nNew hulls:\n', rds.hull, '\n')
print('\nNew digits:\n', rds.alldigs, '\n')
print('\nNew attractor digits:\n', rds.tidy_attractor, '\n')

#Report expansions of (-6, 7) and (0, 0) in new DS
u = vector(ZZ, [6, -7])
v = vector(ZZ, [0, 0])

print('\nRadix expansion of ', u, 'in the new DS is:\n', rds.expand(u)[1], '\n')
print('\nRadix expansion of ', v, 'in the new DS is:\n', rds.expand(v)[1], '\n')

#Calculate the expansion of (-3, 0)+ A*(0, 2)+A^2*(-1, 2)
v1 = vector(ZZ, 2, [-3, 0])
v2 = vector(ZZ, 2, [0, 2])
v3 = vector(ZZ, 2, [-1, 2])

carry, w1 = rds.F(v1)
carry, w2 = rds.F(v2+carry)
carry, w3 = rds.F(v3+carry)
diglist = [w1, w2, w3]+rds.expand(carry)[1]

print("Expansion of ", v1, "+A*", v2, "+A^2*", v3, " = ", rds.assemble(diglist), "=", diglist)
