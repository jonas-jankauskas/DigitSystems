reset()

attach('digsys.sage')

M = matrix(QQ, 2, 2, [[0, -1],[1, -1/2]])
N = matrix(ZZ, 2, 2, [[0, 1],[0, 0]])	

#print(M)

rds = RotationDS(M, 'normal')
tds = TwistedDS(M, M, N)

#newhull = {0: [(2, 0), (-2, 1), (-2, -1)], 1: [(0, 1), (0, -1), (2, 0)]}
#rds.hull = {r: [vector(ZZ, v) for v in newhull[r]] for r in rds.res.keys()}
load()#rds._build_digits_()

u = vector(ZZ, [1, 2])
v = vector(ZZ, [-3, 4])

w = vector(u.list()+v.list())


