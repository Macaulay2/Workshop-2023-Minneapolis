load "./_init.m2"

b = matrix(QQ,{{1,2},{2,0}})
beta = gwClass(b)

g = matrix(QQ,{{5,6,7},{6,0,1},{7,1,0}})
gamma = gwClass(g)


-- For checking upper triangularity

b1 = matrix(QQ,{{5,6,1},{6,2,0},{1,0,0}})
b2 = matrix(GF(17), {{4,3,1,7},{3,2,5,0},{1,5,0,0},{7,0,0,0}})
g = matrix(QQ,{{5,6,7},{6,0,1},{7,1,0}})

