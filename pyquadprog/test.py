from PyQuadProg import PyQuadProg
from numpy import array

G = array([[11., 0], [0, 12]])
CE = array([[1.], [1]])
CI = array([[1.], [1]])

g0 = array([10., 20])
ce0 = array([-5.])
ci0 = array([10])


qp = PyQuadProg(G, g0, CE, ce0, CI, ci0)

print array(qp.x.getArray())
