import numpy as np

n = 10
M = np.random.rand(n,n)
Mt = np.matrix.transpose(M)

coherence = np.dot(M,Mt)

c = {x: {y: coherence[x,y] for y in range(n)} for x in range(n)}
clist = np.array(coherence).reshape(-1,).tolist()

max()




