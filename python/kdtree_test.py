
import pykdtree
import pypcdio
import numpy as np


pcd = pypcdio.pcdread("/Users/Krister/data/PointClouds/iowaClay_high_res.xyzm")
nan_idx = np.argwhere(np.isnan(pcd.points))
x = np.delete(pcd.points, nan_idx[:,0], 0)
kdt = pykdtree.KDTree()
kdt.assign(x)
idx, dist = kdt.searchKNN(x[0],20)
result = np.concatenate([x[idx,:], np.expand_dims(dist,1)], axis=1)
print(result)