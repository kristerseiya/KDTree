
import pykdtree
import pypcdio
import numpy as np
import linecache
import os
import tracemalloc
from sklearn.neighbors import KDTree as skKDTree


def display_top(snapshot, key_type='lineno', limit=10):
    snapshot = snapshot.filter_traces((
        tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
        tracemalloc.Filter(False, "<unknown>"),
    ))
    top_stats = snapshot.statistics(key_type)

    print("Top %s lines" % limit)
    for index, stat in enumerate(top_stats[:limit], 1):
        frame = stat.traceback[0]
        print("#%s: %s:%s: %.1f KiB"
              % (index, frame.filename, frame.lineno, stat.size / 1024))
        line = linecache.getline(frame.filename, frame.lineno).strip()
        if line:
            print('    %s' % line)

    other = top_stats[limit:]
    if other:
        size = sum(stat.size for stat in other)
        print("%s other: %.1f KiB" % (len(other), size / 1024))
    total = sum(stat.size for stat in top_stats)
    print("Total allocated size: %.1f KiB" % (total / 1024))


pcd = pypcdio.pcdread("/Users/Krister/data/PointClouds/dual_clutch_gearbox.ply")
# nan_idx = np.argwhere(np.isnan(pcd.points))
# x = np.delete(pcd.points, nan_idx[:,0], 0)
#kdt = pykdtree.KDTree()
#kdt.assign(pcd.points)
skkdt = skKDTree(pcd.points)
tracemalloc.start()
#kdt.searchKNN(pcd.points[0],10000)
skkdt.query([pcd.points[0]],10000)
snapshot = tracemalloc.take_snapshot()
display_top(snapshot)

# result = np.concatenate([x[idx,:], np.expand_dims(dist,1)], axis=1)
# print(result)