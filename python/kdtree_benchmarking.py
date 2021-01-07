from pypcdio import pcdread
import numpy as np
from sklearn.neighbors import KDTree as skKDTree
from scipy.spatial import KDTree as scKDTree
from pykdtree import KDTree as myKDTree
import timeit

def sklearn_kdtree_test(data):
    kdt = skKDTree(data)
    return kdt

def my_kdtree_test(data):
    kdt = myKDTree()
    kdt.assign(data)
    return kdt

pcd1 = pcdread("/Users/Krister/data/PointClouds/dual_clutch_gearbox.ply")
pcd2 = pcdread("/Users/Krister/data/PointClouds/human_skeleton.ply")
pcd3 = pcdread("/Users/Krister/data/PointClouds/gom.stl")

data = np.concatenate([pcd1.points, pcd2.points, pcd3.points], 0)
