import logging
from scipy.io.matlab import loadmat
from scipy.sparse import lil_matrix
from rescal import rescal_als

# Set logging to INFO to see RESCAL information
logging.basicConfig(level=logging.INFO)

# Load Matlab data and convert it to dense tensor format
T = loadmat('data/alyawarradata.mat')['Rs']
X = [lil_matrix(T[:, :, k]) for k in range(T.shape[2])]

# Decompose tensor using RESCAL-ALS
A, R, fit, itr, exectimes = rescal_als(X, 100, init='nvecs', lambda_A=10, lambda_R=10)
