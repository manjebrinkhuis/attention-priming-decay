import numpy as np


def add_regressors(matrix, regressors):
    """ Add a regressor to a design matrix """
    M = np.array(matrix)
    return np.vstack((M, regressors))


def array2dummies(arr):
    """ Returns matrix """
    values = sorted(list(set(arr)))
    values = np.array([[v] for v in values])
    return (arr == values)


def n_back_array(arr, n_back):
    """ Returns boolean array comparing each index to the index - n_back """
    arr = np.array(arr)
    pre = np.zeros(n_back)
    rep = arr[n_back:] == arr[:len(arr)-n_back]
    return np.concatenate((pre, rep))


def n_back_series(arr, max_back):
    """ Returns matrix containing 1-back to max_back regressors. """
    M = np.zeros((0, len(arr)))
    for i in range(max_back):
        rep = n_back_array(arr, i+1)
        M = add_regressors(M, rep)
    return M


# A set of functions to perform a general linear model
# analysis on any type of data (so not just fMRI).
def glm(Y, X):
    """ GLM - Returns b-weights for each row in X with timeseries Y """
    X = np.matrix(X).T
    Y = np.matrix(Y)
    B = np.linalg.pinv(X.T*X)*X.T*Y.T
    return np.squeeze(np.array(B))


def contrast(Y, B, X, contrasts=[]):
    """ """

    X = np.matrix(X).T
    N, H = X.shape
    XX = np.linalg.pinv(X.T*X)
    predicted = np.sum(X.dot(B), axis=0)

    SSe = (np.array(Y - predicted)**2).sum()
    MSe = SSe / (N-H)
    C = np.array(contrasts)
    explained = np.sum(B*C)
    unexplained = MSe*(C.T.dot(XX).dot(C))
    t = explained / unexplained[0, 0]**.5
    unexplained

    return t
