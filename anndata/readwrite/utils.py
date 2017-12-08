import numpy as np

def save_sparse_csr(X, key='X'):
    from scipy.sparse.csr import csr_matrix
    X = csr_matrix(X)
    key_csr = key + '_csr'
    return {key_csr + '_data': X.data,
            key_csr + '_indices': X.indices,
            key_csr + '_indptr': X.indptr,
            key_csr + '_shape': np.array(X.shape)}


def load_sparse_csr(d, key='X'):
    from scipy.sparse.csr import csr_matrix
    key_csr = key + '_csr'
    d[key] = csr_matrix((d[key_csr + '_data'],
                         d[key_csr + '_indices'],
                         d[key_csr + '_indptr']),
                        shape=d[key_csr + '_shape'])
    del d[key_csr + '_data']
    del d[key_csr + '_indices']
    del d[key_csr + '_indptr']
    del d[key_csr + '_shape']
    return d
