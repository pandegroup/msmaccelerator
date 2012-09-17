
"""
SinghalSampler.py

This file is part of MSMAccelerator.

To Do:
-- Test eigenvector calculation for correctness
-- Confirm that pivoting from ssaCalc is correct
-- Install and test UMFPACK routine
-- Confirm sparse correctness

Written August 2012, by TJ Lane <tjlane@stanford.edu>
"""

import sys

import numpy as np
import scipy, scipy.io
from scipy.sparse import linalg as spla
from scipy.sparse.linalg.dsolve import umfpack

try:
    ctx = umfpack.UmfpackContext()
    UMFPACK_SUPPORT = True
except ImportError as e:
    UMFPACK_SUPPORT = False


DEBUG = True  # set this to do a bunch of internal checking
              # but at significant computational cost


class SinghalSampler(object):
    """
    This sampler minimizes the statistical error in one or more
    eigenmodes (values and/or vectors) in the transition matrix.

    To initialize the object, you have to pass two lists describing
    which processes you care about and their respective weights. The
    algorithm can minimize the error in the eigenvalues, or both the
    eigenvalues and the eigenvectors.

    Further, three priors are available:

    -- Jeffery's prior, alpha = 1/n 
    -- uniform prior,   alpha = 1
    -- sparse prior,    alpha = 0   (allows for sparse data structures)

    The choice of prior likely makes a small difference, but choosing the
    sparse prior should result in increases in memory efficiency but
    losses in computational speed. Use dense linear algebra unless
    absolutely necessary.
    

    Parameters
    ----------
    mode_indices : array_like, int
        the indices of the modes you want to minimize the error of
        (zero indexed)
        
    mode_weights : array_like, float
        the respective weight you want to algorithm to give to each mode,
        higher weight means it will work proportionally harder to minimize
        error in that mode
        
    prior : str
        the Dirichlet prior to employ - see above - one of 'Jefferys',
        'uniform', 'sparse'

    eigenvector_weight : float
        A float between 0.0 and 1.0 that indicates how much to weight the
        eigenvectors relative to the eigenvalues in the error analysis.
        0.0 indicates only the eigenvalues will be used, 1.0 is only the
        eigenvectors, and intermediate values are a linear combination of
        the two with the given weight.


    Notes
    -----
    Main use of this sampler is called through the sample() method.

    References
    ----------
    :: [1] Singhal, NS and Pande, VS J. Chem. Phys. 126, 244101 (2007)
    """


    def __init__(self, mode_indices, mode_weights,
                 prior='Jefferys', eigenvector_weight=0.0):

        # store and check the requested modes
        self.mode_indices = np.array( mode_indices )
        self.mode_weights = np.array( mode_weights )
        self.mode_weights /= self.mode_weights.sum()
        if not len(mode_indices) == len(mode_weights):
            raise ValueError("Each weight must represent a single eigenmode index")
        
        if (eigenvector_weight < 0.0 or eigenvector_weight > 1.0):
            raise ValueError("Parameter `eigenvector_weight` must be between"
                             "0.0 and 1.0. Passed was: %f" % eigenvector_weight)
        else:
            self.eigenvector_weight = eigenvector_weight
        if self.eigenvector_weight == 0.0:
            self.use_vectors = False
        else:
            self.use_vectors = True

        # check the prior
        priors = ['Jefferys', 'uniform', 'sparse']
        if prior not in priors:
            raise ValueError("Prior must be one of: %s" % str(priors))
        else:
            self.prior = prior


    def sample(self, counts_matrix, force_dense=False):
        """
        Retrieve the multinomial distribution for the optimal
        sampling strategy. This is the main routine that should
        be called in order to do adpative sampling.

        Parameters
        ----------
        counts_matrix : sparse matrix
            the counts for the MSM of interest
        force_dense:
            force use of dense linear algebra - may increase performace
            for small MSMs using the sparse prior

        Returns
        -------
        multinomial : array_like, float
            an array for the multinomial distribution over states

        Notes
        -----
        The algorithm tries to be smart about maintaining sparse data
        structures where possible. This is only possible when `prior`
        is set to sparse, however (when the sampler is initialized).
        If you're running into memory or speed issues, try the sparse
        prior.
        """

        # if the sparse prior is used try sparse linear alg.
        self.sparse = False
        if self.prior == 'sparse':
            if not force_dense:
                if UMFPACK_SUPPORT:
                    self.sparse = True
                else:
                    raise Exception("Could not access UMFPACK via scipy, which provides "
                                    "necessary support for sparse linear algebra. "
                                    "Try installing it -- see MSMAccelerator Docs.")

        if DEBUG: print "SinghalSampler: prior=%s, sparseLA=%s" % (self.prior, self.sparse)

        # store some internals for later
        self.C = counts_matrix
        self.N = self.C.shape[0]

        # estimate the transition probability matrix (expectation of posterior)
        self.Pbar, self.w = self._calc_Pbar(counts_matrix)

        # now cacluate the eigenvalues/vectors (left) of interest
        if self.sparse:
            v, y = spla.eigs(self.Pbar.T, k=np.max(self.mode_indices)+1, which='LR')
        else: # if dense
            v, y = np.linalg.eig(self.Pbar.T)

        # if eigenvalues are imag, then we didn't have a good counts matrix
        if np.sum( np.abs( np.imag(v) )) > 10.**-6:
            print >> sys.stderr, "WARNING in SinghalSampler:"
            print >> sys.stderr, "\tThe eigenvalues of the transition matrix are imaginary!"
            print >> sys.stderr, "\tThis means the model did not satisfy detailed balance."

        v = np.real(v)
        perm = np.argsort(v)[::-1]
        self.values = v[perm][self.mode_indices]
        self.vectors = np.real( y[:,perm] )[:,self.mode_indices]

        del v # we have stored these in self now, free them up
        del y

        # iterate over all of the eigenmodes and determine their error
        multinomial = np.zeros( self.N )
        for index, mode in enumerate(self.mode_indices):
            if DEBUG:
                print "SinghalSampler: Evaluating error in mode: %d, %f" % \
                      (mode, self.values[index])

            # calculate the partial derivatives d(lambda) / d(T_{ij})
            Abar = self.Pbar - scipy.sparse.eye(self.N, self.N) * self.values[index]
            s, S = self._eigmode_partials(Abar, index)

            # calculate the error factors
            sample_weight = np.zeros( self.N )
            for i in range( self.N ): # eq. (25) in [1]

                if self.sparse:
                    pi = self.Pbar.getrow(i).toarray().flatten()
                    K = scipy.sparse.dia_matrix((pi, 0), shape=self.Pbar.shape) - \
                        self._sparse_outer_product(pi, pi)
                    qi = (s[i,:] * (K * s[i,:].T)).toarray()[0,0]
                else:
                    pi = self.Pbar[i,:]
                    K = np.diag(pi)-np.outer(pi, pi)
                    qi = np.dot( np.dot( s[i,:].T, K ), s[i,:] )
                    
                if qi < 0.0:   # There is inherent error in the approximations used to 
                    qi = 0.0   # calculate the partials -- we may get a prediction 
                    if DEBUG:  # of negative variance decrease. It should be small.
                        print "SinghalSampler: Warning, q_i value is negative: %f" % qi

                norm = (self.w[i]**2 + 3.0*self.w[i] + 2.0)
                sample_weight[i] = qi / norm

            multinomial += (1.0 - self.eigenvector_weight) * \
                           sample_weight * self.mode_weights[index]

            # calculate the eigenvector covariance matrix (eq. 24 in [1])
            if self.use_vectors:
                assert S != None

                # let Ri be the analogs of qi (but matrix instead of scalar)
                Ri = np.dot( np.dot( S[i,:,:], K ), S[i,:,:].T )

                # simply treat all eigenvector components equallly for the
                # purposes of determining the error
                sample_weight = np.sum( Ri / norm, axis=1 )

                if len(np.where( sample_weight < 0.0 )) != 0:   
                    sample_weight[ sample_weight < 0.0 ] = 0.0
                    if DEBUG:
                        print len(np.where( sample_weight < 0.0 ))
                        print("SinghalSampler: Warning, %d elements of 'sample_weight'"
                              " in the eigenvector calculation were less than"
                              " zero" % len(np.where( sample_weight < 0.0 )) )
                            
                multinomial += self.eigenvector_weight * \
                               sample_weight * self.mode_weights[index]

        multinomial /= multinomial.sum()
        return multinomial


    def _calc_Pbar(self, C):
        """
        Returns the P-bar matrix from [1], which is the expectation
        of the Dirichlet posterior of the transition matrix.
        """

        # choose alpha values based on the prior
        if self.prior == 'Jefferys':
            alpha = np.ones(self.N) / float(self.N)
            
        elif self.prior == 'uniform':
            alpha = np.ones(self.N)
            
        elif self.prior == 'sparse': # just normalize and return
            if self.sparse:
                if not scipy.sparse.issparse(C):
                    C = scipy.sparse.lil_matrix(C)
                w = np.array( C.sum(axis=1).astype(float) ).flatten()
                Winv = scipy.sparse.dia_matrix((1.0/w, 0), shape=self.C.shape)
                P = C * Winv
                assert scipy.sparse.issparse(P)
                return P.tocsr(), w
            else:
                alpha = np.zeros(self.N)
            
        # now evaluate the expectation of the posterior
        P = C.copy().astype(float)
        
        if scipy.sparse.issparse(P):
            P = P.todense()
        P += alpha
        w = P.sum(axis=1)
        P /= w # normalize

        return P, w
    
    
    def _eigmode_partials(self, Abar, mode_index):
        """
        Computes the first-order Taylor approximation to the eigenvalue and
        eigenvector sensitivities.
        
        Parameters
        ----------
        Abar : matrix
            The average matrix
        
        Returns
        -------
        s : array_like, 1D
            the partials for the eigenvalue
        S : array_like, 2D
            the partials for the eigenvector

        Notes
        -----
        Straight from appendices in [1]. Some functionality loosely copied from
        the deprecated SSA library.
        """

        ek = np.zeros(self.N)
        ek[-1] = 1.0
        #z = np.zeros(self.N)  # turns out we may not need this

        if self.sparse:
            assert scipy.sparse.isspmatrix(Abar)

            P, L, U, Q = self._sparse_LU_decomp(Abar.todense())
            
            # solve for the projection vectors x and x_a from [1] Appendix A
            try:
                x = spla.gmres(U, ek)[0]
            except LinAlgError as e:          # if U is singular
                x = np.ones(self.N)           # introduce a stand-in

            xa = spla.gmres(L.T * P.T, ek)[0]
            
            # get the eigenvalue partials
            s = self._sparse_outer_product(xa, x, norm=True)

        else: # AKA, if dense
            if scipy.sparse.isspmatrix(Abar):
                Abar = Abar.todense()

            # perform a well-pivoted LU decomposition
            P, L, U, Q = self._dense_LU_decomp(Abar)

            # solve for the projection vectors x and x_a from [1] Appendix A
            try:
                x = scipy.linalg.solve_triangular(U, ek)
            except LinAlgError as e:          # if U is singular
                x = np.ones(self.N)           # introduce a stand-in
                
            #xa = self._nontrivial_backwards_substitution(L.T)
            LT = np.matrix(L, np.float64, copy=False).T
            PT = np.matrix(P, np.float64, copy=False).T
            xa = scipy.linalg.solve( LT * PT, ek ) # TJL: This is *not* what is in [1]
                                                   # I copied it from the ssaCalculator
                                                   # code, not sure if it's right. The
                                                   # equation in [1]  *does* present issues

            # get the eigenvalue partials
            s = np.outer(xa, x) / np.dot(xa, x)

        # next, the eigenvector partials
        if self.use_vectors:
            S = self._calculate_eigenvector_partials(P, L, U, s, mode_index)
        else: # no eigenvector calculation desired, skip it
            S = None
                    
        return s, S


    def _calculate_eigenvector_partials(self, P, L, U, s, mode_index):
        """
        Calculated the eigenvalue partial derivatives, dV_{lambda} / dp_{ij}

        Parameters
        ----------
        P, L, U : matrix
            the LU factors (and permutation matrix) of Abar

        Returns
        -------
        S : 3D array
            the eigenvalue partials, dV_{lambda} / dp_{ij}, such that
            
                      S[i,j,k] = dV_{lambda}[j] / dp_[i,k]

            thus, the shape of S is N x N x N, and the ordering of S
            is consistent with [1]
        """

        S = np.zeros( (self.N,)*3 ) # the array to store results in

        if self.sparse:
            # TJL: I'm not even sure a sparse implementation is possible
            # this algorithm will eat up a lot of memory no matter what...
            raise NotImplementedError("No sparse support for eigenvector"
                                      "calculations currently.")
        
        else: # dense

            # eq. B8 (bottom) in [1]
            b = self.vectors[:,mode_index].copy()
            d = self._solve_augmented_LU(P, L, U, b, self.vectors[:,mode_index], 0)

            # eq. B8 (top) in [1]
            z = np.zeros(self.N)
            for i in range(self.N):
                b = z.copy()
                b[i] = 1.0
                c_i = self._solve_augmented_LU(P, L, U, b, self.vectors[:,mode_index], 0)

                # eq. B9 in [1]
                #S[i,:,:] = np.outer(1.0/self.vectors[:,mode_index], c_i) + np.dot(s, d)
                #print S[i,:,:]
                for j in range(self.N):
                    S[i,:,j] = c_i / self.vectors[j,mode_index] + np.dot(s, d)

                # here is the "safe" code to use for the check
                if DEBUG:
                    for j in range(self.N):
                        test = c_i / self.vectors[j,mode_index] + np.dot(s, d)
                        assert np.all( S[i,:,j] == test )

        return S

        
    def _dense_LU_decomp(self, A):
        """
        Performs an LU decomposition with a pivot that places the
        smallest element of U in the lower right corner of the
        upper triangle. This smallest element (usually close to zero?)
        is then replaced with unity, which removes the possibility of
        a trivial solution later.
        """
        
        # first decomposition
        P, L, U = scipy.linalg.lu(A)

        # find the element along the matrix diagonal that has the smallest abs. val.
	diag = np.abs( np.diagonal(U) )
	smallest_index = diag.argmin()

        # generate a permutation matrix
        Q = np.identity(self.N)

        # exchange smallest_index row with dim-1 row
        # multiplying A by this matrix on both sides (Q.A.Q) will ensure that the
        # smallest element of U will be in the lower right corner.
        if smallest_index+1 != self.N :
            newrowx = np.copy( Q[self.N-1,:] )
            newrowy = np.copy( Q[smallest_index,:] )
            Q[smallest_index,:] = newrowx
            Q[self.N-1,:] = newrowy
			
            # recompute the decomposition
            Q = np.matrix(Q)
            A = np.matrix(A)
            P, L, U = scipy.linalg.lu( Q*A*Q )

        U[self.N-1, self.N-1] = 1.0 # the last element of U
	
        return P, L, U, Q


    def _sparse_LU_decomp(self, A):
        """
        Performs an LU decomposition with a pivot that places the
        smallest element of U in the lower right corner of the
        upper triangle.
        """

        # first decomposition
        UMF = umfpack.UmfpackContext()
        L, U, P, q, r, do_recip = UMF.lu(A)

        # find the element along the matrix diagonal that has the smallest abs. val.
        print type(U)
	diag = np.abs( np.diagonal(U) )
	smallest_index = diag.argmin()

        # generate a permutation matrix
        Q = np.identity(self.N)

        # exchange smallest_index row with dim-1 row
        # multiplying A by this matrix on both sides (Q.A.Q) will ensure that the
        # smallest element of U will be in the lower right corner.
        if smallest_index+1 != self.N :
            newrowx = np.copy( Q[self.N-1,:] )
            newrowy = np.copy( Q[smallest_index,:] )
            Q[smallest_index,:] = newrowx
            Q[self.N-1,:] = newrowy
			
            # recompute the decomposition
            L, U, P, q, r, do_recip = umfpack.lu( Q*A*Q )

        U = U.tolil()
        U[self.N-1, self.N-1] = 1.0
        U = U.tocsr()
                
        return P, L, U, Q
    
            
    def _sparse_outer_product(self, a, b, norm=False):
        """
        Compute the normalized outer product in a smart way
        so as to save memory. Should be fairly speedy, too.
        """

        shp = (len(a), len(b))
        
        a_nonzero_inds = np.where( a != 0.0 )[0]
        b_nonzero_inds = np.where( b != 0.0 )[0]

        dim1 = len(a_nonzero_inds)
        dim2 = len(b_nonzero_inds)
        X_nonzero_inds = np.zeros(( 2, dim1*dim2 ), dtype=np.int)
        X_nonzero_inds[0,:] = np.repeat(a_nonzero_inds, dim2)
        X_nonzero_inds[1,:] = np.tile(b_nonzero_inds, dim1)

        X = scipy.sparse.lil_matrix(shp)
        X[ X_nonzero_inds ] = 1.0

        A = scipy.sparse.dia_matrix((a, 0), shape=shp)
        B = scipy.sparse.dia_matrix((b, 0), shape=shp)

        R = ((X * A).T * B).T
        R = R.tocsr()

        if DEBUG:
            err = np.abs( np.sum( ( R.toarray() - np.outer(a, b) ) ) )
            if err > 10.**-6:
                print >> sys.stderr, "WARNING in SinghalSampler:"
                print >> sys.stderr, "\tLarge error in outer prod:", err
        if norm:
            R /= np.dot(a, b)

        return R


    def _solve_augmented_LU(self, P, L, U, b, augmentation, b_aug):
        """
        Solves the augemented (non-square) equation

            LUx = Pb

        We solve, in order:
            (1) Ly = Pb (via. forward substitution)
            (2) Ux = y  (back substitution)

        This method is specifically designed to solve the equations
        described in Appendix B of [1], in the manner outlined there.

        Parameters
        ----------
        P, L, U : matrix
            the NON-augmented LU factors
        b : array
            the RHS of the equation to solve
        augmentation : array
            the vector by which to augment A
        b_aug : float
            the element augmentation to b
        
        Returns
        -------
        x : the solution vector

        Notes
        -----
        Should work well for both sparse and dense matrices.
        TJL: Currently only dense is implemented, though.
        """

        # extract some size info
        n = U.shape[0]
        assert (n, n) == L.shape
        assert (n, n) == P.shape
        assert (n, n) == U.shape
        assert n == len(b)
        assert n == len(augmentation)

        # initialize data structures
        y = np.zeros(n)
        x = np.zeros(n)

        # let PbA = P * b augmented with b_aug
        Pb = np.dot(P, b)
        PbA = np.append( Pb, b_aug )
        del Pb

        # now, augment L and U
        LA = np.vstack(( L, augmentation ))
        UA = np.vstack(( augmentation, U ))
        
        # forward substitution (Ly = Pb)
        for i in range(n):
            y[i] = PbA[i] - np.sum(LA[i,:] * y) / LA[i,i]

        # back substitution (Ux = y)
        for i in range(n)[::-1]:
            x[i] = y[i] - np.sum(UA[i,:] * x) / UA[i,i]
            
        #print x; sys.exit(1)
        return x


    def _nontrivial_backwards_substitution(self, U):
        """
        Solves Ux = 0, where U is an upper trianglular matrix,
        subject to the constraint that x[-1] = 1.0. Does this by
        backward substitution.
        
        Returns: x, an array
        """

        if scipy.sparse.issparse(U):
            U = U.tolil()
            sparse = True
        else:
            sparse = False

        n = U.shape[0]
        assert n == U.shape[1]

        # initialize (the constraint)
        x = np.zeros(n)
        x[-1] = 1.0
        test = x.copy()

        for i in range(n-1)[::-1]:
            x[i] = -1.0 * np.sum(U[i,:] * x) / U[i,i]

        if DEBUG:
            assert np.all( x != np.zeros(n) )
            if not sparse:
                np.testing.assert_allclose( np.dot(U, x), test, atol=10.**-5 )
            else:
                np.testing.assert_allclose( U * x, test, atol=10.**-5 )
        
        return x

    
    
def test():
    """
    Perform a small series of tests on randomly generated data, and ensure
    that all disparate methods give the same answer.

    Benchmarked against SSA code from MSMBuilder, r1081, and successfully
    reproduced that reference on 8/10/12.                          -- TJL
    """

    print "WARNING: test() function is not meant for production runs."

    C = scipy.io.mmread('test_C.mtx')
    C += C.T

    indices = np.array([0, 1, 2, 3, 5])
    weights = np.array([1.0, 0.5, 0.1, 0.0, 0.3])

    sampler1 = SinghalSampler(indices, weights, prior='sparse')
    sampler1s = SinghalSampler(indices, weights, prior='sparse')
    sampler2 = SinghalSampler(indices, weights, prior='Jefferys', eigenvector_weight=0.5)
    sampler2s = SinghalSampler(indices, weights, prior='sparse', eigenvector_weight=0.5)
    
    print "\n\nANS:"
    print sampler1.sample(C, force_dense=True)
    #sp_err = np.sum( np.abs( sampler1s.sample(C, force_dense=True) - \
    #                         sampler1s.sample(C) )) / float(C.shape[0])
    #print "Error in sparse:", sp_err
    print sampler2.sample(C)
    #print sampler2s.sample(C)

    return


if __name__ == '__main__':
    test()
