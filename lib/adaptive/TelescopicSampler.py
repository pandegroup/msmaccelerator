
# TJL UNFINISHED CODE -- DO NOT ATTEMPT TO USE

class dirichlet(object):
    """ A collection of useful functions for dealing with Dirichlet distributions.
        While these functions were generated with a mind for applications to MSMs,
        they should be generally applicable. """


    def __init__(self):
        raise NotImplementedError()


    def mv_beta( array ):
        return np.product( special.gamma( array ) ) / special.gamma( np.sum( array ) )

    def max_likelihood_params( counts ):
        """ Returns the most likely parameters given a set of counts
        """
        C_sums = np.asarray( counts.sum(axis=1) ).flatten()
        D = sparse.dia_matrix( (1.0 / C_sums ,0), C.shape).tocsr()
        T = D.dot(C)
        return T

    def max_posterior_params( counts ):
        """ Returns the parameters at the mode of the posterior

        Note: formula for the mode (x_i)

              x_i = ( counts_i + alpha_i - 1 ) / ( \sum_i ( counts_i + alpha_i ) - N )

        where N is the size of the counts matrix

        Arguments:
          counts - the observed counts
        """
       
        # check this for speed later 
        N = counts.shape[0]
        D = np.asarray( counts.sum(axis=1) ).flatten() + self.alpha - N
        T = (counts + self.alpha - 1.0) * (1.0 / D)
        return T

    
class TelescopicSampler(dirichlet):
    """
    Combines data from many forcefields to provide a best estimate
    of the transition matrix parameters for some gold standard
    forcefield.


    Argument Notes    

    uniform_prior:   sets the prior to be uniform. If false, then takes all
                     Dirichlet alpha values to be 1/N


    Important Latent Variables

    alpha :          the Dirichlet prior parameter, which is taken to be a
                     constant across all FFs 

    """


    def __init__(self, uniform_prior=False):
        raise NotImplementedError()
        self.num_ff    = 0         
        self.uniform_prior = uniform_prior
        self.ff_counts = {}        # dict to hold all counts data
        self.gold      = ''        # gold standard FF's name
        self.Z         = {}        # the coupling parameters
        self.Z_hist    = {}        # the history of the z parameters

    def add_ff( self, name, counts, is_gold=False ):

        if name in ff_counts.keys():
            print "WARNING: %s already in ff_counts for telescopic_estimator" % name
        else:
            self.num_ff += 1

        self.ff_counts[ name ] = counts

        if is_gold:
            self.gold = name


    def update_counts( self, name, counts ):
        assert name in ff_counts.keys()
        self.ff_counts[ name ] = counts

    def _initialize_alpha(self):
        if self.uniform_prior:
            self.alpha = 1.0
        else:
            self.alpha = 1.0 / self.ff_counts[ self.gold ].shape[0]

    def estimate_gold_counts( self ):

        assert self.num_ff > 0 
 
        self._initialize_alpha()
        C_hat = self.ff_counts[ self.gold ] # start with the obs. gold counts

        for ff in self.ff_counts.keys():
            self._estimate_z( name )
            C_hat += self.ff_counts[name] * self.Z[name]

        C = self.max_posterior_params( C_hat )
        return C

    def _add_Zhist( self, z, name ):
        if name in self.Z_hist.keys():
            self.Z_hist[name].append(z)
        else:
            self.Z_hist[name] = [z]

    def _estimate_z( self, name, verbose=False ):
        """ Provides the MLE estimate for the coupling parameters z_i
            for the specified forcefield 'name' """

        # first, maximize the likelihood of FF 'X' with the given prior
        T_X = self.max_posterior( self.ff_counts[name] )
        N_X = np.sum( self.ff_counts[name], axis=1)       # the num. observations 
        C_X = T_X * N_X                                   # estimated counts

        # next, optimize the log likelihood
        N     = counts.shape[0]
        z_opt = np.zeros( N )

        guess = 1.0 / C_X.sum()
        for i in range(N):

            C_X_i = C_X[i,:]
            C_Y_i = self.ff_counts[self.gold][i,:]

            def objective( z ):
                """ the log-likelihood objective function - minize me """
                LL = np.sum( special.gammaln( C_Y_i + z * C_X_i + 1 ) ) \
                             - special.gammaln( z * C_X_i.sum() )
                return -LL

            def D_objective( z ):
                """ The (semi)-analyical derivative of the objective function """
                dLL = np.sum( C_X_i * special.psi( C_Y_i + z * C_X_i + 1 )) \
                              - C_X_i.sum() * special.psi( z * C_X_i.sum() )
                return -dLL

            z_opt[i], info = optimize.minimize( objective, guess, full_output=True,
                                                method='BFGS', jac=D_objective )
            if verbose: print info
            guess = z_opt[i] # set the previous answer as the initial guess for the next state
            
        self.Z[name] = z
        self._add_Zhist( z, name )
