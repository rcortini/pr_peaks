import numpy as np

class Searcher :
    def __init__(self,index,site,td) :
        self.index = index
        self.site = site
        self.td = td

def init_searchers(omega_t,site_taus) :
    """
    Initializes the searchers in the system. Note that this function will only
    be invoked from within the main simulation loop, so the "searchers" list that
    will be created is only visible internally, as a convenient data structure to
    keep track of what's going on, and not to the external world to see.
    """
    # get parameters from the omega_t matrix
    N,n = omega_t.shape
    m = omega_t.sum(axis=1)
    # init the searchers
    searchers = []
    for i in xrange(N) :
        searchers_system = []
        # "searcher_sites" is a vector of indices containing the indices
        # of the sites occupied by the searchers
        searcher_sites = np.where(omega_t[i])[0]
        for s in xrange(m[i]) :
            site = searcher_sites[s]
            td = np.random.exponential(scale=site_taus[site])
            searcher = Searcher(s,site,td)
            searchers_system.append(searcher)
        searchers.append(searchers_system)
    return searchers

def run_chair_simulation(nsteps,omega_t_initial,T,site_taus,seed=None) :
    # init the random number generator if it was passed
    if seed is not None :
        np.random.seed(seed)
    # make an internal copy of the initial omega_t matrix
    omega_t = omega_t_initial.copy()
    # no need to pass N,n,m through the arguments of the function
    N,n = omega_t.shape
    m = omega_t.sum(axis=1)
    # init searchers
    searchers = init_searchers(omega_t,site_taus)
    # cycle on time
    for step in xrange(1,nsteps+1) :
        # cycle on the systems
        for i in xrange(N) :
            # cycle on the searchers
            for s in xrange(m[i]) :
                searcher = searchers[i][s]
                # if the searcher has to stay longer on the site, skip it
                if step>searcher.td :
                    # if not, get the elements corresponding to the transition matrix
                    Tstar = T[searcher.site,:] * (~omega_t[i,:])
                    Tstar /= Tstar.sum()
                    # now get the next site
                    next_site = np.random.choice(n,p=Tstar)
                    # update omega matrix
                    omega_t[i,searcher.site] = False
                    omega_t[i,next_site] = True
                    # update searcher
                    searcher.site = next_site
                    searcher.td = step + np.random.exponential(scale=site_taus[searcher.site])
    return omega_t

def init_omega_t(N,n,mu,sigma=None,seed=None) :
    # init the random number generator if it was passed
    if seed is not None :
        np.random.seed(init_seed)
    # init searcher numbers
    if sigma is not None :
        m = np.random.normal(loc=mu,scale=sigma,size=N).astype(np.int32)
        # ensures that no system has more searchers than available sites and that
        # every system has at least one searcher
        m[m>n] = n
        m[m<=1] = 1
    else :
        m = mu * np.ones(N).astype(np.int32)
    # init omega matrix
    omega_t = np.zeros((N,n),dtype=bool)
    # fill with initial occupancy
    for i in xrange(N) :
        omega_t[i,np.random.choice(n,m[i],replace=False)] = True
    return omega_t

class JumpingModel :
    def __init__ (self,T,site_taus) :
        self.T = T
        self.site_taus = site_taus
        self.omega_t = {}
        self.occupancy = {}
    def run(self,nsteps,mu,sigma,omega_t_initial) :
        self.omega_t[mu] = run_chair_simulation(nsteps,
                                                omega_t_initial,
                                                self.T,
                                                self.site_taus)
        self.occupancy[mu] = self.omega_t[mu].sum(axis=0)

def H_to_L(model,Hsites,Lsites) :
    mus = model.occupancy.keys().sort()
    nmus = len(mus)
    model.avH = {}
    model.avL = {}
    model.H_to_L = np.zeros(nmus)
    for i,mu in enumerate(mus) :
        model.avH[mu] = model.occupancy[mu][Hsites].mean()
        model.avL[mu] = model.occupancy[mu][Lsites].mean()
        model.H_to_L[i] = model.avH[mu]/model.avL[mu]
