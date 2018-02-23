import numpy as np

class Searcher :
    def __init__(self,index,site,td) :
        self.index = index
        self.site = site
        self.td = td

def init_searchers(omega_t,site_taus,deterministic=False) :
    """
    Initializes the searchers in the system. Note that this function will only
    be invoked from within the main simulation loop, so the "searchers" list that
    will be created is only visible internally, as a convenient data structure to
    keep track of what's going on, and not to the external world to see.
    """
    # get parameters from the omega_t matrix
    n = omega_t.size
    m = omega_t.sum()
    # init the searchers
    searchers = []
    # "searcher_sites" is a vector of indices containing the indices
    # of the sites occupied by the searchers
    searcher_sites = np.where(omega_t)[0]
    for s in xrange(m) :
        site = searcher_sites[s]
        if deterministic :
            td = site_taus[site]
        else :
            td = np.random.exponential(scale=site_taus[site])
        searcher = Searcher(s,site,td)
        searchers.append(searcher)
    return searchers

def run_chair_simulation(nsteps,omega_t_initial,T,site_taus,
                         seed=None,teq=0,tsample=1,deterministic=False) :
    # init the random number generator if it was passed
    if seed is not None :
        np.random.seed(seed)
    # make an internal copy of the initial omega_t matrix
    omega_t = omega_t_initial.copy()
    # no need to pass N,n,m through the arguments of the function
    n = omega_t.size
    m = omega_t.sum()
    # init jumping matrix
    J = np.zeros((n,n)).astype(np.int32)
    # init searchers
    searchers = init_searchers(omega_t,site_taus,deterministic=deterministic)
    # init sampling matrix
    nsamples = (nsteps-teq)/tsample
    samples = np.zeros((nsamples,n),dtype=bool)
    i_sample = 0
    # cycle on time
    for step in xrange(1,nsteps+1) :
        # cycle on the searchers
        for s in xrange(m) :
            searcher = searchers[s]
            # if the searcher has to stay longer on the site, skip it
            if step>searcher.td :
                # if not, get the elements corresponding to the transition matrix
                Tstar = T[searcher.site,:] * (~omega_t)
                Tstar /= Tstar.sum()
                # now get the next site
                next_site = np.random.choice(n,p=Tstar)
                # update omega matrix
                omega_t[searcher.site] = False
                omega_t[next_site] = True
                # update jumping matrix
                J[searcher.site,next_site] += 1
                # update searcher
                searcher.site = next_site
                if deterministic :
                    next_td = site_taus[searcher.site]
                else :
                    next_td = np.random.exponential(scale=site_taus[searcher.site])
                searcher.td = step + next_td
        # update samples
        if step>teq and (step-teq)%tsample==0 :
            samples[i_sample,:] = omega_t
            i_sample += 1
    return omega_t, J, samples

def init_omega_t(n,mu,sigma=None,seed=None) :
    # init the random number generator if it was passed
    if seed is not None :
        np.random.seed(init_seed)
    # init searcher numbers
    if sigma is not None :
        m = np.random.normal(loc=mu,scale=sigma).astype(np.int32)
        # ensures that the system has more searchers than available sites and that
        # it has at least one searcher
        if m>n : m=n
        if m<1 : m=1
    else :
        m = mu
    # init omega vector
    omega_t = np.zeros(n,dtype=bool)
    # fill with initial occupancy
    omega_t[np.random.choice(n,m,replace=False)] = True
    return omega_t

class JumpingModel :
    def __init__ (self,T,site_taus) :
        self.T = T
        self.site_taus = site_taus
        self.omega_t = {}
        self.J = {}
        self.samples = {}
        self.theta = {}
    def run(self,nsteps,mu,sigma,omega_t_initial,
           seed=None,teq=0,tsample=1,deterministic=False) :
        self.omega_t[mu],self.J[mu],self.samples[mu] = \
                run_chair_simulation(nsteps,omega_t_initial,self.T,self.site_taus,
                                    seed=seed,teq=teq,tsample=tsample,deterministic=deterministic)
        self.theta[mu] = self.samples[mu].sum(axis=0)/float(self.samples[mu].sum())

def H_to_L(model,Hsites,Lsites) :
    mus = model.samples.keys()
    mus.sort()
    nmus = len(mus)
    model.avH = {}
    model.avL = {}
    model.H_to_L = np.zeros(nmus)
    for i,mu in enumerate(mus) :
        occupancy = model.samples[mu].sum(axis=0)
        model.avH[mu] = occupancy[Hsites].mean()
        model.avL[mu] = occupancy[Lsites].mean()
        model.H_to_L[i] = model.avH[mu]/model.avL[mu]
