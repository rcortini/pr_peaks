import numpy as np

def init_omega_t(n,mu,seed=None) :
    # init the random number generator if it was passed
    if seed is not None :
        np.random.seed(seed)
    # init omega vector
    omega_t = np.zeros(n,dtype=bool)
    # fill with initial occupancy
    omega_t[np.random.choice(n,mu,replace=False)] = True
    return omega_t

def run_chair_simulation(nsteps,omega_t_initial,H,site_taus,boost,
                         seed=None,teq=0,tsample=1) :
    # init the random number generator if it was passed
    if seed is not None :
        np.random.seed(seed)
    # make an internal copy of the initial omega_t vector
    omega_t = omega_t_initial.copy()
    # no need to pass N,n,m through the arguments of the function
    n = omega_t.size
    m = omega_t.sum()
    # init jumping matrix
    J = np.zeros((n,n)).astype(np.int32)
    # init searcher_sites
    searcher_sites = -1 * np.ones(n,dtype=np.int32)
    searcher_sites[omega_t] = np.arange(m)
    # init searcher_times
    searcher_times = np.zeros(m)
    for s in xrange(m) :
        searcher_id = searcher_sites[s]
        if searcher_id != -1 :
            searcher_times[searcher_id] = site_taus[s]
    # init sampling matrix
    nsamples = (nsteps-teq)/tsample
    samples = np.zeros((nsamples,n),dtype=bool)
    i_sample = 0
    # cycle on time
    for step in xrange(1,nsteps+1) :
        # cycle on the searchers
        for s in xrange(n) :
            searcher_id = searcher_sites[s]
            if searcher_id != -1 :
                # if the searcher has to stay longer on the site, skip it
                if step>=searcher_times[searcher_id] :
                    # now get the next site
                    next_site = np.random.choice(np.where(~omega_t)[0])
                    # update omega matrix
                    omega_t[s] = False
                    omega_t[next_site] = True
                    # update jumping matrix
                    J[s,next_site] += 1
                    # update searcher_sites
                    searcher_sites[s] = -1
                    searcher_sites[next_site] = searcher_id
                    # update searcher_times
                    next_td = site_taus[next_site]
                    for contact in H[next_site] :
                        contact_searcher_id = searcher_sites[contact]
                        if contact_searcher_id != -1 :
                            next_td *= boost
                            searcher_times[contact_searcher_id] = step + next_td
                            break
                    searcher_times[searcher_id] = step + next_td
        # update samples
        if step>teq and (step-teq)%tsample==0 :
            samples[i_sample,:] = omega_t
            i_sample += 1
    return omega_t, J, samples

class JumpingModel :
    def __init__ (self,H,site_taus,boost) :
        self.H = H
        self.site_taus = site_taus
        self.boost = boost
        self.J = {}
        self.samples = {}
        self.theta = {}
    def run(self,nsteps,omega_t_initial,seed=None,teq=0,tsample=1) :
        mu = omega_t_initial.sum()
        omega_t,self.J[mu],self.samples[mu] = \
            run_chair_simulation(nsteps,
                                 omega_t_initial,
                                 self.H,
                                 self.site_taus,
                                 self.boost,
                                 seed=None,teq=0,tsample=1)
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
