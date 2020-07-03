import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as alg
import pickle
import scipy, scipy.optimize
import scipy.special as special
import pandas as pd
import clipboard
#Get the i-th bit from 'x'
def get_bit(x, i):
    return (x >> i)&1
#Returns a new integer with the i-th bit = 1
def set_bit(value, i):
    return value | (1<<i)
#Returns a new integer with the i-th bit = 1
def clear_bit(value, i):
    return value & ~(1<<i)
#Counts the number of '1's in the binary representation of 'value' 
def count_bits(value, n):
    cnt = 0
    for i in range(n):
        cnt += (value >> i)%2
    return cnt
#Computes the vector local z-magnetizations on the state 'state'
def sigma_z(state, n_particles):
    mgn = np.zeros(n_particles)
    for i in range(n_particles):
        if(get_bit(state,i)==0):
            mgn[i] = + 1
        else:
            mgn[i] = - 1
    return mgn
#Computes the total z-magnetization on the state 'state'
def magnetization_z(state, n_particles):
    mgn = 0
    for i in range(n_particles):
        if(get_bit(state,i)==0):
            mgn = mgn + 1
        else:
            mgn = mgn - 1
    return mgn
#Returns all magnetizations from all the possible 'n_particles' many body states
def create_magnetizations(n_particles):
    sz = 2**n_particles
    mgns = np.zeros(sz)
    for state in range(sz):
        mgns[state] = magnetization_z(state, n_particles)
    return mgns 
#Local spin for each state and site
def get_spin_map(n_particles):
    full_sz = 1 << n_particles
    bit_map = np.zeros((full_sz, n_particles))
    for state in range(full_sz):
        bit_map[state,:] = np.array([(state>>i) & 1 for i in range(L)])
    return 1-2*bit_map
#Given x_i returns d_ij = x_i-x_j
def diff_matrix(x):
    temp = np.repeat(x.reshape(1,-1), len(x), axis=0)
    return (temp - x.reshape(-1,1)).T
#Simple curve intersection -> best estimate given by data
def curve_intersection(f,g):
    delta = np.abs(f-g)
    return np.argmin(delta)
def random_heseinberg_diag_sector(n_particles,h_strength, n_samples, times, verbosity=1, verbosity_s=1, pbc=True, sector_magn=0, compute_Ct=False):
    ####################################################
    ##############      HAMILTONIAN       ##############
    ####################################################

    #recovers the magnetizations for all possible many body states of 'n_particles'
    mgns = create_magnetizations(n_particles)
    #takes all the states (sector) with a given magnetization. Default = 0
    sector = np.argwhere(mgns==sector_magn).reshape(-1).astype('int')
    #amount of states in the chosen sector
    sz = len(sector)
    #dictionary to map sectpr to indices
    sect_dict =  {sector[i]: i for i in range(sz)} 
    #create the Hamiltonian using only sector states
    H0 = np.zeros((sz,sz))
    #decide whether use PBC or OBC
    if(pbc):
        max_site = n_particles
    else:
        max_site = n_particles-1
    spin_map =  get_spin_map(n_particles)
    #for each state in the sector it contains the local magnetizations
    spin_sector = spin_map[sector,:]
    #I need to address both the site and the index -> same values as in 'sect_dict'
    for idx,state in enumerate(sector):
        #build the zz interaction part of the Hamiltonian
        for site in range(max_site):
            next_site = (site+1)%n_particles
            #if site and the next have same magnetization we have a positive contribution (assuming J_zz==1)
            if(get_bit(state, site)==get_bit(state,next_site)):
                H0[idx, idx] +=  1/4
            #otherwise is negative
            else:
                H0[idx, idx] += -1/4

        #build the xy interaction part of the Hamiltonian
        for i in range(max_site):
            #next particle site
            j = (i+1)%n_particles
            #if the state bits i and i+1 are different the contribution is non-zero
            if(get_bit(state,i) != get_bit(state,j)):
                #i need to flip the two spins
                mask = (1<<i)+(1<<j)
                #the 'bra' state with the flipped spins
                #since flipping two spins doesn't change the total magnetization
                #i stay in the same spin sector ([H,S_z]=0)
                state_p = state^mask
                #recovers from dictionary its index
                idx_p = sect_dict[state_p]
                #positive contribution
                H0[idx, idx_p] += 2/4
    #number of different fields to try
    n_fields = len(h_strength)
    #indices in the spectrum to look -> following original paper we focus on the middle one third
    start_idx =int(sz/3)
    stop_idx = int(2*sz/3)

    ####################################################
    ##############     GLOBAL RESULTS     ##############
    ####################################################

    #spectrum for all samples and fields
    E = np.zeros((n_fields, n_samples, stop_idx-start_idx))
    #local magnetizations
    m = np.zeros((n_fields, n_samples, stop_idx-start_idx, L))
    
    ####################################################
    ##############    STATE RESULTS   ##############
    ####################################################
    #alternating sign state
    #-+-+-+...-+
    state0 = 0
    for i in range(n_particles):
        if(i%2==1):
            state0 += 1<<i

    psi0 = np.zeros(sz)
    psi0[sect_dict[state0]]=1
    #psi0 = np.ones(sz)/np.sqrt(sz)
    C_t = np.zeros((n_fields, n_samples, t_steps))
    ipr = np.zeros((n_fields, n_samples))
    M_wave = np.zeros((n_fields, n_samples))
    wave_cfs = np.exp(1.0j*2*np.pi*np.arange(n_particles)/n_particles)
    if(compute_Ct):
        c_path = None
        path = None
    for k, strg in enumerate(h_strength):
        if(verbosity > 0 and k % verbosity == 0):
            print(("%2.f" % (100*k/n_fields)) + "%")

        for s in range(n_samples):
            if(n_particles > 10 and verbosity_s > 0 and s % verbosity_s == 0):
                print(("Sample: %2.f" % (100*s/n_samples)) + "%")
            #the random field
            z_field = np.random.rand(n_particles)*(2*strg)-strg
            #start from the deterministic hamiltonian
            H = H0.copy()
            #loops for all state in the sector and their index
            for idx,state in enumerate(sector):
                #the function sigma_z recovers all local magnetizations
                H[idx, idx] += 0.5*np.sum(sigma_z(state, n_particles)*z_field) #1/2 sum_i sigma_i^z h_i
            spectrum, vecs = alg.eigh(H)
            c_vecs = np.conj(vecs) #actually, vecs are real since H is real and symmetric
            #<m|s_j^z|n>
            sjz_mn = np.einsum('sm,sn,sj', c_vecs, vecs, spin_sector, optimize=['einsum_path',(0,2),(0,1)]) #indices are nmj
            #psi0 in the eigenbasis
            coefs = np.matmul(c_vecs.T,psi0)
            #ipr of the state
            ipr[k,s] = 1/np.sum(np.abs(coefs)**4)
            if(compute_Ct):
                if(n_particles > 10 and verbosity_s > 0 and s % verbosity_s == 0):
                    print("Time evolution...")
                for t_idx,t in enumerate(times):   
                    #time evolution
                    u_nt = np.exp(-1.0j*spectrum*t)
                    if(c_path is None):
                        c_path, _ = np.einsum_path('jlm, jln, l', sjz_mn, sjz_mn, np.conj(u_nt))

                    #<m|\sigma^z (0)\sigma^z(t)|n>    
                    two_point = np.einsum('jlm, jln, l', sjz_mn, sjz_mn, np.conj(u_nt), optimize=c_path)/n_particles #indices are nm
                    if(path is None):
                        path, _ = np.einsum_path('n,m,nm,n', coefs, np.conj(coefs), two_point, u_nt)
                    C_t[k,s,t_idx] = np.abs(np.einsum('n,m,nm,n', coefs, np.conj(coefs), two_point, u_nt, optimize=path))
            #stores spectrum
            E[k,s,:] = spectrum[start_idx:stop_idx]

            #probabilities
            p = np.abs(vecs)**2

            #computes local magnetizations
            local_m = np.matmul(spin_sector.T, p).T
            m[k, s, :, :] = local_m[start_idx:stop_idx,:]

            #wave of spin part (see paper)
            f_njk = np.einsum('sn,sj,sk',p, spin_sector, spin_sector) #(L,L,sz)
            f_nj = np.einsum('sn,sj', p, spin_sector) #(L,sz)
            M_wave[k,s] = np.mean(1-np.abs(np.matmul(wave_cfs, f_nj))**2/np.einsum('jkn,j,k', f_njk, wave_cfs, np.conj(wave_cfs))).real

            
    return E,m,times,ipr, M_wave, C_t

LOAD_SUMMARY = False
SAVE_SUMMARY = False
compute_Ct = False #leave false for just computing static quantities -> much faster
t_steps = 50
times = np.logspace(-2,4,t_steps)

if(LOAD_SUMMARY):
    summary = pickle.load(open('summary.pkl', 'rb'))
    if(len(summary)==8):
        [h_strg, L_s, n_samples, dm_s, r_s, IPR_s, M_s, C_s] = summary
    else:
        [h_strg, L_s, n_samples, dm_s, r_s, IPR_s, M_s] = summary
    print('Samples present:', n_samples)
else:
    #actual work
    h_strg = [ 0.6, 1.0, 2.0, 2.7, 3.6, 5.0, 8.0 ]
    #h_strg = np.linspace(2.7, 3.7, 10)
    print('Field strengths to probe:')
    print(h_strg)
    L_s = [8,10,12,14]
    n_samples = [10000,10000,1000,50]
    dm_s = []
    r_s = []
    IPR_s = []
    M_s = []
    C_s = []
    for idx, L in enumerate(L_s):
        print('L:',L)
        E,m,times,ipr, M_wave, C_t = random_heseinberg_diag_sector(L,h_strg,n_samples[idx],times,pbc=True,verbosity=1, compute_Ct=compute_Ct)
        deltas = np.abs(np.diff(E, axis=2))
        #these are delta^n and delta^(n+1) where n is the eigenstate
        seq1, seq2 = deltas[:,:,0:-1],deltas[:,:,1:]
        #ratios
        r = np.minimum(seq1, seq2)/np.maximum(seq1,seq2)
        #average ratios (we average over the realizations of disorder and on the eigenstates)
        avg_r = np.mean(r, axis=(1,2))
        avg_dm = np.mean(np.abs(np.diff(m, axis=2)), axis=(1,2,3))
        r_s.append(avg_r)
        dm_s.append(avg_dm)
        IPR_s.append(ipr)
        M_s.append(M_wave)
        C_s.append(C_t)
        print('--------------------')
    if(SAVE_SUMMARY):
        summary = [h_strg, L_s, n_samples, dm_s, r_s, IPR_s, M_s, C_s]
        pickle.dump(summary, open('summary.pkl', 'wb'))




if(compute_Ct):
    f_t = lambda t,a,tau,w1,phi,b,z,c,eta,w2, : a*np.exp(-t/tau)*np.cos(w1*t+phi)+b*(t**(-z))*(1+c*(t**(-eta))*np.sin(w2*t+phi))
    fig, ax = plt.subplots(figsize=(16,10), nrows=len(L_s))
    try:
        len(ax)
    except:
        ax = [ax]
    curve_fits =  []
    error_fits = []
    for idx, L in enumerate(L_s):
        avg_Ct = np.mean(C_s[idx], axis=1)
        for idx2, h in enumerate(h_strg):
            ax[idx].plot(times, avg_Ct[idx2,:], '-o', label='h: %2.1f' % h)
            try:
                p_opt, p_cov = scipy.optimize.curve_fit(f_t, times, avg_Ct[idx2,:], bounds=(0,[1,np.inf,np.inf,2*np.pi,1,np.inf,1,np.inf,np.inf]), maxfev=10000)
                p_err = np.sqrt(np.diag(p_cov))
                curve_fits.append(np.concatenate(([L,h], np.array(p_opt))))
                error_fits.append(np.concatenate(([L,h], p_err)))
            except:
                pass
        ax[idx].set_xscale('log')
        ax[idx].set_yscale('log')
        ax[idx].grid()
        ax[idx].legend()
        ax[idx].set_xlabel('t')
        ax[idx].set_ylabel('$C_t, L=%i$'%L)
    plt.savefig('Images/corr.png')
    plt.show()
    curve_fits = pd.DataFrame(data=curve_fits, columns=['L','h','a','tau','w1','phi','b','z','c','eta','w2'])
    error_fits = pd.DataFrame(data=error_fits, columns=['L','h','a','tau','w1','phi','b','z','c','eta','w2'])
    #clipboard.copy(curve_fits.to_latex(index=False)+'\n\r' +error_fits.to_latex(index=False))
    #exit()

#IPR plot
fig, ax = plt.subplots()
for idx, L in enumerate(L_s):
    ax.plot(h_strg, np.mean(IPR_s[idx],axis=1)/scipy.special.binom(L,L/2), '-o', label='L: %i' %L)
plt.grid()
plt.xlabel('$h$')
plt.ylabel('$Inverse Participation Ratio$')
plt.legend()
#plt.title('IPR for a $T=\\infty$ state')7
plt.savefig('Images/ipr.png')
plt.show()

#DeltaM plot
avg_dm = np.array(dm_s)
fig, ax = plt.subplots()
for idx2, h in enumerate(h_strg):
    ax.plot(L_s, np.log(avg_dm[:,idx2]) , '-o', label='h: %2.1f' % h)
plt.grid()
plt.xlabel('$L$')
plt.ylabel('$\log < |m^{(n+1)}_{i \\alpha}-m^{(n)}_{i \\alpha}|>$')
plt.legend()
plt.savefig('Images/dm.png')
plt.show()

#Spin wave
fig, ax = plt.subplots()
for idx, L in enumerate(L_s):
    ax.plot(h_strg, np.mean(M_s[idx],axis=1), '-o', label='L: %i' %L)
plt.grid()
plt.xlabel('$h$')
plt.ylabel('$< f^{(n)}_\\alpha >$')
plt.legend()
plt.savefig('Images/m_factor.png')
plt.show()

#DeltaE ratios
fig, ax = plt.subplots()
for idx, L in enumerate(L_s):
    ax.plot(h_strg, r_s[idx], '-o', label='L: %i' %L)
plt.grid()
plt.xlabel('$h$')
plt.ylabel("$< r^{(n)}_\\alpha >$")
plt.legend()
plt.savefig('Images/energy_ratio.png')
plt.show()

