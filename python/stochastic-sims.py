
# coding: utf-8

# In[46]:

#!/usr/bin/env python
get_ipython().magic('matplotlib inline')

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi
import pandas as pd
import scipy.stats as st 
import resource
resource.setrlimit(resource.RLIMIT_NOFILE, (1000,-1)) #allow many plots

RV217=pd.read_csv('data/RV217Mono.csv')

V_det=1 #copy/mL, V_det


# In[47]:

#all participants
params=pd.read_csv('data/estimatedIndividualParameters_corr2.csv')

#parameters we want to look at
pzl=['initT_mode','dS_mode','lBt0_mode','dI_mode','n_mode']
pzn=[r'$t_{det}$',r'$\delta_S$',r'$\log_{10}\beta$',r'$\delta_I$',r'$n$']

#set fixed parameters
tau=params['tau_mode'][0]
lp=params['lp_mode'][0]; p=10**lp
V0=params['V0_mode'][0] #deterministic viral load
aS=params['aS_mode'][0] #susecptible cell birth

g = 23.; #virus clearance rate [1/day]


# In[48]:

#look at results from best fits
table1=pd.DataFrame()
table1['id']=params['id']
for pz in pzl:
    table1[pz]=params[pz]


# In[49]:

num_ppts=len(params[param_names[0]])

num_ppts


# In[50]:

#random colors
get_colors = lambda n: list(map(lambda i: "#" + "%06x" % np.random.randint(0, 0xFFFFFF),range(n)))
cz_list=get_colors(num_ppts) # sample return:  ['#8af5da', '#fbc08c', '#b741d0', '#e599f1', '#bbcb59', '#a2a6c0']


# In[109]:


#function that solves stochastically using tau-leap method
def sto_sim(t0,tF,vol,ppt,verbose):

    #the functions to model the same process stochastically with the hybrid stochastic model
    #updating the rates of events
    num_rates=6; num_states=3; rl=np.zeros(num_rates); T=np.zeros([num_states,num_rates])
    def update_rates(X,ti):
        S,P,V=X; #state vars actually spelled out
        pr=np.random.poisson(p)        
        rl[0] = aS*vol;              T[:,0] =[1,0,0];  #constant production 
        rl[1] = dS*S;                T[:,1] =[-1,0,0]  #density dependent susceptible death
        rl[2] = (1-tau)*Bt0/vol*S*V; T[:,2] =[-1,0,-1] #unproductive cell infection
        rl[3] = tau*Bt0/vol*S*V;     T[:,3] =[-1,1,-1] #productive cell infection
        rl[4] = dI*P*(P/vol)**n;     T[:,4] =[0,-1,pr]  #productive cell death
        rl[5] = g*V;                 T[:,5] =[0,0,-1]  #viral clearance
        return rl,T

    idz,initT,dS,lBt0,dI,n=table1.iloc[ppt,:]
    
    Bt0 = 10**(lBt0); #infectivity rate [uL/cells-day], re-exponentiated

    ti=t0; sol=[]; tt=[] #initialize
    
    #calculate reproductive number (approx)
    R0=tau*Bt0*aS*p/(dS*g*dI)
    
    S_0 = aS*vol/dS #assumes T cells at equilibrium, and more in LN
    P_0 = 1. #start with a single infected cell
    V_0 = 1. #no virus yet

    x=np.array([S_0,P_0,V_0]); #initial state matrix
    
    dt=0.01; #time step
    cross_time=0
    extinct_time=0
    while ti < tF:
        x[x<1]=0 #make sure no negative numbers

        S,P,V=x; #state vars actually spelled out

        VL=V/vol #viral load
        
        sol.append(x) #the list of states
        tt.append(ti) #the list of times

        if P+V==0 and ti>0:
            extinct_time=ti
            if verbose:
                print(ppt,'went extinct! R0=',R0)
            break

        if VL>V_det/1e3: #more than a certain copies/mL
            cross_time=ti
            if verbose:
                print(ppt,'max VL reached, R0=',R0)
            break

        rl,T = update_rates(x,ti) #make new matrices
        events = np.random.poisson(rl*dt) #calculate events
        dx = np.sum(T*events,1) #calculate change in state
        x=x+dx #update state variable
        ti=ti+dt #update time

    if verbose:
        print(ppt,'made it past tF! R0=',R0)

    return np.array(tt),np.array(sol),R0,cross_time,extinct_time


# In[110]:

#simulate stochastic, 1 each just to visualize in comparison to ODE
#calculate all MLE R0

#volume of sim
vol=5e6; #use 5L
t0=0; tF=100;

R0l=[]
plt.figure(figsize=(8,3))
cts=[]; ets=[]
#for ppt in range(num_ppts):

#choose N random sims
NN=5
for ppt in np.random.randint(0,num_ppts,NN):
    tt,sol,R0,ct,et=sto_sim(t0,tF,vol,ppt,verbose=True)       
    R0l.append(R0)
    cts.append(ct)
    ets.append(et)
    
    plt.subplot(131)
    plt.step(tt,(sol[:,0]-sol[0,0])/sol[0,0]*100,color=cz_list[ppt])
    plt.ylabel('% deviation from S0')
    plt.xlabel('time days after t_0')
    
    plt.subplot(132)
    plt.step(tt,sol[:,1],color=cz_list[ppt])
    plt.ylabel('infected cells')
    plt.xlabel('time days after t_0')

    VL=sol[:,2]/vol*1000 #viral load as usual units copies/mL in plasma
    plt.subplot(133)
    plt.step(tt[VL>0],np.log10(VL[VL>0]),color=cz_list[ppt])
    plt.ylim([-1,0.2])
    plt.ylabel('viral load log10(copies/mL)')
    plt.xlabel('time days after t_0')
    plt.axhline(np.log10(thresholdVL),ls='--',color='k')

plt.tight_layout()
plt.savefig('figures/sim_sto.pdf',dpi=600)

#also plot R0
plt.figure(figsize=(2,3))
plt.scatter(np.linspace(0.7,1.3,len(R0l)),R0l,alpha=0.5,c=cz_list)
plt.boxplot(R0l,widths=0.7)
plt.xticks([1],[r'$\mathcal{R}_0$'],fontsize=12)
plt.axhline(1,color='k',ls='--')
#plt.ylim([0,16])
#plt.yticks(range(16))
plt.tight_layout()
plt.savefig('figures/R0.pdf',dpi=600)


# In[111]:

#simulate stochastic with replicates for each individual

num_replicates=5
max_reps=20
t0=0; tF=30;

R0l=[]; cts=[]; ets=[]; ets_ind=[]
print('ppt','total_reps')
for ppt in range(num_ppts):
    good_reps=0
    tot_reps=0
    ets_ind.append([])
    while (good_reps < num_replicates) and (tot_reps < max_reps):
        tt,sol,R0,ct,et=sto_sim(t0,tF,vol,ppt,verbose=False)       
        R0l.append(R0); cts.append(ct); ets.append(et)
        if et>0:
            ets_ind[ppt].append(et)
            good_reps+=1
        tot_reps+=1
            
    print(ppt,tot_reps)



# In[112]:

#plot extinction times
plt.figure(figsize=(7,3))

plt.subplot(121)
ets=np.array(ets)
plt.hist(ets[ets>0],np.arange(0,max(ets)+1),color='navy',ec='white')
plt.ylabel('counts')
plt.xlabel('extinction time (days)')

#plot breaktrhough times
plt.subplot(122)
cts=np.array(cts)
plt.hist(cts[cts>0],np.arange(0,max(cts)+1),color='navy',ec='white')
plt.xlabel('detectable time (days)')

plt.tight_layout()
plt.savefig('figures/timesofinterest.pdf',dpi=600)



# In[114]:

#correlate with R0

plt.figure(figsize=(7,3))

R0l=np.array(R0l)

# extinct times
plt.subplot(121)
plt.scatter(R0l[ets>0],ets[ets>0],c=cz_list)
plt.xlabel(r'$\mathcal{R}_0$',fontsize=14)
plt.ylabel('extinction time (days)')
#plt.loglog()
#plt.ylim([1e-2,10])

# breaktrhough times
plt.subplot(122)
plt.scatter(R0l[cts>0],cts[cts>0],c=cz_list)
plt.xlabel(r'$\mathcal{R}_0$',fontsize=14)
plt.ylabel('breakthrough time (days)')
#plt.loglog()

plt.tight_layout()
plt.savefig('figures/R0correlations.pdf',dpi=600)



# In[125]:

#make a plot of each individual's p(t_det), p(t_sto), p(t_net)
BINS=np.arange(0,30,2)

NN=5
for ppt in np.random.randint(0,num_ppts,NN):
#for ppt in range(num_ppts):

    plt.figure(figsize=(7,1.5))
    tdet_mu=params['initT_mean'].iloc[ppt] #mean t_det
    tdet_sd=params['initT_sd'].iloc[ppt] #mean t_det

    tdet_dist=np.random.normal(tdet_mu,tdet_sd,[100])
    
    counts,bins=np.histogram(tdet_dist,bins=BINS)
    ptdet=counts/sum(counts)
    
    plt.subplot(131)
    #plt.hist(tdet_dist)
    #plt.step(bins[1:],ptdet,label='det')
    plt.plot(bins[1:],ptdet,label='deterministic',color='blue')
    plt.fill_between(bins[1:],np.zeros(len(ptdet)),ptdet,alpha=0.5,color='blue')
    plt.xlabel('t_det')
    plt.ylabel('density')
    #plt.title('deterministic',fontsize=8)
    #plt.xticks(np.arange(min(BINS),1,5))

    if len(ets_ind[ppt])>1:
        
        tsto_dist=np.array(ets_ind[ppt])
        counts,bins=np.histogram(tsto_dist,bins=BINS)
        ptsto=counts/sum(counts)
        
        plt.subplot(132)
        plt.plot(bins[1:],ptsto,label='stochastic',color='orange')
        plt.fill_between(bins[1:],np.zeros(len(ptsto)),ptsto,alpha=0.5,color='orange')
        #plt.step(bins[1:],ptsto,label='sto')
        #plt.title('stochastic',fontsize=8)
        plt.xlabel('t_sto')
        #plt.xticks(np.arange(min(BINS),1,5))

        draw_det=np.random.choice(bins[:-1], 10000, p=ptdet)
        draw_sto=np.random.choice(bins[:-1], 10000, p=ptsto)
        counts,bins=np.histogram(draw_det+draw_sto,bins=BINS)
        ptnet=counts/sum(counts)
        #convo=np.convolve(ptdet, ptsto, mode='full')
    
        plt.subplot(133)
        plt.plot(bins[1:],ptnet,label='net',color='green')
        plt.fill_between(bins[1:],np.zeros(len(ptnet)),ptnet,alpha=0.5,color='green')
        #plt.title('net',fontsize=8)
        #plt.annotate()
        #plt.title('mode = '+str(st.mode(ptnet)))
        plt.xlabel('t_0')
        #plt.xticks(np.arange(min(BINS),1,5))

    #plt.legend()
    plt.tight_layout()
    plt.savefig('figures/individual_ests/'+str(params['id'][ppt])+'.pdf',dpi=600)


# In[ ]:



