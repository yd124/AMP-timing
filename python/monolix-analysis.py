#DBR 9/26/19
#coding: utf-8
#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi
import pandas as pd
import scipy.stats as st 
import resource
resource.setrlimit(resource.RLIMIT_NOFILE, (1000,-1)) #allow many plots

#upload RV217 data
RV217=pd.read_csv('data/RV217Mono.csv')

#all participants, change this filename!
params=pd.read_csv('data/estimatedIndividualParameters_corr2.csv')
num_ppts=len(params['id'])

#random colors
get_colors = lambda n: list(map(lambda i: "#" + "%06x" % np.random.randint(0, 0xFFFFFF),range(n)))
cz_list=get_colors(num_ppts) 

#parameters we want to look at
pzl=['initT_mode','dS_mode','lBt0_mode','dI_mode','n_mode']
pzn=[r'$t_{det}$',r'$\delta_S$',r'$\log_{10}\beta$',r'$\delta_I$',r'$n$']

#set fixed parameters
tau=params['tau_mode'][0]
lp=params['lp_mode'][0]; p=10**lp
V0=params['V0_mode'][0] #deterministic viral load
aS=params['aS_mode'][0] #susecptible cell birth
g = 23.; #virus clearance rate [1/day]

#plot boxplots of best fits parameters
plt.figure(figsize=(7,2.5))
pind=0
table1=pd.DataFrame()
table1['id']=params['id']
for pz in pzl:
    plt.subplot(1,5,pind+1)
    plt.boxplot(params[pz],widths=0.7)
    plt.xticks([1],[pzn[pind]],fontsize=12)
    pind+=1
    table1[pz]=params[pz]
plt.tight_layout()
plt.savefig('figures/estimates.pdf',dpi=600)

#print out table
pd.DataFrame.to_csv(table1,'table1.csv')

#only plot the histogram of t_det
plt.figure(figsize=(4,3))    
tmin=70
plt.hist(params[pzl[0]],bins=np.arange(0,70,2),color='gray',ec='white')
plt.xlim([0,tmin])
plt.xlabel(r'$t_{det}$' + ' (days before first positive)')
plt.ylabel('counts')
plt.xticks(np.arange(0,tmin,5),np.arange(0,tmin,5))
plt.tight_layout()
plt.savefig('figures/tdet2.pdf',dpi=600)

print('95% interval of tdet = ',np.quantile(params[pzl[0]],0.0275),np.quantile(params[pzl[0]],0.975))

#ODE model
def ode_sim(X,t,aS,dS,Bt0,tau,dI,n,p,g):
    dY = np.zeros(3);
    S=X[0]; P=X[1]; V=X[2];
    dY[0] = aS - dS*S - Bt0*S*V     #susceptible cells [cells/uL]
    dY[1] = tau*Bt0*S*V - dI*P*P**n #productively infected cells [cells/uL]
    dY[2] = p*P - g*V - Bt0*S*V     #productive virus [virus/uL]
    return dY

#simulate deterministic all on same plots
t0=0 #true infection time
tF=70 #stop at 100 days
ts=np.linspace(t0,tF,1e5)
plt.figure(figsize=(10,3))
VL_peak=[]
VL_setpt=[]
for ppt in range(num_ppts):
    #get Fabian's MLE parameters
    idz,initT,dS,lBt0,dI,n=table1.iloc[ppt,:]

    Bt0 = 10**(lBt0); #infectivity rate [uL/cells-day], re-exponentiated

    S_0 = aS/dS #assumes T cells at equilibrium
    V_0 = V0/1e3 #1e-4 #no virus yet
    P_0 = V_0*tau*g/p #start with a single infected cell

    X0=[S_0,P_0,V_0]; #initial state matrix

    sol=spi.odeint(ode_sim, [S_0,P_0,V_0], ts, (aS,dS,Bt0,tau,dI,n,p,g), mxstep=10000) 
    VL=np.log10(sol[:,2]*1e3) #viral load as usual units copies/mL
    
    VL_peak.append(np.max(VL))
    VL_setpt.append(np.mean(VL[-100:-1]))

    plt.subplot(131)
    plt.step(ts,sol[:,0]-sol[0,0],color=cz_list[ppt])
    plt.ylabel('deviation from S0 (cells/mL)')
    plt.xlabel('time days after first positive')
    
    plt.subplot(132)
    plt.step(ts,sol[:,1],color=cz_list[ppt])
    plt.ylabel('infected cells (per mL)')
    plt.xlabel('time days after first positive')

    plt.subplot(133)
    plt.step(ts,VL,color=cz_list[ppt])
    plt.ylabel('viral load log10(copies/mL)')
    plt.xlabel('time days after first positive')
    plt.axhline(np.log10(1),ls='--',color='k')

plt.xlabel('time days after t_det')
plt.ylim([np.log10(0.02),10])
plt.xlim([-5,70])
plt.ylabel('viral load log10(copies/mL)')
plt.tight_layout()
plt.savefig('figures/ode_sim.pdf',dpi=600)


#sensitivity analysis
plt.figure(figsize=(10,5))
corrz=[[],[]]
for v in range(5):
    plt.subplot(2,5,1+v)
    p_sense=table1.iloc[:,v+1]
    
    if v==1:
        p_sense=p_sense*100
        
    plt.scatter(p_sense,VL_peak)
    if v==0:
        plt.ylabel('peak VL')
    corrz[0].append(st.pearsonr(table1.iloc[:,v+1],VL_peak)[0])
    
    plt.subplot(2,5,6+v)
    plt.scatter(p_sense,VL_setpt)
    plt.xlabel(pzn[v])
    if v==0:
        plt.ylabel('set point VL')
    corrz[1].append(st.pearsonr(table1.iloc[:,v+1],VL_setpt)[0])

plt.tight_layout()
plt.savefig('figures/ode_corr.pdf',dpi=600)

#plot correlations from sensitivity analysis
plt.bar(np.arange(5),corrz[0],width=0.2,label='peakVL')
plt.bar(np.arange(5)+0.2,corrz[1],width=0.2,label='setptVL')
plt.legend()
plt.ylabel('pearson correlation coeff')
plt.xticks(np.arange(5),pzn)
plt.ylim([-1,1])
plt.tight_layout()
plt.savefig('figures/ode_corr_BAR.pdf',dpi=600)

#now plot all the fits
#plot model fitting results
nx=10; ny=5;
fig,axarr=plt.subplots(ny,nx,sharey=True,sharex=True,
                       figsize=(nx*1.25,ny*1.25),
                       gridspec_kw={'wspace':0.1,'hspace':0.1})

sim_t=np.linspace(0,50*7,1e4);

#plot fits together
pptind=0
first_pos=[]
for ppt in params['id']:
    
    fitdat=RV217[RV217['ID']==ppt] #loop through each individual    

    data_t=np.array(fitdat['days']) #rescale
    data_V=np.array(fitdat['log10VL'])
    
    first_pos.append(data_V[0])
    idz,initT,dS,lBt0,dI,n=np.array(table1[table1['id']==ppt])[0]

    Bt0 = 10**(lBt0); #infectivity rate [uL/cells-day], re-exponentiated

    S_0 = aS/dS #assumes T cells at equilibrium
    V_0 = V0/1e3 #convert to uL
    P_0 = V_0*tau*g/p #start with a single infected cell

    sol=spi.odeint(ode_sim, [S_0,P_0,V_0], sim_t, (aS,dS,Bt0,tau,dI,n,p,g), mxstep=10000) 
    sim_V=np.log10(sol[:,2]*1e3) #viral load as usual, reconvert to copies/mL

    ax=axarr[int(pptind/nx)][pptind%nx]
    ax.scatter(data_t/7,data_V,c='gray',marker='o',lw=0,s=45,alpha=0.8)
    ax.plot((sim_t-initT)/7,sim_V,lw=2,color=cz_list[pptind],label=ppt)    
    ax.annotate(str(ppt),[-8,7.7],fontsize=8)
    pptind+=1
    #plt.axvline(1,color='k',ls='--')

while pptind<nx*ny:
    axarr[int(pptind/nx)][pptind%nx].remove()
    pptind+=1
    
plt.ylim([1,9])
plt.yticks(range(0,11,2))
plt.savefig('figures/FabMLEfits_fullX.pdf',dpi=600)

plt.xlim([-10,10])
plt.xticks(range(-10,11,5))
plt.savefig('figures/FabMLEfits.pdf',dpi=600)

#deterministic time vs first positive viral load
plt.scatter(first_pos,table1['initT_mode'],c=cz_list)
plt.xlabel('first positive viral load log10(copies/mL)')
plt.ylabel('deterministic time ($t_{det}$)')
plt.tight_layout()
plt.savefig('figures/firstpos_vs_tdet.pdf',dpi=600)

#correlate parameters with multiple infections!
single = [40512,20368,40353,20511,30190,40061,10066,20314,20225,20631,40511,40250,40257,40007,40577,40168,30924,20509,20507,40231,30112,10428,40094,40265]
multi =[10463,10220,30812,40100,40123,20502,30124,40436,40363]

plt.figure(figsize=(10,2))

sl=[[],[],[],[],[]]
ml=[[],[],[],[],[]]
for ppt in table1['id']:
    dat=np.array(table1[table1['id']==ppt])
    for pind in range(5):
        if ppt in single:
            sl[pind].append(dat[0,pind+1])
        if ppt in multi:
            ml[pind].append(dat[0,pind+1])
    
for pind in range(5):
    plt.subplot(151+pind)
    
    plt.boxplot([sl[pind],ml[pind]],widths=0.7)
    Ut,Up=st.mannwhitneyu(sl[pind],ml[pind])
    plt.xticks([1,2],['1','2+'])
    plt.xlim([0,3])
    #plt.ylabel(pzn[pind])
    plt.title(pzn[pind]+', p='+str(round(Up,2)),fontsize=10)
    
    if pind==2:
        plt.xlabel('# founders')
    
plt.tight_layout()
plt.savefig('figures/founders.pdf',dpi=600)


#compare with BEAST estimates
beast_est=pd.read_csv('data/timingclicks.csv')

plt.figure(figsize=(7,6))

dcomp=[]
for ppt in table1['id']:
    dat=np.array(table1[table1['id']==ppt])

    tdet=-dat[0,1]
    
    if ppt in list(beast_est['ID']):
        plt.scatter(beast_est['daysprefirstpos'][beast_est['ID']==ppt],tdet,
                    label=str(ppt),s=100,alpha=0.8,marker='D')
        dcomp.append([float(beast_est['daysprefirstpos'][beast_est['ID']==ppt]),tdet])
        
plt.xlabel('BEAST estimate (days relative to first positive)')
plt.ylabel('VMEM estimate (days relative to first positive)')
plt.xlim([-50,3])
plt.ylim([-50,3])
xy=np.arange(-60,3)
plt.plot(xy,xy,ls='--',color='gray')
plt.plot(xy,np.zeros(len(xy)),ls='--',color='gray')
plt.plot(np.zeros(len(xy)),xy,ls='--',color='gray')
plt.legend(bbox_to_anchor=(1.1,1))

#plt.tight_layout()
plt.savefig('figures/dcomp.pdf',dpi=600)

print(st.pearsonr(np.array(dcomp)[:,0],np.array(dcomp)[:,1]))
print(st.spearmanr(np.array(dcomp)[:,0],np.array(dcomp)[:,1]))

plt.boxplot(np.array(dcomp)[:,0]-np.array(dcomp)[:,1])



