
To do
================

all:

- look into low dose intra-v/r challenge data in monkeys - we're in the window?
dan:
- fix python code with new monolix
- make it true python, add to github, more of pipeline?
- make a few plots
  - tdet vs first pos VL
  - t_sto vs t_det heat map
  - t_sto distribution
- quantify histogram estimates in table (mean/median/95%CI)
- double check on multiple founders (bimodal tdet)
- make table for best parameters for each individual
- compare results with Morgane! add to final figure

bryan:

- exclusion criterion
  - 2 histograms, peak, set-point
  - CV criterion
- estimate a few VL's from early APTIMA and add last negative

```
#python code for histogram exclusion by peak
pk=[]
for ppt in params['id']:
    fitdat=RV217[RV217['ID']==ppt] #loop through each individual    
    mxx=np.max(fitdat['log10VL'])
    pk.append(mxx)
    if mxx<5:
        print(ppt,mxx)
plt.hist(pk)
```

Fabian:

- last negative for censored data, refit in Monolix
- check on VL threshold for estimates, rerun Monolix with different Vdet =[0.01,0.1,1,50]
- Table for AIC and model selection
- make supplementary figure for identifiability by correlation


Old stuff:

- R scripts to exhaustively fit models while narrowing down parameters to find the optimal estimator of infection time?
- Dan will look a bit into diversity stuff and see if that can be incorporated easily.

# DONE Summer 2019

- All (Bryan) learn about SAEM.
- Bryan can explore the correlation problem.
- Any of us can use the Monolix based on Fabian's initial set up to play around.
- Dan will take the parameter estimates and start simulating data, both control and VRC01 cohorts.
- Dan will send back undersampled data with the aim of sending data that looks like real trial data with the infection time blinded.
- Fabian can start to try to do the estimation in the control and VRC01 cohorts?
- Bryan will think about, or maybe even start formulating an approach to uncertainty quantification?



# DONE Mid Jan 2019

- Fabian can try to set up the Monolix to fit the data? I think we’ve decided on 3 models?
    1. Basic with power law death rate
    2. My monkey model
    3. Monkey model + memory CD8 compartment ("smooth monkey" I’m calling it in my head)

# DONE Mid Dec 2018

- Bryan can update GitHub so that the data are only up to 12 weeks? And maybe drop any individual with fewer than 5 data points?
