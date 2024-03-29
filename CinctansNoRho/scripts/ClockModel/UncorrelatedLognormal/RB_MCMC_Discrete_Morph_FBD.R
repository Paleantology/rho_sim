## ----global_options, eval = TRUE, include=TRUE---------------------------
name = "mk_gamma_FBD"

moves = VectorMoves()
monitors = VectorMonitors()

## ---- include=TRUE, eval = TRUE------------------------------------------
morpho_f <- readDiscreteCharacterData("data/feedingCharacters.nex")
morpho_nf <- readDiscreteCharacterData("data/nonFeedingCharacters.nex")
taxa <- readTaxonData("data/cincta_fas_fuzzy.tsv")


## ---- include=TRUE, eval = TRUE------------------------------------------
    taxa <- morpho_f.names()
    n_taxa <- taxa.size()
    num_branches <- 2 * n_taxa - 2

    source("scripts/ClockModel/StrictClock/Basic_FBD_model.R")


## ---- include=TRUE, eval = TRUE------------------------------------------
    alpha_morpho ~ dnUniform( 0, 1E6 )
    rates_morpho := fnDiscretizeGamma( alpha_morpho, alpha_morpho, 4 )
    #Moves on the parameters to the Gamma distribution.
    moves.append(mvScale(alpha_morpho, lambda=1, weight=2.0))

    ucln_mean ~ dnExponential(2.0)
    ucln_sigma ~ dnExponential(3.0)
    ucln_var := ucln_sigma * ucln_sigma
    ucln_mu := ln(ucln_mean) - (ucln_var * 0.5)
    moves.append( mvScale(ucln_mean, lambda=1.0, tune=true, weight=4.0))
    moves.append( mvScale(ucln_sigma, lambda=0.5, tune=true, weight=4.0))
    for(i in 1:num_branches){
       branch_rates[i] ~ dnLnorm(ucln_mu, ucln_sigma)
       moves.append( mvScale(branch_rates[i], lambda=1, tune=true, weight=2.))
    }

## ---- include=TRUE, eval = TRUE------------------------------------------
alpha_morpho ~ dnUniform( 0, 1E6 )
alpha_morpho2 ~ dnUniform( 0, 1E6 )

rates_morpho := fnDiscretizeGamma( alpha_morpho, alpha_morpho, 4 )
rates_morpho2 := fnDiscretizeGamma( alpha_morpho2, alpha_morpho2, 4 )

#Moves on the parameters to the Gamma distribution.
moves.append(mvScale(alpha_morpho, lambda=1, weight=2.0))
moves.append(mvScale(alpha_morpho2, lambda=1, weight=2.0))


## ---- include=TRUE, eval = TRUE------------------------------------------
n_max_states <- 2
idx = 1
morpho_f_bystate <- morpho_f.setNumStatesVector()
for (i in 1:n_max_states) {
nc = morpho_f_bystate[i].nchar()
# for non-empty character blocks
if (nc > 0) {
    # make i-by-i rate matrix
    q[idx] <- fnJC(i)
# create model of evolution for the character block
    m_morph[idx] ~ dnPhyloCTMC( tree=fbd_tree,
                                Q=q[idx],
                                nSites=nc,
                                branchRates=branch_rates,
                                siteRates=rates_morpho,
                                type="Standard")

    # attach the data
  m_morph[idx].clamp(morpho_f_bystate[i])

    # increment counter
    idx = idx + 1

}
}

n_max_states <- 4
morpho_nf_bystate <- morpho_nf.setNumStatesVector()
for (i in 1:n_max_states) {
nc = morpho_nf_bystate[i].nchar()
# for non-empty character blocks
if (nc > 0) {
    # make i-by-i rate matrix
    q[idx] <- fnJC(i)
# create model of evolution for the character block
    m_morph[idx] ~ dnPhyloCTMC( tree=fbd_tree,
                                Q=q[idx],
                                nSites=nc,
                                branchRates=branch_rates,
                                siteRates=rates_morpho2,
                                type="Standard")

    # attach the data
  m_morph[idx].clamp(morpho_nf_bystate[i])
    # increment counter
    idx = idx + 1

}
}
## ---- include=TRUE, eval = TRUE------------------------------------------
    mymodel = model(fbd_tree)


## ---- include=TRUE, eval = TRUE------------------------------------------
    monitors.append(mnModel(filename="output/" +name + "log", printgen=10))


## ---- include=TRUE, eval = TRUE------------------------------------------
    monitors.append(mnFile(filename="output/" + name + "trees", printgen=10, fbd_tree))


## ---- include=TRUE, eval = TRUE------------------------------------------
    monitors.append(mnScreen(printgen=100))


## ---- include=TRUE, eval = TRUE------------------------------------------
    mymcmc = mcmc(mymodel, monitors, moves, nruns=2, combine="mixed")

#ss_analysis = powerPosterior(mymodel, monitors, moves, "output/ucln/" + name + "/ss", cats=20, alpha=0.3)
#ss_analysis.burnin(generations=1000,tuningInterval=100)
#ss_analysis.run(generations=50000)

#ss = steppingStoneSampler("output/" + name + "/ucln_ss", "power", "likelihood", TAB)
#ss.marginal()
### ---- include=TRUE, eval = TRUE------------------------------------------
    mymcmc.run(generations=100000, tuningInterval=200)


## ---- include=TRUE, eval = TRUE------------------------------------------
    q()


## ------------------------------------------------------------------------
