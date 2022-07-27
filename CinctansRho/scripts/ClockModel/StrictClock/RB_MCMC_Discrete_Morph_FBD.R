
taxa <- readTaxonData("data/cincta_fas_fuzzy.tsv")
taxa
morpho_f <- readDiscreteCharacterData("data/feedingCharacters.nex")
morpho_nf <- readDiscreteCharacterData("data/nonfeedingCharacters.nex")
moves = VectorMoves()
monitors = VectorMonitors()

     n_taxa <- taxa.size()

     num_branches <- 2 * n_taxa - 2


     moves = VectorMoves()

     monitors = VectorMonitors()

      # Diversification Rates based on Echinodermata
      speciation_rate ~ dnExponential(1.471);
      # NOTE: If it gets stuck in this script, then set origination & extinction to 1.0
      moves.append(mvScale(speciation_rate, lambda=0.01, weight=5));
      moves.append(mvScale(speciation_rate, lambda=0.10, weight=3));
      moves.append(mvScale(speciation_rate, lambda=1.00, weight=1));

      turnover ~ dnUnif(0.9, 1.05);
      moves.append(mvSlide(turnover, delta=0.01, weight=5));
      moves.append(mvSlide(turnover, delta=0.10, weight=3));
      moves.append(mvSlide(turnover, delta=1.00, weight=1));
      extinction_rate := turnover*speciation_rate;
      diversification := speciation_rate - extinction_rate;

      # old extinction stuff. We should not use this, as extinction should not be independent of origination!
      #extinction_rate ~ dnExponential(1.471);
      #moves.append(mvScale(extinction_rate, lambda=0.01, weight=5));
      #moves.append(mvScale(extinction_rate, lambda=0.10, weight=3));
      #moves.append(mvScale(extinction_rate, lambda=1.00, weight=1));
      #turnover := extinction_rate/speciation_rate;

      # Fossil Sampling Rates based on collection occupied by Echinodermata
      psi ~ dnExponential(3.892);
      completeness := psi/(extinction_rate+psi);
      moves.append(mvScale(psi, lambda=0.01, weight=5));
      moves.append(mvScale(psi, lambda=0.10, weight=3));
      moves.append(mvScale(psi, lambda=1.00, weight=1));

      # Proportional Taxon Sampling of Youngest Time Slice
      rho <- 0.506;	# 'extant' sampling.

     # Establish Basal Divergence Time
     origin_time ~ dnUnif(7.3, 12.11);
     moves.append(mvSlide(origin_time, delta=0.01, weight=5));
     moves.append(mvSlide(origin_time, delta=0.10, weight=3));
     moves.append(mvSlide(origin_time, delta=1.00, weight=1));


     fbd_dist = dnFBDP(originAge=origin_time, lambda=speciation_rate, mu=extinction_rate, psi=psi, rho=rho, timeline=timeline, taxa=taxa, condition="sampling")



     outgroup = clade("Ctenocystis");
     ingroup = clade("Gyrocystis_platessa","Gyrocystis_testudiformis","Gyrocystis_cruzae","Gyrocystis_badulesiensis","Gyrocystis_erecta","Progyrocystis_disjuncta","Protocinctus_mansillaensis","Elliptocinctus_barrandei","Elliptocinctus_vizcainoi","Sucocystis_theronensis","Sucocystis_bretoni","Lignanicystis_barriosensis","Undatacinctus_undata","Sucocystis_acrofera","Undatacinctus_quadricornuta","Undatacinctus_melendezi","Asturicystis_jaekeli","Sotocinctus_ubaghsi","Trochocystites_bohemicus","Trochocystoides_parvus","Ludwigicinctus_truncatus","Graciacystis_ambigua","Asturicystis_havliceki","Nelegerocystis_ivantzovi","Rozanovicystis_triangularis","Davidocinctus_pembrokensis");

 constraints = v(ingroup)


     fbd_tree ~ dnConstrainedTopology(fbd_dist, constraints=constraints)


     moves.append(mvFNPR(fbd_tree, weight=15.0))

     moves.append(mvCollapseExpandFossilBranch(fbd_tree, origin_time, weight=6.0))

     moves.append(mvNodeTimeSlideUniform(fbd_tree, weight=40.0))

     moves.append(mvRootTimeSlideUniform(fbd_tree, origin_time, weight=5.0))



     intervals = readDataDelimitedFile(file="data/cincta_fas_fuzzy.tsv", header=true)




 # Setup the fossil tip sampling #

 # Use a for loop to create a uniform distribution on the occurence time for each fossil #

 # The boundaries of the uniform distribution are specified in the tsv file #

 fossils = fbd_tree.getFossils()
 for(i in 1:fossils.size())
 {
     t[i] := tmrca(fbd_tree, clade(fossils[i]))

     a_i = fossils[i].getMinAge()
     b_i = fossils[i].getMaxAge()

     F[i] ~ dnUniform(t[i] - b_i, t[i] - a_i)
     F[i].clamp( 0 )
 }

 # Add a move to sample the fossil times #
 moves.append( mvFossilTimeSlideUniform(fbd_tree, origin_time, weight=5.0) )


     num_samp_anc := fbd_tree.numSampledAncestors()

     pruned_fbd_tree := fnPruneTree(fbd_tree, prune=v("Asturicystis_havliceki","Nelegerocystis_ivantzovi","Rozanovicystis_triangularis","Davidocinctus_pembrokensis"))


     clock_morpho ~ dnExponential(1.0)

     moves.append( mvScale(clock_morpho, lambda=0.01, weight=4.0) )
     moves.append( mvScale(clock_morpho, lambda=0.1,  weight=4.0) )
     moves.append( mvScale(clock_morpho, lambda=1,    weight=4.0) )



     alpha_morpho ~ dnUniform( 0, 1E6 )

     rates_morpho := fnDiscretizeGamma( alpha_morpho, alpha_morpho, 4 )

     #Moves on the parameters of the Gamma distribution.

     moves.append(mvScale(alpha_morpho, lambda=1, weight=2.0))





     n_max_states <- 3
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
                                         siteRates=rates_morpho,
                                         branchRates=clock_morpho,
                                         type="Standard")

             # attach the data
     	    m_morph[idx].clamp(morpho_f_bystate[i])

             # increment counter
             idx = idx + 1
     idx
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
                                         siteRates=rates_morpho,
                                         branchRates=clock_morpho,
                                         type="Standard")

             # attach the data
     	    m_morph[idx].clamp(morpho_nf_bystate[i])

             # increment counter
             idx = idx + 1
     idx
     }
     }



     mymodel = model(fbd_tree)



     monitors.append(mnModel(filename="output/cinc6_dated.log", printgen=10))



     monitors.append(mnFile(filename="output/cinc6_dated.trees", printgen=10, pruned_fbd_tree))



     monitors.append(mnScreen(printgen=10, num_samp_anc, origin_time))



     mymcmc = mcmc(mymodel, monitors, moves)


     mymcmc.run(generations=1000000)



     q()
