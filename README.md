---
author: Nicola F. Müller
level: Professional
title: AIM in StarBeast2 v0.0.1 Tutorial
subtitle: Inferring species trees and gene flow from multiple loci
beastversion: 2.5.0
tracerversion: 1.6.0
figtreeversion: 1.4.2
---


# Background



----

# Programs used in this Exercise 

### BEAST2 - Bayesian Evolutionary Analysis Sampling Trees 2

BEAST2 ([http://www.beast2.org](http://www.beast2.org)) is a free software package for Bayesian evolutionary analysis of molecular sequences using MCMC and strictly oriented toward inference using rooted, time-measured phylogenetic trees. This tutorial is written for BEAST v{{ page.beastversion }} {% cite BEAST2book2014 --file MASCOT-Tutorial/master-refs %}. 


### BEAUti2 - Bayesian Evolutionary Analysis Utility

BEAUti2 is a graphical user interface tool for generating BEAST2 XML configuration files.

Both BEAST2 and BEAUti2 are Java programs, which means that the exact same code runs on all platforms. For us it simply means that the interface will be the same on all platforms. The screenshots used in this tutorial are taken on a Mac OS X computer; however, both programs will have the same layout and functionality on both Windows and Linux. BEAUti2 is provided as a part of the BEAST2 package so you do not need to install it separately.

### TreeAnnotator

TreeAnnotator is used to summarise the posterior sample of trees to produce a maximum clade credibility tree. It can also be used to summarise and visualise the posterior estimates of other tree parameters (e.g. node height).

TreeAnnotator is provided as a part of the BEAST2 package so you do not need to install it separately.


### Tracer

Tracer ([http://tree.bio.ed.ac.uk/software/tracer](http://tree.bio.ed.ac.uk/software/tracer)) is used to summarise the posterior estimates of the various parameters sampled by the Markov Chain. This program can be used for visual inspection and to assess convergence. It helps to quickly view median estimates and 95% highest posterior density intervals of the parameters, and calculates the effective sample sizes (ESS) of parameters. It can also be used to investigate potential parameter correlations. We will be using Tracer v{{ page.tracerversion }}.


### FigTree

FigTree ([http://tree.bio.ed.ac.uk/software/figtree](http://tree.bio.ed.ac.uk/software/figtree)) is a program for viewing trees and producing publication-quality figures. It can interpret the node-annotations created on the summary trees by TreeAnnotator, allowing the user to display node-based statistics (e.g. posterior probabilities). We will be using FigTree v{{ page.figtreeversion }}.

### R

R.....
R is needed to further analyse and plot the output of the AIM analysis

----

# Practical: Joint inference of the species history and gene flow of 5 anopheles mosquitos by using AIM and coupled MCMC

{% cite Mueller348391 --file AIM-Tutorial/master-refs.bib %}

In this tutorial, we wil infer the species history of 5 different anophoeles mosquitos species by using AIM, which is short for *Approximate Isolation with Migration*. AIM is part of the StarBeast2 package. This model can be used when we want to jointly infer the species histories for multiple loci and gene flow between extant and ancestral species 

The aim is to:

-  Learn how to jointly infer the species tree and gene flow of multiple loci from multiple species
-  Get to know how to choose the set-up of such an analysis
-  Learn how to process the output of an AIM analysis

## About the data

The data was previously used to infer the species history of the anopheles gambie comples in {% cite fontaine2015extensive --file AIM-Tutorial/master-refs.bib %}. It is comprised of 27 loci, each of a length of about 1000 bp, from the left arm of the 3rd chromosome. 



## Setting up an analysis in BEAUti

### Download StarBeast2 and CoupledMCMC
First, we have to download the packages StarBeast2 and CoupledMCMC by using the BEAUTi package manager. Go to _File >> Manage Packages_ and download the package and CoupledMCMC StarBeast2. 

<figure>
	<a id="fig:example1"></a>
	<img style="width:50%;" src="figures/StarBeastDownload.png" alt="">
	<figcaption>Figure 1: Download the StarBeast2 and CoupledMCMC packages.</figcaption>
</figure>



### Loading the template
Next, we have to load the BEAUTi template from _File_, select _Template >> AIM_.


### Loading the different loci

The sequences for the different loci can be found in the _data_ folder name can be either drag and dropped into BEAUti. To speed up the setup later, we can press _Link Site Models_ and _Link Clock Models_


### Get species corresponding to the different individuals (Taxon sets)

To assign the different individuals to different species, press the _Guess_ button. Next, use everything before first and press the _OK_ button.

<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/TaxonSet.png" alt="">
	<figcaption>Figure 3: Guess the species of each individual.</figcaption>
</figure>

### Specify the Site Model (Site Model)

Since we Linked all the Site Models of the different loci together when loading the sequence data, we only have to set up the site models once. We will be using an HKY model that allows for different relative rates of transversions and transitions. Additionally, we should make sure that the _estimate_ button for the Substitution rates is klicked to allow for rate variation across different loci.  To reduce the number of parameters we have to estimate, we can set Frequencies to Empirical. Next, we go back to the _Partitions_ field and press _Unlink Site Models_. Now each loci will have it's own site model. 

<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/SiteModel.png" alt="">
	<figcaption>Figure 4: Set the site model.</figcaption>
</figure>


### Set the clock model (Clock Model)

Since we have all sequences sampled in the present and no calibration, we have to information to estimate the clock rate.  


### Specify the priors (Priors)
Now, we need to set the priors of the effective population sizes and the migration rates. Next, we can change the prior to a Log Normal prior with M=0 and S=1. Since we have only a few samples per location, meaning little information about the different effective population sizes, we will need an informative prior.
Next, we have to set the dimension of the migration rate parameter. The exponential distribution as a prior on the migration rate puts much weight on lower values while not prohibiting larger ones. For migration rates, a prior that prohibits too large values while not greatly distinguishing between very small and very very small values (such as the inverse uniform) is generally a good choice.
Next, we have to set a prior for the clock rate. Since we only have a narrow time window of less than a year and only 24 sequences, there isn't much information in the data about the clock rate. We have however a good idea about it for Influenza A/H3N2 Hemagglutinin. We can therefore set the prior to be normally distributed around 0.005 substitution per site and year with a variance of 0.0001. (At this point we could also just fix the rate)

<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/Priors.png" alt="">
	<figcaption>Figure 6: Set up of the prior distributions.</figcaption>
</figure>


### Specify the MCMC chain length (MCMC)

Here we can set the length of the MCMC chain and after how many iterations the parameter and trees a logged. For this dataset, 2 million iterations should be sufficient. In order to have enough samples but not create too large files, we can set the logEvery to 2500, so we have 801 samples overall. Next, we have to save the `*.xml` file under _File >> Save as_.

<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/MCMC.png" alt="">
	<figcaption>Figure 7: save the *.xml.</figcaption>
</figure>

### Set up the xml to run two chains
In order to setup the analysis to run with coupled MCMC, we have to open the  `*.xml` and change one line in the xml.
To do so, go to the line with:
```
<run id="mcmc" spec="MCMC" chainLength="200000000" numInitializationAttempts="10">
```
To have a run with coupled MCMC, we have to replace that one line with:
```
<run id="mcmc" spec="beast.coupledMCMC.CoupledMCMC" logHeatedChains="true" chainLength="100000000" storeEvery="1000000" deltaTemperature="0.025" chains="2" resampleEvery="10000">
```
* `logHeatedChains="true"` also logs the log files of the heated chains if true.
* `chainLength="100000000"` defines for how many iterations the chains is run
* `deltaTemperature="0.025"` defines the temperature difference between the chain *n* and chain *n-1*.
* `chains="2"` defines the number of parallel chains that are run. The first chain is the one that explores the posterior just like a normal MCMC chain. All other chains are what's called *heated*. This means that MCMC moves of those chains have a higher probability of being accepted. While these heated chains don't explore the posterior properly, they can be used to propose new states to the one cold chain.   




Here we can set the length of the MCMC chain and after how many iterations the parameter and trees a logged. For this dataset, 2 million iterations should be sufficient. In order to have enough samples but not create too large files, we can set the logEvery to 2500, so we have 801 samples overall. Next, we have to save the `*.xml` file under _File >> Save as_.

<figure>
<a id="fig:example1"></a>
<img style="width:70%;" src="figures/MCMC.png" alt="">
<figcaption>Figure 7: save the *.xml.</figcaption>
</figure>


### Run the Analysis using BEAST2
Run the `*.xml` using BEAST2 or use finished runs from the *precooked-runs* folder. The analysis should take about 6 to 7 minutes. 

### Analyse the log file using Tracer

First, we can open the `*.log` file in tracer to check if the MCMC has converged. The ESS value should be above 200 for almost all values and especially for the posterior estimates. The burnin taken by Tracer is 10%, but for this analysis 1% is enough.

<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/LogPosterior.png" alt="">
	<figcaption>Figure 8: Check if the posterior converged.</figcaption>
</figure>



Next, we can have a look at the inferred effective population sizes. New York is inferred to have the largest effective population size before Hong Kong and New Zealand. This tells us that two lineages that are in the New Zealand are expected to coalesce quicker than two lineages in Hong Kong or New York.

<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/LogNe.png" alt="">
	<figcaption>Figure 9: Compare the different inferred effective population sizes.</figcaption>
</figure>

In this example, we have relatively little information about the effective population sizes of each location. This can lead to estimates that are greatly informed by the prior. Additionally, there can be great differences between median and mean estimates. The median estimates are generally more reliable since they are less influence by extreme values. 

<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/MeanMedian.png" alt="">
	<figcaption>Figure 10: Differences between Mean and Meadian estimates.</figcaption>
</figure>

We can then look at the inferred migration rates. The migration rates have the label b_migration.*, meaning that they are backwards in time migration rates. The highest rates are from New York to Hong Kong. Because they are backwards in time migration rates, this means that lineages from New York are inferred to be likely from Hong Kong if we're going backwards in time. In the inferred phylogenies, we should therefore make the observation that lineages ancestral to samples from New York are inferred to be from the Hong Kong backwards. A very good explanaition on the differences between forward and backwards migration rates can be found in the following blog post by Peter Beerli [http://popgen.sc.fsu.edu/Migrate/Blog/Entries/2013/3/22_forward-backward_migration_rates.html](http://popgen.sc.fsu.edu/Migrate/Blog/Entries/2013/3/22_forward-backward_migration_rates.html)

<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/LogMigration.png" alt="">
	<figcaption>Figure 11: Compare the inferrred migration rates.</figcaption>
</figure>

### Make the MCC tree using TreeAnnotator
Next, we want to summarize the trees. This we can do using treeAnnotator. Open the programm and then set the options as below. You have to specify the _Burnin precentage_, the _Node heights_, _Input Tree File_ and the _Output File_ after clicking _Run_ the programm should summarize the trees.

<figure>
	<a id="fig:example1"></a>
	<img style="width:50%;" src="figures/TreeAnnotator.png" alt="">
	<figcaption>Figure 12: Make the maximum clade credibility tree.</figcaption>
</figure>

### Check the MCC tree using FigTree
We can now open the MCC tree using FigTree. The output contains several things. Each node has several traits. Among them are those called Hong_Kong, New_York and New_Zealand. The value of those traits is the probability of that node being in that location as inferred using MASCOT. 


<figure>
	<a id="fig:example1"></a>
	<img style="width:100%;" src="figures/HongKongLabels.png" alt="">
	<img style="width:100%;" src="figures/NewZealandLabels.png" alt="">
	<img style="width:100%;" src="figures/NewYorkLabels.png" alt="">
	<figcaption>Figure 13: Compare the inferred node probabilities.</figcaption>
</figure>

We can now check if lineages ancestral to samples from New York are actually inferred to be from Hong Kong, or the probability of the root being in any of the locations. It should here be mentioned that the inference of nodes being in a particular location makes some simplifying assumptions, such as that there are no other locations where lineages could have been. To actually visualize the colors of node, you can go to _Appearance >> Colour by_ and select *max*. This will color each node and the branch ancestral to this node by the location that was inferred to be most often the most likely location.

----

# Useful Links

- MASCOT source code: [https://github.com/nicfel/Mascot](https://github.com/nicfel/Mascot)
- [Bayesian Evolutionary Analysis with BEAST 2](http://www.beast2.org/book.html) {% cite BEAST2book2014 --file MASCOT-Tutorial/master-refs.bib %}
- BEAST 2 website and documentation: [http://www.beast2.org/](http://www.beast2.org/)
- Join the BEAST user discussion: [http://groups.google.com/group/beast-users](http://groups.google.com/group/beast-users) 

----

# Relevant References

{% bibliography --cited --file MASCOT-Tutorial/master-refs %}

