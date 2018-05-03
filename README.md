# JAGS-Source-Tracker

_Author: Zeph Turner, under advisement of Dr. Kim Roth, Juniata College_

Estimating sources of metagenomic data using JAGS through R.

Included files:

- jsg_functions: A set of R functions that allow the use of a Bayesian hierarchical model of metagenomic source mixing to estimate the sources of an input sink sample. Includes functions for creating datasets of simulated sources mixed in a sample, estimating the sources of a sample from an OTU table and metadata map, and visualizing results.
- JSG_tutorial: An RMD file and the HTML output it creates describing how to call the functions in jsg_functions. View the HTML file for a quick overview of how to use jsg_functions.
- JAGS_source_tracker_theory: A paper describing the theory behind the JAGS source tracker application.
- example_biom_data and its contents: A folder containing an example .biom file to show how to manipulate .biom files to be used with jsg_functions.

This project was based on the previous work of Knights et al., described in their 2011 paper, [Bayesian community-wide culture-independent microbial source tracking](https://www.ncbi.nlm.nih.gov/pubmed/21765408). The GitHub page for SourceTracker2, a different piece of software that implements the same metagenomic mix model implemented by JAGS Source Tracker, can be found [here](https://github.com/biota/sourcetracker2). 

The code here uses [Just Another Gibbs Sampler (JAGS)](http://mcmc-jags.sourceforge.net/) to implement the metagenomic mixing model for data analysis.