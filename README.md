# M2AS---Covid-19

This repository contains two files:

## contacts.m

This script is used to estimate the reduction of the infection rate after isolating a given proportion v of the population in epidemiological models derived from SIR models (Kermack-McKendrick). 

This script simulates the isolation of a proportion v of the population on synthetic population of N individuals (with I infectious individuals), with and without isolation. It generates a batch of simulations which determine the number of infections in a population with and without isolation. Contacts are randomly created between individuals as edges linking two vertices (individuals). The number of infectious contacts is then computed (number of edges linking exactly one susceptible and one infected individual). For each simulation, the number of infectious contact is computed with and without isolation, and the values obtained are represented on a figure. The figure obtained illustrate that the number of infectious contacts with isolation corresponds to   the number of infectious contacts without isolation multiplied by (1-v)^2.

## contactsIllustration.m

This script generates a figure of contacts in order to illustrate the reduction of infectious contacts when a proportion v of the population is isolated. 

N individuals are created and are randomly assigned a susceptible (blue) or infectious (red) status. Contacts are randomly drawn as edges linking two individuals.

The left figures shows contacts without isolation. Infectious contacts are represented in red color, non-infectious contacts (two susceptible or two infectious) in dotted grey segments.

The right figures show the same graph after isolating some individuals (greyed individuals). Grey dashed edges represent contacts that were infectious without isolation, but that are discarded because of isolation.
