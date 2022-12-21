# M2AS---Covid-19

This repository contains two files:

# contacts.m

This script is used to estimate the reduction of the infection rate after isolating a given proportion v of the population in epidemiological models derived from SIR models (Kermack-McKendrick). 

It generates a batch of simulations which determine the number of infections in a synthetic population of N individuals (with I infectious individuals), with and without isolation. Contacts are randomly created between individuals as edges linking two vertices (individuals). The number of infectious contacts is then computed (number of edges linking exactly one susceptible and one infected individual). For each simulation, the number of infectious contact is computed with and without isolation, and the values obtained are represented on a figure. The figure obtained illustrate that the number of infectious contacts with isolation corresponds to   the number of infectious contacts without isolation multiplied by (1-v)^2.

# contactsIllustration.m
