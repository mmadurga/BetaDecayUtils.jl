# BetaDecayUtils

[![Build Status](https://github.com/mmadurga/BetaDecayUtils.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mmadurga/BetaDecayUtils.jl/actions/workflows/CI.yml?query=branch%3Amain)

## β decay utilities for low energy nuclear physics.

Small collection of functions to calculate typical beta-decay observables and other properties. 

Log(f) calculation adapted from Juhani Kantele's "Handbook of Nuclear Spectrometry", Academic Press Limited, London.

Neutron penetrability calculation adapted from "Theoretical Nuclear Physics" by Blatt and Weisskopf.

### Coulomb Shift

`ecoulomb(Z,N)`

calculate the coulomb correction for a shell model calculation of isotope (Z,N)

###

### Beta Decay Q value from shell model

`calculateqbetashellmodel(groundstateparent,groundstatechild,zparent,aparent,zchild,achild)`

calculate the beta-decay Q value using shell model calculation of the ground state energies

groundstateparent:absolute energy of the parent ground state in MeV
groundstatechild:absolute energy of the child ground state in MeV 
zparent: Z of the parent
aparent: A of the parent
zchild: Z of the child
aparent: A of the child

###


### Child Activity

`childActivity(x,A,λ)`

`A`: initial activity

`λ`: decay probability (ln2/T12)

###

### Grandchild activity

`grandChildActivity(x,A,λ,μ)`

`A`: initial daughter activity

`λ`: child decay probability 

`μ`: grandchild decay probability

###

### chainActivity(x,A,λ,n)

`chainActivity(x,A,λ,n)`

activity of the n'th member of a decay chain.

`A`:initial parent nuclei. If initial `activity` is desired as input use `activity/λ[1]`

`λ`:vector containing  all decay probabilities in the chain

`n`:decay curve of n'th member

###

### Logarithm of the Fermi function

`logf(z,Qᵦ,Eₓ)`

Calculate the log10 of the Fermi function for allowed beta decay

`z`: atomic number of the parent

`Qᵦ`: β decay Q value in MeV

`Eₓ`: daughter state energy relative to the ground state energy in MeV

###

### Halflife from BGT distribution


`calculateT12(z,Qᵦ,Eₓ,BGT)`

calculate halflife of the beta decay of an isotope given feedings to excited states

`Qᵦ`: β decay Q value in MeV

`Eₓ`: vector of daughter states relative to the ground state energy in MeV

`BGT`: vector of BGT values

`z`: atomic number of the parent

### Partial branching ratios from BGT distribution

`calculateIb(z,Qᵦ,Eₓ,BGT)`

calculate branching ratios of the beta decay of an isotope given feedings to excited states

`Qᵦ`: β decay Q value in MeV

`Eₓ`: vector of daughter states relative to the ground state energy in MeV

`BGT`: vector of BGT values

`z`: atomic number of the parent

###

### Log(ft) from partial branching ratios (Iᵦ)

`logftfromib(z,t₁₂,Qᵦ,Eₓ,Iᵦ)`

calculate logft of a given transition to an excitated state

`z`: atomic number of the parent

`t₁₂`: decay halflife

`Qᵦ`: β decay Q value in MeV

`Eₓ`: energy of daughter state relative to the ground state energy in MeV

`Iᵦ`: partial branching values


### Log(ft) from BGT distribution

`logftfrombgt(bgt)`

calculate the logft for a given `BGT` (not quenched)

###

### Neutron penetrability as a function of the neutron angular momentum

`nPenetrability(x,mass::Vector,Lorb)`

calculates the neutron penetrability p(x,Lorb).

`x` is the excitation energy above `Sₙ`, `Lorb` is the neutron angular momentum

`mass[1]` is the recoil, `mass[2]` is the neutron mass.


###