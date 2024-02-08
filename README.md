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

### Daughter Activity

`daughterActivity(x,A,λ)`

`A`: initial activity

`λ`: decay probability (ln2/T12)

###



### Granddaughter activity

`grandDaughterActivity(x,A,λ,μ)`

`A`: initial daughter activity

`λ`: daughter decay probability 

`μ`: granddaughter decay probability


### Halflife from BGT distribution


`calculateT12(z,Qᵦ,Eₓ,BGT)`

calculate halflife of the beta decay of an isotope given feedings to excited states

`Qᵦ`: β decay Q value in MeV

`Eₓ`: vector of daughter states relative to the ground state energy in MeV

`BGT`: vector of BGT values

`z`: parent Z 

### Partial branching ratios from BGT distribution

`calculateIb(z,Qᵦ,Eₓ,BGT)`

calculate branching ratios of the beta decay of an isotope given feedings to excited states

`Qᵦ`: β decay Q value in MeV

`Eₓ`: vector of daughter states relative to the ground state energy in MeV

`BGT`: vector of BGT values

`z`: parent Z 

###

### Log(ft) from partial branching ratios (Iᵦ)

`logftfromib(z,t₁₂,Qᵦ,Eₓ,Iᵦ)`

calculate logft of a given transition to an excitated state

`Z` of the parent, `Qᵦ` and `Eₓ` in MeV, `t₁₂` in seconds, `Iᵦ` absolute value


### Log(ft) from BGT distribution

`logftfrombgt(bgt)`

calculate the logft for a given `BGT` (not quenched)

###

### nPenetrability

`nPenetrability(x,mass::Vector,Lorb)`

calculates the neutron penetrability p(x,Lorb). Given a general reduced width γ one can calculate the partial width as

Γ = p(x,Lorb)

`x` is the excitation energy above `Sₙ`, `Lorb` is the neutron amngular momentum

`mass[1]` is the recoil, `mass[2]` is the neutron mass


###