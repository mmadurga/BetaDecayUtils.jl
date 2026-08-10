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

`calculateqbetashellmodel(groundstateparent,groundstatechild,zparent,nparent,zchild,nchild)`

calculate the beta-decay Q value using shell model calculation of the ground state energies

`groundstateparent`:absolute energy of the parent ground state in MeV

`groundstatechild`:absolute energy of the child ground state in MeV 

`zparent`: Z of the parent

`nparent`: N of the parent

`zchild`: Z of the child

`nparent`: N of the child

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

### BGTs from partial branching ratios

`calculateBGT(z,Qᵦ,T₁₂,Eₓ,Iᵦ)`

calculate beta decay strength of an isotope given partial branching ratios to excited states

`z:` Atomic number of the parent

`Qᵦ`: β decay Q value in MeV

`T₁₂`: beta decay half-life

`Eₓ`: vector of daughter states relative to the ground state energy in MeV

`Iᵦ`: vector of Iᵦ values


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

`nPenetrability(x,AM1,AM2,Lorb; R=1.4*(AM1+AM2)^(1/3))`

calculates the neutron penetrability p(x,Lorb).

`x` is the neutron energy above `Sₙ`, `Lorb` is the neutron angular momentum, and `R` is the channel radius in fm.

`AM1` is the recoil, `AM2` is the neutron mass in au.


###

### E2 Electromagnetic decay halflife

`calculateBE2T12(A,Eᵧ,BE2,units=nothing)`

calculate the halflife of a E2 nuclear electromagnetic transition (does not include correction for conversion electrons)

`A`: nuclear mass number
`Eᵧ`: transition energy in MeV
`BE2`: transition probability in e2fm4 (for `BE2` in w.u. define units=`"wu"`)


###

### Single particle resonance width.

`Gamma(E,AM1,AM2,Lorb; R=1.4*(AM1+AM2)^(1/3), warn_validity=false)`

calculates the Bohr--Mottelson low-energy square-well neutron single-particle width

`Gamma = 2*P*(ħc)^2/(μ*c^2*R^2)*(2*Lorb-1)/(2*Lorb+1)`

in MeV. `E` is the neutron energy above Sₙ, `AM1` and `AM2` are the residue and neutron masses in au, `Lorb ≥ 1` is the neutron angular momentum, and `R` is the channel radius in fm. The approximation is derived for `kR << 1`; set `warn_validity=true` to warn when `kR ≥ 1`. It is a depth-independent single-particle estimate rather than an exact complex-pole width, and a physical partial width may require a spectroscopic factor.

The former `Gamma(E,V0,AM1,AM2,Lorb; R=...)` and `Gamma(E,V0,AM1,AM2,Lorb,R)` forms are deprecated compatibility wrappers. `V0` is ignored.

###
