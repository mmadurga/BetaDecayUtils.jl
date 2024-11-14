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

`nPenetrability(x,mass::Vector,Lorb)`

calculates the neutron penetrability p(x,Lorb).

`x` is the excitation energy above `Sₙ`, `Lorb` is the neutron angular momentum

`mass[1]` is the recoil, `mass[2]` is the neutron mass.


###

### Gaussian fit of gamma photopeak

`gammafit(data,xlow,xhigh,param::Vector,n=nothing,lowerbounds=nothing,upperbounds=nothing)`

Fit a gamma line in a histogram using a gaussian distribution with linear background. 
Returns a three element list containing: parameters from fit; standard error of the parameters; fitting function f(x,param)

Current version requires `LsqFit` to be loaded

`data`:           histogram in the `[energy counts]` format

`xlow`:           low energy cut for the fit

`xhigh`:          high energy cut for the fit

`param`:          initial fit parameters vector (must be 3*n+1). Format: `[background-constant,background-slope,area1,centroid1,sigma1,area2,centroid2,sigma2,...]`

`n`:              number of peaks to fit (default=1)

`lowerbounds`:    lower parameter bounds (optional). same format as param

`upperbounds`:    upper parameter bounds (optional). same format as param

Example: fit a gamma line at 988 keV

`xlow,xhigh,param=980,1000,[600,0.05,1000,988,0.8]`

`p,s,f = gammafit(data,xlow,xhigh,param)`

to print the result of the fit parameters

`for (i,val) in enumerate(p)`
    `println("P$i = ",val,"($(s[i]))")`
`end`

to plot the result of the fit use the returned function f with the optimized parameters p

`plot(e->f(e,p),xlow,xhigh)`

###