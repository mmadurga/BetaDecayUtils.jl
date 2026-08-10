module BetaDecayUtils

export  calculateqbetashellmodel,
        logf,
        logftfromib,
        calculateT12,
        calculateT12fullFermi,
        childActivity, 
        grandChildActivity,
        chainActivity,
        nPenetrability,
        logftfrombgt,
        calculateIb,
        calculateBGT,
        ecoulomb,
        calculateBE2T12,
        Gamma

import LsqFit,SpecialFunctions,QuadGK

### β decay utilities

"""

ecoulomb(Z,N)

calculate the coulomb correction for a shell model calculation of isotope (Z,N)

"""
function ecoulomb(Z,N)
    A = Z + N; T=abs(Z-N)/2
    Rc = ℯ^(1.5/A)*A^(1/3)*(0.946-0.573*(2*T/A)^2)
    return 0.7*(Z*(Z-1)-0.76*(Z*(Z-1))^(2/3))/(Rc)
end

"""

calculateqbetashellmodel(groundstateparent,groundstatechild,zparent,nparent,zchild,nchild)

calculate the beta-decay Q value using shell model calculation of the ground state energies

groundstateparent:absolute energy of the parent ground state in MeV

groundstatechild:absolute energy of the child ground state in MeV 

zparent: Z of the parent

nparent: N of the parent

nchild: Z of the child

nparent: N of the child

"""
function calculateqbetashellmodel(groundstateparent,groundstatechild,zparent,nparent,zchild,nchild)
    return groundstateparent + ecoulomb(zparent,nparent) - groundstatechild - ecoulomb(zchild,nchild) + 0.78
end

"""

childActivity(x,A,λ)

A is initial activity), λ is the decay probability (ln2/T12)

"""
function childActivity(x,A,λ::Float64)
    return A*exp(-λ*x)  #A is initial activity, λ is the decay probability (ln2/T12)
end


"""

grandChildActivity(x,A,λ,μ)

A:initial activity, λ:child, μ:grandchild

"""
function grandChildActivity(x,A,λ::Float64,μ::Float64)
    return (μ)/(μ-λ)*A*(exp(-λ*x)-exp(-μ*x)) #A:initial activity, λ:child, μ:grandchild
end


"""

chainActivity(x,A,λ,n)

activity of the n'th member of a decay chain.

A:initial parent nuclei. If initial activity is desired as input use activity/λ[1]

λ:vector containing  all decay probabilities in the chain

n:decay curve of n'th member

"""
function chainActivity(t,A,λ,n)

    sum=0

    for i in 1:n

        product=λ[i]

        for j in 1:n if j!=i product*=λ[j]/(λ[j]-λ[i]);  end end
    
        sum+=product*exp(-λ[i]*t)

    end
    
    return A*sum

end


"""

logf(z,Qᵦ,Eₓ)

Calculate the log10 of the Fermi function for allowed beta decay

z: atomic number of the parent

Qᵦ: β decay Q value in MeV

Eₓ: daughter state energy relative to the ground state energy in MeV

"""
function logf(z,Qᵦ,Eₓ)  
    coeff = [ -17.2       7.9015    -2.54        0.28482;
    3.31368   -2.06273    0.703822   -0.075039;
   -0.364018   0.387961  -0.142528    0.016;
    0.0278071 -0.026519   0.0098854  -0.00113772
]  ;
zDaughter = z + 1
evalCoeff = [
coeff[1,1] + log(zDaughter) * coeff[1,2] + coeff[1,3]*log(zDaughter)^2. + coeff[1,4]*log(zDaughter)^3.,
coeff[2,1] + log(zDaughter) * coeff[2,2] + coeff[2,3]*log(zDaughter)^2. + coeff[2,4]*log(zDaughter)^3.,
coeff[3,1] + log(zDaughter) * coeff[3,2] + coeff[3,3]*log(zDaughter)^2. + coeff[3,4]*log(zDaughter)^3.,
coeff[4,1] + log(zDaughter) * coeff[4,2] + coeff[4,3]*log(zDaughter)^2. + coeff[4,4]*log(zDaughter)^3.  
]
βEp = (Qᵦ -  (Eₓ)) * 1000 #convert to keV

return evalCoeff[1] + evalCoeff[2]*log(βEp) + evalCoeff[3]*log(βEp)^2. + evalCoeff[4]*log(βEp).^3.

end


"""

calculateT12(z,Qᵦ,Eₓ,BGT)

calculate halflife of the beta decay of an isotope given feedings to excited states

z: Atomic number of the parent 

Qᵦ: β decay Q value in MeV

Eₓ: vector of daughter states relative to the ground state energy in MeV

BGT: vector of BGT values

"""
function calculateT12(z,Qᵦ,Eₓ::Vector,BGT::Vector)


    coeff = [ -17.2       7.9015    -2.54        0.28482;
           3.31368   -2.06273    0.703822   -0.075039;
          -0.364018   0.387961  -0.142528    0.016;
           0.0278071 -0.026519   0.0098854  -0.00113772
       ]  ;
    zDaughter = z + 1
    evalCoeff = [
    coeff[1,1] + log(zDaughter) * coeff[1,2] + coeff[1,3]*log(zDaughter)^2. + coeff[1,4]*log(zDaughter)^3.,
    coeff[2,1] + log(zDaughter) * coeff[2,2] + coeff[2,3]*log(zDaughter)^2. + coeff[2,4]*log(zDaughter)^3.,
    coeff[3,1] + log(zDaughter) * coeff[3,2] + coeff[3,3]*log(zDaughter)^2. + coeff[3,4]*log(zDaughter)^3.,
    coeff[4,1] + log(zDaughter) * coeff[4,2] + coeff[4,3]*log(zDaughter)^2. + coeff[4,4]*log(zDaughter)^3.  
    ]
    βEp = (Qᵦ .-  (Eₓ)) .* 1000 #convert to keV
    lf = evalCoeff[1] .+ evalCoeff[2].*log.(βEp[findall(βEp.>0)]) .+ evalCoeff[3].*log.(βEp[findall(βEp.>0)]).^2. .+ evalCoeff[4].*log.(βEp[findall(x->x>0,βEp)]).^3.

    D=6144.5
    λ=log(2) .* 10 .^lf .* BGT[findall(βEp.>0)] ./ D

    return log(2) ./ sum(λ)
    
end


"""

calculateT12fullFermi(A,z,Qᵦ,Eₓ,BGT)

Calculate the beta-decay half-life using a relativistic Fermi-function
integral for each energetically allowed transition.

`A` is the mass number used for the finite-size nuclear radius. It may be a
scalar or a vector matching `Eₓ`.

This implements allowed beta-minus decay, with `z+1` as the daughter charge.
Q values and excitation energies are in MeV, and `BGT` uses the same
normalization as `calculateT12`.
"""
function calculateT12fullFermi(A,z,Qᵦ,Eₓ::Vector,BGT::Vector)

    α = 1 / 137.035999084
    mₑ = 0.51099895 # MeV
    λₑ = 386.15926764 # electron Compton wavelength in fm
    zDaughter = z + 1
    γ = sqrt(1 - (α*zDaughter)^2)
    βEp = (Qᵦ .- Eₓ)
    allowed = findall(βEp .> 0)

    f = similar(βEp[allowed], Float64)
    for (j, i) in enumerate(allowed)
        radius = 1.2 * (A isa AbstractVector ? A[i] : A)^(1/3) / λₑ
        w₀ = 1 + βEp[i] / mₑ
        f[j] = QuadGK.quadgk(w -> begin
            p = sqrt(w^2 - 1)
            y = α * zDaughter * w / p
            F = 4 * (2*p*radius)^(2*(γ - 1)) *
                abs(SpecialFunctions.gamma(γ + im*y))^2 * exp(pi*y) /
                SpecialFunctions.gamma(2*γ + 1)^2
            p * w * (w₀ - w)^2 * F
        end, 1, w₀)[1]
    end

    D = 6144.5
    λ = log(2) .* f .* BGT[allowed] ./ D
    return log(2) / sum(λ)

end


"""
calculateIb(z,Qᵦ,Eₓ,BGT)

calculate branching ratios of the beta decay of an isotope given feedings to excited states

z: Atomic number of the parent

Qᵦ: β decay Q value in MeV

Eₓ: vector of daughter states relative to the ground state energy in MeV

BGT: vector of BGT values


"""
function calculateIb(z,Qᵦ,Eₓ::Vector,BGT::Vector)


    coeff = [ -17.2       7.9015    -2.54        0.28482;
           3.31368   -2.06273    0.703822   -0.075039;
          -0.364018   0.387961  -0.142528    0.016;
           0.0278071 -0.026519   0.0098854  -0.00113772
       ]  ;
    zDaughter = z + 1
    evalCoeff = [
    coeff[1,1] + log(zDaughter) * coeff[1,2] + coeff[1,3]*log(zDaughter)^2. + coeff[1,4]*log(zDaughter)^3.,
    coeff[2,1] + log(zDaughter) * coeff[2,2] + coeff[2,3]*log(zDaughter)^2. + coeff[2,4]*log(zDaughter)^3.,
    coeff[3,1] + log(zDaughter) * coeff[3,2] + coeff[3,3]*log(zDaughter)^2. + coeff[3,4]*log(zDaughter)^3.,
    coeff[4,1] + log(zDaughter) * coeff[4,2] + coeff[4,3]*log(zDaughter)^2. + coeff[4,4]*log(zDaughter)^3.  
    ]
    βEp = (Qᵦ .-  (Eₓ)) .* 1000 #convert to keV
    lf = evalCoeff[1] .+ evalCoeff[2].*log.(βEp[findall(βEp.>0)]) .+ evalCoeff[3].*log.(βEp[findall(βEp.>0)]).^2. .+ evalCoeff[4].*log.(βEp[findall(x->x>0,βEp)]).^3.

    D=6144.5
    λ=log(2) .* 10 .^lf .* BGT[findall(βEp.>0)] ./ D

    t₁₂=calculateT12(z,Qᵦ,Eₓ,BGT)

    return λ.*t₁₂/log(2)
    
end

"""
calculateBGT(z,Qᵦ,T₁₂,Eₓ,Iᵦ)

calculate beta decay strength of an isotope given partial branching ratios to excited states

z: Atomic number of the parent

Qᵦ: β decay Q value in MeV

T₁₂: beta decay half-life

Eₓ: vector of daughter states relative to the ground state energy in MeV

Iᵦ: vector of Iᵦ values


"""
function calculateBGT(z,Qᵦ,T₁₂,Eₓ::Vector,Iᵦ::Vector)


    coeff = [ -17.2       7.9015    -2.54        0.28482;
           3.31368   -2.06273    0.703822   -0.075039;
          -0.364018   0.387961  -0.142528    0.016;
           0.0278071 -0.026519   0.0098854  -0.00113772
       ]  ;
    zDaughter = z + 1
    evalCoeff = [
    coeff[1,1] + log(zDaughter) * coeff[1,2] + coeff[1,3]*log(zDaughter)^2. + coeff[1,4]*log(zDaughter)^3.,
    coeff[2,1] + log(zDaughter) * coeff[2,2] + coeff[2,3]*log(zDaughter)^2. + coeff[2,4]*log(zDaughter)^3.,
    coeff[3,1] + log(zDaughter) * coeff[3,2] + coeff[3,3]*log(zDaughter)^2. + coeff[3,4]*log(zDaughter)^3.,
    coeff[4,1] + log(zDaughter) * coeff[4,2] + coeff[4,3]*log(zDaughter)^2. + coeff[4,4]*log(zDaughter)^3.  
    ]
    βEp = (Qᵦ .-  (Eₓ)) .* 1000 #convert to keV
    lf = evalCoeff[1] .+ evalCoeff[2].*log.(βEp[findall(βEp.>0)]) .+ evalCoeff[3].*log.(βEp[findall(βEp.>0)]).^2. .+ evalCoeff[4].*log.(βEp[findall(x->x>0,βEp)]).^3.

    D=6144.5

    return D.*Iᵦ./(10 .^lf*T₁₂)
    
end


""" 

logftfromib(z,t₁₂,Qᵦ,Eₓ,Iᵦ)

calculate logft of a given transition to an excitated state

z: atomic number of the parent

t₁₂: decay halflife

Qᵦ: β decay Q value in MeV

Eₓ: energy of daughter state relative to the ground state energy in MeV

Iᵦ: partial branching values



"""
function logftfromib(z,t₁₂,Qᵦ,Eₓ,Iᵦ)  
    coeff = [ -17.2       7.9015    -2.54        0.28482;
    3.31368   -2.06273    0.703822   -0.075039;
   -0.364018   0.387961  -0.142528    0.016;
    0.0278071 -0.026519   0.0098854  -0.00113772
]  ;
zDaughter = z + 1
evalCoeff = [
coeff[1,1] + log(zDaughter) * coeff[1,2] + coeff[1,3]*log(zDaughter)^2. + coeff[1,4]*log(zDaughter)^3.,
coeff[2,1] + log(zDaughter) * coeff[2,2] + coeff[2,3]*log(zDaughter)^2. + coeff[2,4]*log(zDaughter)^3.,
coeff[3,1] + log(zDaughter) * coeff[3,2] + coeff[3,3]*log(zDaughter)^2. + coeff[3,4]*log(zDaughter)^3.,
coeff[4,1] + log(zDaughter) * coeff[4,2] + coeff[4,3]*log(zDaughter)^2. + coeff[4,4]*log(zDaughter)^3.  
]
βEp = (Qᵦ -  (Eₓ)) * 1000 #convert to keV
lf = evalCoeff[1] + evalCoeff[2]*log(βEp) + evalCoeff[3]*log(βEp)^2. + evalCoeff[4]*log(βEp).^3.

return log10(10^lf*t₁₂/Iᵦ)

end


"""
logftfrombgt(bgt)

calculate the logft for a given BGT (not quenched)

    
"""
function logftfrombgt(bgt)
    
    return log10((6144.5/bgt))

end


_default_neutron_radius(AM1, AM2) = 1.4 * (AM1 + AM2)^(1 / 3)

function _validate_neutron_inputs(E, AM1, AM2, Lorb, R)
    E >= 0 || throw(DomainError(E, "neutron energy must be nonnegative"))
    AM1 > 0 || throw(DomainError(AM1, "recoil mass must be positive"))
    AM2 > 0 || throw(DomainError(AM2, "neutron mass must be positive"))
    R > 0 || throw(DomainError(R, "channel radius must be positive"))
    isinteger(Lorb) && Lorb >= 0 ||
        throw(ArgumentError("Lorb must be a nonnegative integer"))
    Lorb <= 7 || throw(ArgumentError("Lorb must be ≤ 7"))
    return Int(Lorb)
end

"""
    nPenetrability(E, AM1, AM2, Lorb; R=1.4*(AM1+AM2)^(1/3))

Calculate the dimensionless neutron penetrability ``P_ℓ(E,R)`` from Blatt
and Weisskopf and Lane and Thomas. `E` is the neutron center-of-mass energy in
MeV above threshold, and `R` is the channel radius in fm.

`AM1` and `AM2` are the recoil and neutron masses in atomic mass units, and
`Lorb` is the neutron orbital angular momentum.

"""
function nPenetrability(x, AM1, AM2, Lorb; R=_default_neutron_radius(AM1, AM2))
    ℓmax = _validate_neutron_inputs(x, AM1, AM2, Lorb, R)
    E = x * 1000.0       # Energy in keV
    RMAS = AM1*AM2/(AM1+AM2)*931502
    k = sqrt(2*RMAS*E)/197329
    
    σ = k*R

    s = zeros(8)
    p = zeros(8)
    p[1] = σ

    for ℓ in 1:7
        denominator = (ℓ - s[ℓ])^2 + p[ℓ]^2
        s[ℓ + 1] = σ^2 * (ℓ - s[ℓ]) / denominator - ℓ
        p[ℓ + 1] = p[ℓ] * σ^2 / denominator
    end

    return p[ℓmax + 1]

end


"""
calculateBE2T12(A,Eᵧ,BE2,units=nothing)

calculate the halflife of a E2 nuclear electromagnetic transition (does not include correction for conversion electrons)

A: nuclear mass number
Eᵧ: transition energy in MeV
BE2: transition probability in e2fm4 (for BE2 in w.u. define units="wu")


"""
function calculateBE2T12(A,Eᵧ,BE2,units=nothing)
    if units=="wu" BE2e2fm4 = BE2 * 0.0594 * A^(4/3);
    else BE2e2fm4 = BE2
    end
    λ = 1.2239e9 * Eᵧ^5 * BE2e2fm4
    return log(2)/λ
end

"""
    Gamma(E, AM1, AM2, Lorb; R=1.4*(AM1+AM2)^(1/3), warn_validity=false)

Calculate the Bohr--Mottelson low-energy square-well estimate of the neutron
single-particle width,

``Γ_ℓ(E) = 2P_ℓ(E,R) ((ħc)^2/(μc^2 R^2)) (2ℓ-1)/(2ℓ+1)``,

in MeV. `E` is the neutron center-of-mass energy above threshold in MeV,
`AM1` and `AM2` are masses in atomic mass units, `Lorb` is the neutron orbital
angular momentum, and `R` is the channel radius in fm.

This approximation requires `Lorb ≥ 1` and is derived for ``kR ≪ 1``. Set
`warn_validity=true` to emit a warning when ``kR ≥ 1``. The result is a
depth-independent single-particle estimate, not an exact S-matrix pole width;
a physical partial width may also require a spectroscopic factor.

The former signatures `Gamma(E, V0, AM1, AM2, Lorb; R=...)` and
`Gamma(E, V0, AM1, AM2, Lorb, R)` remain as deprecated compatibility wrappers;
`V0` is ignored.
"""
function Gamma(E, AM1, AM2, Lorb;
               R=_default_neutron_radius(AM1, AM2), warn_validity=false)
    ħc = 197.326
    ℓ = _validate_neutron_inputs(E, AM1, AM2, Lorb, R)
    ℓ >= 1 || throw(ArgumentError(
        "the Bohr--Mottelson square-well width approximation requires Lorb ≥ 1"))
    RMAS = AM1*AM2/(AM1+AM2)*931.502
    kR = R * sqrt(2*RMAS*E) / ħc
    if warn_validity && kR >= 1
        @warn "Bohr--Mottelson square-well width is a low-energy approximation" kR
    end
    P = nPenetrability(E, AM1, AM2, ℓ; R=R)
    γ² = ħc^2/(RMAS*R^2) * (2ℓ - 1)/(2ℓ + 1)
    return 2*P*γ²
end

function Gamma(E, V0, AM1, AM2, Lorb; R=_default_neutron_radius(AM1, AM2),
               warn_validity=false)
    Base.depwarn("Gamma(E, V0, AM1, AM2, Lorb) is deprecated because V0 is " *
                 "unused; call Gamma(E, AM1, AM2, Lorb; R=R) instead", :Gamma)
    return Gamma(E, AM1, AM2, Lorb; R=R, warn_validity=warn_validity)
end

function Gamma(E, V0, AM1, AM2, Lorb, R; warn_validity=false)
    Base.depwarn("Gamma(E, V0, AM1, AM2, Lorb, R) is deprecated because V0 is " *
                 "unused; call Gamma(E, AM1, AM2, Lorb; R=R) instead", :Gamma)
    return Gamma(E, AM1, AM2, Lorb; R=R, warn_validity=warn_validity)
end

end
