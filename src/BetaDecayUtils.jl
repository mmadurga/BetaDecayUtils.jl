module BetaDecayUtils

export  calculateqbetashellmodel,
        logf,
        logftfromib,
        calculateT12,
        childActivity, 
        grandChildActivity,
        chainActivity,
        nPenetrability,
        logftfrombgt,
        calculateIb,
        calculateBGT,
        ecoulomb,
        gammafit,
        wavefunction,
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

    D=6147
    λ=log(2) .* 10 .^lf .* BGT[findall(βEp.>0)] ./ D

    return log(2) ./ sum(λ)
    
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

    D=6147
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

    D=6147

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
    
    return log10((6147/bgt))

end


"""
nPenetrability(x,AM1,AM2,Lorb)

calculates the neutron penetrability p(x,Lorb) from Blatt and Weisskopf "Theoretical Nuclear Physics" p. 463

x is the excitation energy above Sₙ, Lorb is the neutron angular momentum

AM1 is the recoil, AM2 is the neutron mass

"""
function nPenetrability(x,AM1,AM2,Lorb)
    
    E = x  * 1000.       #Energy in keV
    R = 1.4  * (AM1 + AM2)^(1/3)
    RMAS = AM1*AM2/(AM1+AM2)*931502
    eo = (197329^2)/(2*RMAS*R^2)
    k = sqrt(2*RMAS*E)/197329
    
    σ = k*R

    s0 = 0
    p0 = σ
    
    s1 = -1/(1+σ^2)
    p1 = σ^3/(1+σ^2)
    
    s2 = (σ^2*(2-s1) / ((2-s1)^2+p1^2) ) - 2 
    p2 = p1*(σ^2/((2-s1)^2+p1^2))

    s3 = (σ^2*(3-s2) / ((3-s2)^2+p2^2) ) - 3 
    p3 = p2*(σ^2/((3-s2)^2+p2^2))

    s4 = (σ^2*(4-s3) / ((4-s3)^2+p3^2) ) - 4
    p4 = p3*(σ^2/((4-s3)^2+p3^2))

    s5 = (σ^2*(5-s4) / ((5-s4)^2+p4^2) ) - 5 
    p5 = p4*(σ^2/((5-s4)^2+p4^2))

    s6 = (σ^2*(6-s5) / ((6-s5)^2+p5^2) ) - 6
    p6 = p5*(σ^2/((6-s5)^2+p5^2))

    s7 = (σ^2*(7-s6) / ((7-s6)^2+p6^2) ) - 7
    p7 = p6*(σ^2/((6-s6)^2+p6^2))

    if (Lorb==0.) 
        return p0
    elseif (Lorb==1.) 
        return p1
    elseif (Lorb==2.) 
        return p2
    elseif (Lorb==3.) 
        return p3
    elseif (Lorb==4.) 
        return p4
    elseif (Lorb==5.)
        return p5
    elseif (Lorb==6.)
        return p6
    elseif (Lorb==7.)
        return p7
    end

end



"""
wavefunction(E,V0,AM1,AM2,Lorb,r)

calculates the inner and outer wavefunctions as a function of the spherical radial coordinate r for neutron emission using a spherical square well potential

E is the excitation energy above Sₙ, V0 is the potential depth, AM1 and AM2 are the residue's and neutron mass in au and Lorb is the neutron angular momentum

"""
function wavefunction(E,V0,AM1,AM2,Lorb,r)
 
    ħc = 197.326

    RMAS = AM1*AM2/(AM1+AM2)*931.502
    Ro = 1.4*(AM1+AM2)^0.3333

    k = sqrt(2 * RMAS * (E + V0)) / ħc
    q = sqrt(2 * RMAS * E) / ħc

    A = sqrt(1/(QuadGK.quadgk(r->SpecialFunctions.sphericalbesselj(Lorb,k*r)^2*r^2,0,Ro)[1]))
    C = A*SpecialFunctions.sphericalbesselj(Lorb,k*Ro)/(SpecialFunctions.sphericalbesselj(Lorb,q*Ro)+SpecialFunctions.sphericalbessely(Lorb,q*Ro))
    
    u1 = A*SpecialFunctions.sphericalbesselj(Lorb,k*r)*r
    o1 = C*(SpecialFunctions.sphericalbesselj(Lorb,q*r)+SpecialFunctions.sphericalbessely(Lorb,q*r))*r
    
    return u1,o1
end

"""
Gamma(E,V0,AM1,AM2,Lorb,r)

calculates the neutron emission single particle width assuming a square well potential.
E is the excitation energy above Sₙ, V0 is the potential depth, AM1 and AM2 are the residue's and neutron mass in au and Lorb is the neutron angular momentum

"""
function Gamma(E,V0,AM1,AM2,Lorb) 
   
    ħc = 197.326
    R = 1.48 * (AM1 + AM2)^(1/3)
    RMAS = AM1*AM2/(AM1+AM2)*931.502

    return 2*ħc^2/(2*RMAS*R)*wavefunction(E,V0,AM1,AM2,Lorb,R)[1]^2*nPenetrability(E,AM1,AM2,Lorb)

end

end