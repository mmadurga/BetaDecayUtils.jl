"""Archived explicit-recurrence implementation of `nPenetrability`."""
function nPenetrability(x,AM1,AM2,Lorb)
    if Lorb > 7
        throw(ArgumentError("Lorb must be ≤ 7"))
    end
    E = x * 1000.
    R = 1.4 * (AM1 + AM2)^(1/3)
    RMAS = AM1*AM2/(AM1+AM2)*931502
    eo = (197329^2)/(2*RMAS*R^2)
    k = sqrt(2*RMAS*E)/197329
    σ = k*R
    
    s0 = 0
    p0 = σ
 
    s1 = -1/(1+σ^2)
    p1 = σ^3/(1+σ^2)
 
    s2 = (σ^2*(2-s1) / ((2-s1)^2+p1^2)) - 2
    p2 = p1*(σ^2/((2-s1)^2+p1^2))
 
    s3 = (σ^2*(3-s2) / ((3-s2)^2+p2^2)) - 3
    p3 = p2*(σ^2/((3-s2)^2+p2^2))
 
    s4 = (σ^2*(4-s3) / ((4-s3)^2+p3^2)) - 4
    p4 = p3*(σ^2/((4-s3)^2+p3^2))
 
    s5 = (σ^2*(5-s4) / ((5-s4)^2+p4^2)) - 5
    p5 = p4*(σ^2/((5-s4)^2+p4^2))
 
    s6 = (σ^2*(6-s5) / ((6-s5)^2+p5^2)) - 6
    p6 = p5*(σ^2/((6-s5)^2+p5^2))
 
    s7 = (σ^2*(7-s6) / ((7-s6)^2+p6^2)) - 7
    p7 = p6*(σ^2/((7-s6)^2+p6^2))
 
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
