function getTangentialPlane!(t1::Array{Float64,1}, t2::Array{Float64,1}, n::Array{Float64,1},ind::Array{Int,1})
    length(n)==2 && return nothing;
    normalize!(n)
    partialsortperm!(ind,abs.(n),1:2,rev=true,initialized=true)
    fill!(t1,0.0)
    t1[ind[1]]=n[ind[2]]
    t1[ind[2]]=-n[ind[1]]
    normalize!(t1)
    t2[:]=cross(n,t1);
    normalize!(t2)
    return nothing;
end
