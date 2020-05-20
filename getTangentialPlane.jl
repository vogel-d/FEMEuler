function getTangentialPlane(n::Array{Float64,1})
    length(n)==2 && return [],[];
    normalize!(n)
    ind=partialsortperm(abs.(n),1:2,rev=true)
    t1=zeros(3);
    t1[ind[1]]=n[ind[2]]
    t1[ind[2]]=-n[ind[1]]
    normalize!(t1)
    t2=cross(n,t1);
    normalize!(t2)
    return t1, t2;
end
#=
function getTangentialPlane(n::Array{Float64,1})
    ind=partialsortperm(abs.(n),1:2,rev=true)
    t1=zeros(3);
    t1[ind[1]]=n[ind[2]]
    t1[ind[2]]=-n[ind[1]]
    if length(n)==2
        t1=[1.0,0.0,0.0]
        t2=[0.0,1.0,0.0]
        push!(n,1.0)
    else
        t2=cross(n,t1);
    end
    normalize!(t1)
    normalize!(t2)
    return t1, t2;
end
=#
