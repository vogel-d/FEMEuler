function getTangentialPlane(n::Array{Float64,1})
    length(n)==2 && push!(n,1.0);
    #length(n)==2 && return [],[];
    ind=partialsortperm(abs.(n),1:2,rev=true)
    t1=zeros(3);
    t1[ind[1]]=n[ind[2]]
    t1[ind[2]]=-n[ind[1]]
    t2=cross(n,t1);
    return t1, t2;
end
