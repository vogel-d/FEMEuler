function initJacobi(size::Tuple{Int,Int}, sk::Tuple{Int,Int})
    phi=Array{Array{AbstractFloat,2},2}(undef,size);
    for k in 1:size[2]
        phi[1,k]=Array{AbstractFloat,2}(undef,sk);
        phi[2,k]=Array{AbstractFloat,2}(undef,sk);
    end
    return phi;
end

function initJacobi(size::Tuple{Int,Int}, sk::Int)
    phi=Array{Array{AbstractFloat,1},2}(undef,size);
    for k in 1:size[2]
        phi[1,k]=Array{AbstractFloat,1}(undef,sk);
        phi[2,k]=Array{AbstractFloat,1}(undef,sk);
    end
    return phi;
end
