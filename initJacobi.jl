function initJacobi(size::Tuple{Int64,Int64}, sk::Tuple{Int64,Int64})
    phi=Array{Array{Float64,2},2}(undef,size);
    for k in 1:size[2]
        for j in 1:size[1]
            phi[j,k]=Array{Float64,2}(undef,sk);
            phi[j,k]=Array{Float64,2}(undef,sk);
        end
    end
    return phi;
end

function initJacobi(size::Tuple{Int64,Int64}, sk::Int64)
    phi=Array{Array{Float64,1},2}(undef,size);
    for k in 1:size[2]
        for j in 1:size[1]
            phi[j,k]=Array{Float64,1}(undef,sk);
            phi[j,k]=Array{Float64,1}(undef,sk);
        end
    end
    return phi;
end
