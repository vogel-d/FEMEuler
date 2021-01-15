import Base.findall
function findall(x::Array{T,1} where T<:Real, v::Array{T,2} where T<:Real,atol=0)
    sol=Array{Int64,1}(undef,0)
    #length(x) != size(v,1) && error("Dimensionen passen nicht zusammen.")

    (rows::Int64, cols::Int64)=size(v);
    for i in 1:cols
        for j in 1:rows
            if !isapprox(x[j],v[j,i],atol=atol)
                break;
            elseif j==rows
                push!(sol,i);
            end
        end
    end

    return sol;
end


function findall(x::Real, v::Array{T,1} where T<:Real)
    sol=Array{Int64,1}(undef,0)
    for i in eachindex(v)
        if isapprox(v[i],x)
            push!(sol,i)
        end
    end

    return sol;
end
