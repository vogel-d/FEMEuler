function transp(m::Array{Array{T,N} where N,2} ) where T
    res=Array{Array{T,N} where N,2}(undef,size(m,2),size(m,1));
    for j in 1:size(m,2)
        for i in 1:size(m,1)
            res[j,i]=m[i,j]
        end
    end
    return res
end
function transp(m::Array{Array{T,1},2} ) where T
    res=Array{Array{T,1},2}(undef,size(m,2),size(m,1));
    for j in 1:size(m,2)
        for i in 1:size(m,1)
            res[j,i]=m[i,j]
        end
    end
    return res
end
function transp(m::Array{Array{T,2},2} ) where T
    res=Array{Array{T,2},2}(undef,size(m,2),size(m,1));
    for j in 1:size(m,2)
        for i in 1:size(m,1)
            res[j,i]=m[i,j]
        end
    end
    return res
end
