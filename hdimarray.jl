function hdimarray(s::Array{Array{Float64,2},1})
    sh=Array{Float64,3}(undef,size(s[1],1),size(s[1],2),length(s));
    for k in 1:length(s)
        for j in 1:size(s[1],2)
            for i in 1:size(s[1],1)
                sh[i,j,k]=s[k][i,j]
            end
        end
    end
    return sh;
end

function hdimarray(s::Array{Array{Float64,2},2})
    sh=Array{Float64,4}(undef,size(s,1),size(s[1],1),size(s[1],2),size(s,2));
    for k in 1:size(s,2)
        for l in 1:size(s,1)
            for j in 1:size(s[1],2)
                for i in 1:size(s[1],1)
                    sh[l,i,j,k]=s[l,k][i,j]
                end
            end
        end
    end
    return sh;
end
