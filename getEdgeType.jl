function getEdgeType(mt::Int, n::Array{Float64,1})
    if mt==3
        if n==[0.0,1.0]
            edgeType=1;
        elseif n==[0.0,-1.0]
            edgeType=1;
        elseif n==[1.0,1.0]
            edgeType=2;
        elseif n==[1.0,0.0]
            edgeType=3;
        elseif n==[-1.0,0.0]
            edgeType=3;
        else
            error("Kein zulässiger Normalenvektor: $(n)");
        end
    elseif mt==4
        if n==[0.0,-1.0]
            edgeType=1;
        elseif n==[1.0,0.0]
            edgeType=2;
        elseif n==[0.0, 1.0]
            edgeType=3;
        elseif n==[-1.0,0.0]
            edgeType=4;
        else
            error("Kein zulässiger Normalenvektor: $(n)");
        end
    end
    return edgeType
end
