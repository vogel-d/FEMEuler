function getOrderBoundary(boundary::SparseVector{Int,Int})
    sortBoundary=collect(1:length(boundary));
    bound=Int[];
    int=Int[];
    d=0;
    anzP=0;
    anzC=0;
    for i in 1:length(boundary)
        if boundary[i]<0
            sortBoundary[i]=0;
            d+=1;
            anzP+=1;
        else
            sortBoundary[i]-=d;
            if boundary[i]==1
                push!(bound,i);
                anzC+=1;
            else
                push!(int,i);
            end
        end
    end
    ordB=sortBoundary[int];
    append!(ordB, sortBoundary[bound]);
    invert!(ordB);
    sortBoundary[int]=ordB[sortBoundary[int]];
    sortBoundary[bound]=ordB[sortBoundary[bound]];
    return sortBoundary, anzP, anzC;
end
