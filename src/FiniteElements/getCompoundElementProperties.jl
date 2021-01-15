function getCompoundElementProperties(type::Symbol, mt::Int, x, y, compoundData::compoundData)
    if mt==4
        phi, divphi, gradphi, cm, nFace, nEdge, nVert=getQuadElementProperties(type);
    elseif mt==3
        phi, divphi, gradphi, cm, nFace, nEdge, nVert=getTriElementProperties(type);
    end
    valPhi=similar(phi,Float64);
    for k=1:length(phi)
        valPhi[k]=phi[k](x,y);
    end
    return valPhi
end



function getCompoundElementProperties(comp::Symbol, compoundData::compoundData{:HexToTris})
    if comp==:RT0
        refBound=Dict([1,2]=>[1,0,0,0,0,0], [2,3]=>[0,1,0,0,0,0], [3,4]=>[0,0,1,0,0,0],
                      [4,5]=>[0,0,0,1,0,0], [5,6]=>[0,0,0,0,1,0], [1,6]=>[0,0,0,0,0,1]);
        edgeTypes=Dict([1,2]=>1, [2,3]=>2, [3,4]=>3, [4,5]=>4, [5,6]=>5, [1,6]=>6);
    elseif comp==:RT1
        @error("refBound not correct")
        refBound=Dict([1,2]=>[1,0,0,0,0,0], [2,3]=>[0,1,0,0,0,0], [3,4]=>[0,0,1,0,0,0],
                      [4,5]=>[0,0,0,1,0,0], [5,6]=>[0,0,0,0,1,0], [1,6]=>[0,0,0,0,0,1]);
        edgeTypes=Dict([1,2]=>1, [2,3]=>2, [3,4]=>3, [4,5]=>4, [5,6]=>5, [1,6]=>6);
    end

    return refBound, edgeTypes
end


function getCompoundElementProperties(comp::Symbol, compoundData::compoundData{:HexToKites})
    if comp==:RT0
        refBound=Dict([1,2]=>[1,0,0,0,0,0], [2,3]=>[0,1,0,0,0,0], [3,4]=>[0,0,1,0,0,0],
                      [4,5]=>[0,0,0,1,0,0], [5,6]=>[0,0,0,0,1,0], [1,6]=>[0,0,0,0,0,1]);
        edgeTypes=Dict([1,2]=>1, [2,3]=>2, [3,4]=>3, [4,5]=>4, [5,6]=>5, [1,6]=>6);
    elseif comp==:RT1
        @error("refBound not correct")
        refBound=Dict([1,2]=>[1,0,0,0,0,0], [2,3]=>[0,1,0,0,0,0], [3,4]=>[0,0,1,0,0,0],
                      [4,5]=>[0,0,0,1,0,0], [5,6]=>[0,0,0,0,1,0], [1,6]=>[0,0,0,0,0,1]);
        edgeTypes=Dict([1,2]=>1, [2,3]=>2, [3,4]=>3, [4,5]=>4, [5,6]=>5, [1,6]=>6);
    end

    return refBound, edgeTypes
end

function getCompoundElementProperties(comp::Symbol, compoundData::compoundData{:RectToKites})
    if comp==:RT0
        refBound=Dict([1,2]=>[1,0,0,0], [2,3]=>[0,1,0,0],
                      [3,4]=>[0,0,1,0], [1,4]=>[0,0,0,1]);
        edgeTypes=Dict([1,2]=>1, [2,3]=>2, [3,4]=>3, [1,4]=>4);
    elseif comp==:RT1
        @error("refBound not correct")
        refBound=Dict([1,2]=>[1,0,0,0,0,0], [2,3]=>[0,1,0,0,0,0], [3,4]=>[0,0,1,0,0,0],
                      [4,5]=>[0,0,0,1,0,0], [5,6]=>[0,0,0,0,1,0], [1,6]=>[0,0,0,0,0,1]);
        edgeTypes=Dict([1,2]=>1, [2,3]=>2, [3,4]=>3, [4,5]=>4, [5,6]=>5, [1,6]=>6);
    end

    return refBound, edgeTypes
end

function getCompoundElementProperties(comp::Symbol, compoundData::compoundData{:RectToTris})
    if comp==:RT0
        refBound=Dict([1,2]=>[1,0,0,0], [2,3]=>[0,1,0,0],
                      [3,4]=>[0,0,1,0], [1,4]=>[0,0,0,1]);
        edgeTypes=Dict([1,2]=>1, [2,3]=>2, [3,4]=>3, [1,4]=>4);
    elseif comp==:RT1
        @error("refBound not correct")
        refBound=Dict([1,2]=>[1,0,0,0,0,0], [2,3]=>[0,1,0,0,0,0], [3,4]=>[0,0,1,0,0,0],
                      [4,5]=>[0,0,0,1,0,0], [5,6]=>[0,0,0,0,1,0], [1,6]=>[0,0,0,0,0,1]);
        edgeTypes=Dict([1,2]=>1, [2,3]=>2, [3,4]=>3, [4,5]=>4, [5,6]=>5, [1,6]=>6);
    end

    return refBound, edgeTypes
end
