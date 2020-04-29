function getCompoundElementProperties(comp::Symbol, compoundData::compoundData{:HexToTris})
    if comp==:RT0
        refBound=Dict([1,2]=>[1,0,0,0,0,0], [2,3]=>[0,1,0,0,0,0], [3,4]=>[0,0,1,0,0,0],
                      [4,5]=>[0,0,0,1,0,0], [5,6]=>[0,0,0,0,1,0], [6,1]=>[0,0,0,0,0,1]);
        edgeTypes=Dict([1,2]=>1, [2,3]=>2, [3,4]=>3, [4,5]=>4, [5,6]=>5, [6,1]=>6);
    elseif comp==:RT1
        @error("refBound not correct")
        refBound=Dict([1,2]=>[1,0,0,0,0,0], [2,3]=>[0,1,0,0,0,0], [3,4]=>[0,0,1,0,0,0],
                      [4,5]=>[0,0,0,1,0,0], [5,6]=>[0,0,0,0,1,0], [6,1]=>[0,0,0,0,0,1]);
        edgeTypes=Dict([1,2]=>1, [2,3]=>2, [3,4]=>3, [4,5]=>4, [5,6]=>5, [6,1]=>6);
    end

    return refBound, edgeTypes
end


function getCompoundElementProperties(comp::Symbol, compoundData::compoundData{:HexToKites})
    if comp==:RT0
        refBound=Dict([1,2]=>[1,0,0,0,0,0], [2,3]=>[0,1,0,0,0,0], [3,4]=>[0,0,1,0,0,0],
                      [4,5]=>[0,0,0,1,0,0], [5,6]=>[0,0,0,0,1,0], [6,1]=>[0,0,0,0,0,1]);
        edgeTypes=Dict([1,2]=>1, [2,3]=>2, [3,4]=>3, [4,5]=>4, [5,6]=>5, [6,1]=>6);
    elseif comp==:RT1
        @error("refBound not correct")
        refBound=Dict([1,2]=>[1,0,0,0,0,0], [2,3]=>[0,1,0,0,0,0], [3,4]=>[0,0,1,0,0,0],
                      [4,5]=>[0,0,0,1,0,0], [5,6]=>[0,0,0,0,1,0], [6,1]=>[0,0,0,0,0,1]);
        edgeTypes=Dict([1,2]=>1, [2,3]=>2, [3,4]=>3, [4,5]=>4, [5,6]=>5, [6,1]=>6);
    end

    return refBound, edgeTypes
end
