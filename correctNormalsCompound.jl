function correctNormalsCompound!(n1::Array{Float64,1},n2::Array{Float64,1},eT1::Int64,eT2::Int64,compoundData::compoundData{:RectToKites})
    if in(eT1,[3,4])
        n1[1]=-n1[1];
        n1[2]=-n1[2];
    else
        #n1=n1;
    end
    if in(eT2,[3,4])
        n2[1]=-n2[1];
        n2[2]=-n2[2];
    else
        #n2=n2;
    end
end


function correctNormalsCompound!(n1::Array{Float64,1},n2::Array{Float64,1},eT1::Int64,eT2::Int64,compoundData::compoundData)
    @warn("wrong normalcorrection")
    if in(eT1,[4,5,6])
        n1[1]=-n1[1];
        n1[2]=-n1[2];
    else
        #n1=n1;
    end
    if in(eT2,[4,5,6])
        n2[1]=-n2[1];
        n2[2]=-n2[2];
    else
        #n2=n2;
    end
end
