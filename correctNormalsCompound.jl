function correctNormalsCompound!(n1::Array{Float64,1},n2::Array{Float64,1},eT1::Int64,eT2::Int64,compoundData::compoundData{:RectToKites})
    if eT1>=3
        n1[1]=-n1[1];
        n1[2]=-n1[2];
    else
        #n1=n1;
    end
    if eT2>=3
        n2[1]=-n2[1];
        n2[2]=-n2[2];
    else
        #n2=n2;
    end
end


function correctNormalsCompound!(n1::Array{Float64,1},n2::Array{Float64,1},eT1::Int64,eT2::Int64,compoundData::compoundData)
    if eT1>=4
        n1[1]=-n1[1];
        n1[2]=-n1[2];
    else
        #n1=n1;
    end
    if eT2>=4
        n2[1]=-n2[1];
        n2[2]=-n2[2];
    else
        #n2=n2;
    end
end
