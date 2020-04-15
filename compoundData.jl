struct compoundData
    nSubCells::Int64;
    getSubCells::Symbol;
    assembledPhi::Dict{Symbol,Array{Array{Float64,2},1}};
end

function createCompoundData()
    compoundData(0,:zero,Dict());
end

function createCompoundData(method::Symbol,femElements::Set{Symbol})
    if method==:HexToKites
        return compoundData(6,:getSubCellsHexToKites!,assemblePhiHexToKites(femElements))
    end
end

function assemblePhiHexToKites(femElements::Set{Symbol})
    assembledPhi=Dict();
    for type in femElements
        if type==:RT0
            phi2times = (1/3)*
                        [0  3  0 0 0  0  0  3  0 0 0  0;
                         3  0  0 0 0  0  3  0  0 0 0  0;
                         0 -2 -1 0 1  2  0 -2 -1 0 1  2;
                        -2  0  2 1 0 -1 -2  0  2 1 0 -1];

            phi1 = phi2times[:,7:12];

            phi2 = phi2times[:,6:11];

            phi3 = phi2times[:,5:10];

            phi4 = phi2times[:,4:9];

            phi5 = phi2times[:,3:8];

            phi6 = phi2times[:,2:7];

            assembledPhi[type]=[phi1, phi2, phi3, phi4, phi5, phi6];
        elseif type==:DG0
            phi2times = ones(Float64,1,12);

            phi1 = phi2times[:,7:12];

            phi2 = phi2times[:,6:11];

            phi3 = phi2times[:,5:10];

            phi4 = phi2times[:,4:9];

            phi5 = phi2times[:,3:8];

            phi6 = phi2times[:,2:7];

            assembledPhi[type]=[phi1, phi2, phi3, phi4, phi5, phi6];
        else
            error("Tell me how to assemble compound ansatzfunctions in $type")
        end
    end
    return assembledPhi;
end


function testest()
    lal=getfield(Main,:println)
    lal("was geht")
end
