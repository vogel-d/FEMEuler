function assemblePhiHexToKites(femElements::Set{Symbol})
    @warn("not right for quadrilaterals, check betas and order of ansatzfunctions")
    assembledPhi=Dict();
    for type in femElements
        if type==:RT0
            #phi2times includes the factors one compound ansatzfunction needs
            #from each subelement (coloumns) with its ansatzfunctions (rows)
            #assemblation of every individual phi shows cyclic behaviour
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

            assembledPhi[type]=[phi1];
        else
            error("Tell me how to assemble HexToKite compound ansatzfunctions in $type")
        end
    end
    return assembledPhi;
end


function assemblePhiHexToTris(femElements::Set{Symbol})
    assembledPhi=Dict();
    for type in femElements
        if type==:RT0
            #phi2times includes the factors one compound ansatzfunction needs
            #from each subelement (coloumns) with its ansatzfunctions (rows)
            #assemblation of every individual phi shows cyclic behaviour
            phi2times = (1/6)*
                        [6  6  0  0  0  0  0 0 0 0 0 0 6  6  0  0  0  0  0 0 0 0 0 0;
                         0 -5 -4 -3 -2 -1  0 1 2 3 4 5 0 -5 -4 -3 -2 -1  0 1 2 3 4 5;
                         5  0 -5 -4 -3 -2 -1 0 1 2 3 4 5  0 -5 -4 -3 -2 -1 0 1 2 3 4];

            phi1 = phi2times[:,5:16];
            phi2 = phi2times[:,3:14];
            phi3 = phi2times[:,13:24];
            phi4 = phi2times[:,11:22];
            phi5 = phi2times[:,9:20];
            phi6 = phi2times[:,7:18];

            assembledPhi[type]=[phi1,phi2,phi3,phi4,phi5,phi6];

        elseif type==:DG0
            phi2times = ones(Float64,1,24);

            phi1 = phi2times[:,5:16];

            assembledPhi[type]=[phi1];

        else
            error("Tell me how to assemble HexToTri compound ansatzfunctions in $type")
        end
    end
    return assembledPhi;
end
