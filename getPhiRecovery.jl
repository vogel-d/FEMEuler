macro recovery(recoverySpace, recoverySpaceVec)
    fname=:getPhiRecovery!
    quote
        rS = $(esc(recoverySpace))
        rSV = $(esc(recoverySpaceVec))
        if rS==:R1
            function $(esc(fname))(phi::Array{Float64,1},xyz::Array{Float64,1})

                phi[:]=[1.0, xyz[1], xyz[2]]

                return nothing
            end
        elseif rS==:R2
            function $(esc(fname))(phi::Array{Float64,1},xyz::Array{Float64,1})

                phi[:]=[1.0, xyz[1], xyz[2], xyz[1]*xyz[2], xyz[1]^2, xyz[2]^2]

                return nothing
            end
        end

        if rSV==:VecR1
            function $(esc(fname))(phi::Array{Float64,2},xyz::Array{Float64,1})

                phi[:]=[1.0 xyz[1] xyz[2] 0.0 0.0 0.0;
                     0.0 0.0 0.0 1.0 xyz[1] xyz[2]]

                return nothing
            end
        elseif rSV==:VecR1S
            function $(esc(fname))(phi::Array{Float64,2},xyz::Array{Float64,1})

                phi[:]=[1.0 xyz[1] xyz[2] 0.0 0.0 0.0 0.0 0.0 0.0;
                     0.0 0.0 0.0 1.0 xyz[1] xyz[2] 0.0 0.0 0.0;
                     0.0 0.0 0.0 0.0 0.0 0.0 1.0 xyz[1] xyz[2]]

                return nothing
            end
        elseif rSV==:VecR2
            function $(esc(fname))(phi::Array{Float64,2},xyz::Array{Float64,1})

                phi[:]=[1.0 xyz[1] xyz[2] xyz[1]*xyz[2] xyz[1]^2 xyz[2]^2 0.0 0.0 0.0 0.0 0.0 0.0;
                     0.0 0.0 0.0 0.0 0.0 0.0 1.0 xyz[1] xyz[2] xyz[1]*xyz[2] xyz[1]^2 xyz[2]^2]

                return nothing
            end
        elseif rSV==:VecR2S
            function $(esc(fname))(phi::Array{Float64,2},xyz::Array{Float64,1})

                phi[:]=[1.0 xyz[1] xyz[2] xyz[1]*xyz[2] xyz[1]^2 xyz[2]^2 0.0 0.0 0.0 0.0  0.0  0.0  0.0 0.0 0.0 0.0  0.0  0.0;
                     0.0 0.0 0.0 0.0  0.0  0.0  1.0 xyz[1] xyz[2] xyz[1]*xyz[2] xyz[1]^2 xyz[2]^2 0.0 0.0 0.0 0.0  0.0  0.0;
                     0.0 0.0 0.0 0.0  0.0  0.0  0.0 0.0 0.0 0.0  0.0  0.0  1.0 xyz[1] xyz[2] xyz[1]*xyz[2] xyz[1]^2 xyz[2]^2]

                return nothing
            end
        end
    end
end



function getPhiRecoveryLength(recoverySpace::Symbol)
    if recoverySpace==:R1
        nW=3;
    elseif recoverySpace==:R2
        nW=6;
    elseif recoverySpace==:VecR1
        nW=6;
    elseif recoverySpace==:VecR1S
        nW=9;
    elseif recoverySpace==:VecR2
        nW=12;
    elseif recoverySpace==:VecR2S
        nW=18;
    end
    return nW
 end
