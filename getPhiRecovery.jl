macro recovery(recoverySpace, recoverySpaceVec)
    fname=:getPhiRecovery!
    quote
        rS = $(esc(recoverySpace))
        rSV = $(esc(recoverySpaceVec))
        if rS==:DG1
            function $(esc(fname))(phi::Array{Float64,1},xyz::Array{Float64,1})

                #phi[:]=[1.0, xyz[1], xyz[2]]
                phi[:]=[(1-xyz[1])*(1-xyz[2]), xyz[1]*(1-xyz[2]), xyz[1]*xyz[2], (1-xyz[1])*xyz[2]]

                return nothing
            end
        elseif rS==:DGLin
            function $(esc(fname))(phi::Array{Float64,1},xyz::Array{Float64,1})

                phi[:]=[1.0, xyz[1]-0.5, xyz[2]-0.5]

                return nothing
            end
        elseif rS==:DGQuad
            function $(esc(fname))(phi::Array{Float64,1},xyz::Array{Float64,1})

                phi[:]=[1.0, xyz[1]-0.5, xyz[2]-0.5, (xyz[1]-0.5)*(xyz[2]-0.5), (xyz[1]-0.5)^2, (xyz[2]-0.5)^2]

                return nothing
            end
        elseif rS==:R2
            function $(esc(fname))(phi::Array{Float64,1},xyz::Array{Float64,1})

                phi[:]=[1.0, xyz[1], xyz[2], xyz[1]*xyz[2], xyz[1]^2, xyz[2]^2]

                return nothing
            end
        end

        if rSV==:VecDG1
            function $(esc(fname))(phi::Array{Float64,2},xyz::Array{Float64,1})

                phi[:]=[(1-xyz[1])*(1-xyz[2]) 0.0    0.0    xyz[1]*(1-xyz[2])  0.0    0.0    xyz[1]*xyz[2]  0.0    0.0    (1-xyz[1])*xyz[2]  0.0    0.0;
                     0.0   (1-xyz[1])*(1-xyz[2])  0.0    0.0    xyz[1]*(1-xyz[2])  0.0    0.0    xyz[1]*xyz[2]  0.0    0.0    (1-xyz[1])*xyz[2]  0.0;
                     0.0   0.0    (1-xyz[1])*(1-xyz[2])  0.0    0.0    xyz[1]*(1-xyz[2])  0.0    0.0    xyz[1]*xyz[2]  0.0    0.0    (1-xyz[1])*xyz[2]]

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
    if recoverySpace==:DG1
        #nW=3;
        nW=4
    elseif recoverySpace==:DGLin
        nW=3;
    elseif recoverySpace==:DGQuad
        nW=6;
    elseif recoverySpace==:R2
        nW=6;
    elseif recoverySpace==:VecDG1
        #nW=6;
        nW=12
    elseif recoverySpace==:VecR1S
        nW=9;
    elseif recoverySpace==:VecR2
        nW=12;
    elseif recoverySpace==:VecR2S
        nW=18;
    end
    return nW
 end
