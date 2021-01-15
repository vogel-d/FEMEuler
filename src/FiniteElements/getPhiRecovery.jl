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
        elseif rS==:DGTri
            function $(esc(fname))(phi::Array{Float64,1},xyz::Array{Float64,1})

                phi[:]=[1.0, xyz[1]-0.5, xyz[2]-0.5, (xyz[1]-0.5)*(xyz[2]-0.5), (xyz[1]-0.5)^2, (xyz[2]-0.5)^2, (xyz[1]-0.5)^2*(xyz[2]-0.5), (xyz[1]-0.5)*(xyz[2]-0.5)^2, (xyz[1]-0.5)^3, (xyz[2]-0.5)^3]

                return nothing
            end
        end

        if rSV==:VecDG1S
            function $(esc(fname))(phi::Array{Float64,2},xyz::Array{Float64,1})

                phi[:]=[(1-xyz[1])*(1-xyz[2])   xyz[1]*(1-xyz[2])     xyz[1]*xyz[2]   (1-xyz[1])*xyz[2]     0.0                     0.0                   0.0             0.0                   0.0                     0.0                   0.0             0.0;
                        0.0                     0.0                   0.0             0.0                   (1-xyz[1])*(1-xyz[2])   xyz[1]*(1-xyz[2])     xyz[1]*xyz[2]   (1-xyz[1])*xyz[2]     0.0                     0.0                   0.0             0.0;
                        0.0                     0.0                   0.0             0.0                   0.0                     0.0                   0.0             0.0                   (1-xyz[1])*(1-xyz[2])   xyz[1]*(1-xyz[2])     xyz[1]*xyz[2]   (1-xyz[1])*xyz[2]]
                #=
                phi[:]=[(1-xyz[1])*(1-xyz[2]) 0.0    0.0    xyz[1]*(1-xyz[2])  0.0    0.0    xyz[1]*xyz[2]  0.0    0.0    (1-xyz[1])*xyz[2]  0.0    0.0;
                     0.0   (1-xyz[1])*(1-xyz[2])  0.0    0.0    xyz[1]*(1-xyz[2])  0.0    0.0    xyz[1]*xyz[2]  0.0    0.0    (1-xyz[1])*xyz[2]  0.0;
                     0.0   0.0    (1-xyz[1])*(1-xyz[2])  0.0    0.0    xyz[1]*(1-xyz[2])  0.0    0.0    xyz[1]*xyz[2]  0.0    0.0    (1-xyz[1])*xyz[2]]
                =#
                return nothing
            end
        elseif rSV==:VecDGLinS
            function $(esc(fname))(phi::Array{Float64,2},xyz::Array{Float64,1})

                phi[:]=[1.0   xyz[1]-0.5     xyz[2]-0.5      0.0      0.0            0.0             0.0      0.0            0.0;
                        0.0   0.0            0.0             1.0      xyz[1]-0.5     xyz[2]-0.5      0.0      0.0            0.0;
                        0.0   0.0            0.0             0.0      0.0            0.0             1.0      xyz[1]-0.5     xyz[2]-0.5]

                return nothing
            end
        elseif rSV==:VecDGQuadS
            function $(esc(fname))(phi::Array{Float64,2},xyz::Array{Float64,1})

                phi[:]=[1.0   xyz[1]-0.5     xyz[2]-0.5      (xyz[1]-0.5)*(xyz[2]-0.5)    (xyz[1]-0.5)^2      (xyz[2]-0.5)^2    0.0      0.0            0.0             0.0                          0.0                 0.0             0.0      0.0            0.0            0.0                          0.0                 0.0;
                        0.0   0.0            0.0             0.0                          0.0                  0.0              1.0      xyz[1]-0.5     xyz[2]-0.5      (xyz[1]-0.5)*(xyz[2]-0.5)    (xyz[1]-0.5)^2      (xyz[2]-0.5)^2  0.0      0.0            0.0            0.0                          0.0                 0.0;
                        0.0   0.0            0.0             0.0                          0.0                  0.0              0.0      0.0            0.0             0.0                          0.0                 0.0             1.0      xyz[1]-0.5     xyz[2]-0.5     (xyz[1]-0.5)*(xyz[2]-0.5)    (xyz[1]-0.5)^2      (xyz[2]-0.5)^2]

                return nothing
            end
        elseif rSV==:VecDGTriS
            function $(esc(fname))(phi::Array{Float64,2},xyz::Array{Float64,1})

                phi[:]=[1.0   xyz[1]-0.5     xyz[2]-0.5      (xyz[1]-0.5)*(xyz[2]-0.5)    (xyz[1]-0.5)^2      (xyz[2]-0.5)^2    (xyz[1]-0.5)^2*(xyz[2]-0.5)   (xyz[1]-0.5)*(xyz[2]-0.5)^2   (xyz[1]-0.5)^3      (xyz[2]-0.5)^3  0.0      0.0            0.0             0.0                          0.0                 0.0                0.0                           0.0                           0.0                 0.0                 0.0      0.0            0.0             0.0                          0.0                 0.0                0.0                           0.0                           0.0                 0.0;
                        0.0   0.0            0.0             0.0                          0.0                  0.0              0.0                           0.0                           0.0                 0.0             1.0      xyz[1]-0.5     xyz[2]-0.5      (xyz[1]-0.5)*(xyz[2]-0.5)    (xyz[1]-0.5)^2      (xyz[2]-0.5)^2     (xyz[1]-0.5)^2*(xyz[2]-0.5)   (xyz[1]-0.5)*(xyz[2]-0.5)^2   (xyz[1]-0.5)^3      (xyz[2]-0.5)^3      0.0      0.0            0.0             0.0                          0.0                 0.0                0.0                           0.0                           0.0                 0.0;
                        0.0   0.0            0.0             0.0                          0.0                  0.0              0.0                           0.0                           0.0                 0.0             0.0      0.0            0.0             0.0                          0.0                 0.0                0.0                           0.0                           0.0                 0.0                 1.0      xyz[1]-0.5     xyz[2]-0.5      (xyz[1]-0.5)*(xyz[2]-0.5)    (xyz[1]-0.5)^2      (xyz[2]-0.5)^2     (xyz[1]-0.5)^2*(xyz[2]-0.5)   (xyz[1]-0.5)*(xyz[2]-0.5)^2   (xyz[1]-0.5)^3      (xyz[2]-0.5)^3]

                return nothing
            end
        elseif rSV==:VecDG1
            function $(esc(fname))(phi::Array{Float64,2},xyz::Array{Float64,1})

                phi[:]=[(1-xyz[1])*(1-xyz[2])   xyz[1]*(1-xyz[2])     xyz[1]*xyz[2]   (1-xyz[1])*xyz[2]     0.0                     0.0                   0.0             0.0              ;
                        0.0                     0.0                   0.0             0.0                   (1-xyz[1])*(1-xyz[2])   xyz[1]*(1-xyz[2])     xyz[1]*xyz[2]   (1-xyz[1])*xyz[2]]
                #=
                phi[:]=[(1-xyz[1])*(1-xyz[2]) 0.0    0.0    xyz[1]*(1-xyz[2])  0.0    0.0    xyz[1]*xyz[2]  0.0    0.0    (1-xyz[1])*xyz[2]  0.0    0.0;
                     0.0   (1-xyz[1])*(1-xyz[2])  0.0    0.0    xyz[1]*(1-xyz[2])  0.0    0.0    xyz[1]*xyz[2]  0.0    0.0    (1-xyz[1])*xyz[2]  0.0;
                     0.0   0.0    (1-xyz[1])*(1-xyz[2])  0.0    0.0    xyz[1]*(1-xyz[2])  0.0    0.0    xyz[1]*xyz[2]  0.0    0.0    (1-xyz[1])*xyz[2]]
                =#
                return nothing
            end
        elseif rSV==:VecDGLin
            function $(esc(fname))(phi::Array{Float64,2},xyz::Array{Float64,1})

                phi[:]=[1.0   xyz[1]-0.5     xyz[2]-0.5      0.0      0.0            0.0;
                        0.0   0.0            0.0             1.0      xyz[1]-0.5     xyz[2]-0.5]

                return nothing
            end
        elseif rSV==:VecDGQuad
            function $(esc(fname))(phi::Array{Float64,2},xyz::Array{Float64,1})

                phi[:]=[1.0   xyz[1]-0.5     xyz[2]-0.5      (xyz[1]-0.5)*(xyz[2]-0.5)    (xyz[1]-0.5)^2      (xyz[2]-0.5)^2    0.0      0.0            0.0             0.0                          0.0                 0.0;
                        0.0   0.0            0.0             0.0                          0.0                  0.0              1.0      xyz[1]-0.5     xyz[2]-0.5      (xyz[1]-0.5)*(xyz[2]-0.5)    (xyz[1]-0.5)^2      (xyz[2]-0.5)^2]

                return nothing
            end
        elseif rSV==:VecDGTri
            function $(esc(fname))(phi::Array{Float64,2},xyz::Array{Float64,1})

                phi[:]=[1.0   xyz[1]-0.5     xyz[2]-0.5      (xyz[1]-0.5)*(xyz[2]-0.5)    (xyz[1]-0.5)^2      (xyz[2]-0.5)^2    (xyz[1]-0.5)^2*(xyz[2]-0.5)   (xyz[1]-0.5)*(xyz[2]-0.5)^2   (xyz[1]-0.5)^3      (xyz[2]-0.5)^3  0.0      0.0            0.0             0.0                          0.0                 0.0                0.0                           0.0                           0.0                 0.0;
                        0.0   0.0            0.0             0.0                          0.0                  0.0              0.0                           0.0                           0.0                 0.0             1.0      xyz[1]-0.5     xyz[2]-0.5      (xyz[1]-0.5)*(xyz[2]-0.5)    (xyz[1]-0.5)^2      (xyz[2]-0.5)^2     (xyz[1]-0.5)^2*(xyz[2]-0.5)   (xyz[1]-0.5)*(xyz[2]-0.5)^2   (xyz[1]-0.5)^3      (xyz[2]-0.5)^3]

                return nothing
            end
        end
    end
end



function getPhiRecoveryLength(recoverySpace::Symbol)
    if recoverySpace==:DG1
        nW=4
    elseif recoverySpace==:DGLin
        nW=3;
    elseif recoverySpace==:DGQuad
        nW=6;
    elseif recoverySpace==:DGTri
        nW=10;
    elseif recoverySpace==:VecDG1S
        nW=12
    elseif recoverySpace==:VecDGLinS
        nW=9;
    elseif recoverySpace==:VecDGQuadS
        nW=18;
    elseif recoverySpace==:VecDGTriS
        nW=30;
    elseif recoverySpace==:VecDG1
        nW=8
    elseif recoverySpace==:VecDGLin
        nW=6;
    elseif recoverySpace==:VecDGQuad
        nW=12;
    elseif recoverySpace==:VecDGTri
        nW=20;
    end
    return nW
 end
