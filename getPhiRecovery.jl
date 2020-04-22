function getPhiRecovery(m::Array{Float64,1},order::Val{1})

    #null(xyz) = 0.0
    phi1(xyz) = 1.0;

    phix(xyz) = (xyz[1]-m[1]);
    phiy(xyz) = (xyz[2]-m[2]);
    phiz(xyz) = (xyz[3]-m[3]);

    phi=[phi1, phix, phiy, phiz]
    #=
    phi=[phi1 phix phiy phiz null null null null null null null null;
         null null null null phi1 phix phiy phiz null null null null;
         null null null null null null null null phi1 phix phiy phiz]
    =#
    return phi
end

function getPhiRecovery(m::Array{Float64,1},order::Val{2})

    #null(xyz) = 0.0
    phi1(xyz) = 1.0;

    phix(xyz) = (xyz[1]-m[1]);
    phiy(xyz) = (xyz[2]-m[2]);
    phiz(xyz) = (xyz[3]-m[3]);

    phixy(xyz) = (xyz[1]-m[1])*(xyz[2]-m[2])
    phixz(xyz) = (xyz[1]-m[1])*(xyz[3]-m[3])
    phiyz(xyz) = (xyz[2]-m[2])*(xyz[3]-m[3])

    phix2(xyz) = (xyz[1]-m[1])^2
    phiy2(xyz) = (xyz[2]-m[2])^2
    phiz2(xyz) = (xyz[3]-m[3])^2

    phi=[phi1 phix phiy phiz phixy phixz phiyz phix2 phiy2 phiz2]
    #=
    phi=[phi1 phix phiy phiz phixy phixz phiyz phix2 phiy2 phiz2 null null null null null null null null null null null null null null null null null null null null;
         null null null null null null null null null null phi1 phix phiy phiz phixy phixz phiyz phix2 phiy2 phiz2 null null null null null null null null null null;
         null null null null null null null null null null null null null null null null null null null null phi1 phix phiy phiz phixy phixz phiyz phix2 phiy2 phiz2]
    =#
    return phi
end
