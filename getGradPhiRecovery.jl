function getGradPhiRecovery(m::Array{Float64,1},recoverySpace::Val{:R1})

    null(xyz) = 0.0
    phi1(xyz) = 1.0;

    gradphi=[null  phi1   null;
             null  null   phi1];

    return gradphi
end

function getGradPhiRecovery(m::Array{Float64,1},recoverySpace::Val{:R1S})

    null(xyz) = 0.0
    phi1(xyz) = 1.0;

    gradphi=[null  phi1   null    null;
             null  null   phi1    null;
             null  null   null    phi1];
    return gradphi
end

function getGradPhiRecovery(m::Array{Float64,1},recoverySpace::Val{:R2})

    null(xyz) = 0.0
    phi1(xyz) = 1.0;

    phix(xyz) = (xyz[1]-m[1]);
    phiy(xyz) = (xyz[2]-m[2]);

    phi2x(xyz) = 2*(xyz[1]-m[1]);
    phi2y(xyz) = 2*(xyz[2]-m[2]);

    gradphi=[null  phi1   null    phiy   phi2x   null;
             null  null   phi1    phix   null    phi2y];

    return gradphi
end

function getGradPhiRecovery(m::Array{Float64,1},recoverySpace::Val{:R2S})

    null(xyz) = 0.0
    phi1(xyz) = 1.0;

    phix(xyz) = (xyz[1]-m[1]);
    phiy(xyz) = (xyz[2]-m[2]);
    phiz(xyz) = (xyz[3]-m[3]);

    phi2x(xyz) = 2*(xyz[1]-m[1]);
    phi2y(xyz) = 2*(xyz[2]-m[2]);
    phi2z(xyz) = 2*(xyz[3]-m[3]);

    gradphi=[null  phi1   null    null    phiy   phiz   null   phi2x   null     null;
             null  null   phi1    null    phix   null   phiz   null    phi2y    null;
             null  null   null    phi1    null   phix   phiy   null    null     phi2z]

    return gradphi
end
