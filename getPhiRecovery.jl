function getPhiRecovery(recoverySpace::Val{:R1})

    #null(xyz) = 0.0
    phi1(xyz) = 1.0;

    phix(xyz) = xyz[1];
    phiy(xyz) = xyz[2];

    phi=[phi1, phix, phiy]

    return phi
end

function getPhiRecovery(recoverySpace::Val{:R1S})

    #null(xyz) = 0.0
    phi1(xyz) = 1.0;

    phix(xyz) = xyz[1];
    phiy(xyz) = xyz[2];
    phiz(xyz) = xyz[3];

    phi=[phi1, phix, phiy, phiz]

    return phi
end

function getPhiRecovery(recoverySpace::Val{:R2})

    #null(xyz) = 0.0
    phi1(xyz) = 1.0;

    phix(xyz) = xyz[1];
    phiy(xyz) = xyz[2];

    phixy(xyz) = xyz[1]*xyz[2]

    phix2(xyz) = xyz[1]^2
    phiy2(xyz) = xyz[2]^2

    phi=[phi1, phix, phiy, phixy, phix2, phiy2]

    return phi
end

function getPhiRecovery(recoverySpace::Val{:R2S})

    #null(xyz) = 0.0
    phi1(xyz) = 1.0;

    phix(xyz) = xyz[1];
    phiy(xyz) = xyz[2];
    phiz(xyz) = xyz[3];

    phixy(xyz) = xyz[1]*xyz[2]
    phixz(xyz) = xyz[1]*xyz[3]
    phiyz(xyz) = xyz[2]*xyz[3]

    phix2(xyz) = xyz[1]^2
    phiy2(xyz) = xyz[2]^2
    phiz2(xyz) = xyz[3]^2

    phi=[phi1, phix, phiy, phiz, phixy, phixz, phiyz, phix2, phiy2, phiz2]

    return phi
end


function getPhiRecovery(recoverySpace::Val{:VecR1})

    null(xyz) = 0.0
    phi1(xyz) = 1.0;

    phix(xyz) = xyz[1];
    phiy(xyz) = xyz[2];

    phi=[phi1 phix phiy null null null;
         null null null phi1 phix phiy]

    return phi
end

function getPhiRecovery(recoverySpace::Val{:VecR1S})

    null(xyz) = 0.0
    phi1(xyz) = 1.0;

    phix(xyz) = xyz[1];
    phiy(xyz) = xyz[2];
    phiz(xyz) = xyz[3];

    phi=[phi1 phix phiy phiz null null null null null null null null;
         null null null null phi1 phix phiy phiz null null null null;
         null null null null null null null null phi1 phix phiy phiz]

    return phi
end

function getPhiRecovery(recoverySpace::Val{:VecR2})

    #null(xyz) = 0.0
    phi1(xyz) = 1.0;

    phix(xyz) = xyz[1];
    phiy(xyz) = xyz[2];

    phixy(xyz) = xyz[1]*xyz[2]

    phix2(xyz) = xyz[1]^2
    phiy2(xyz) = xyz[2]^2

    phi=[phi1 phix phiy phixy phix2 phiy2 null null null null null null;
         null null null null null null phi1 phix phiy phixy phix2 phiy2]

    return phi
end

function getPhiRecovery(recoverySpace::Val{:VecR2S})

    #null(xyz) = 0.0
    phi1(xyz) = 1.0;

    phix(xyz) = xyz[1];
    phiy(xyz) = xyz[2];
    phiz(xyz) = xyz[3];

    phixy(xyz) = xyz[1]*xyz[2]
    phixz(xyz) = xyz[1]*xyz[3]
    phiyz(xyz) = xyz[2]*xyz[3]

    phix2(xyz) = xyz[1]^2
    phiy2(xyz) = xyz[2]^2
    phiz2(xyz) = xyz[3]^2

    phi=[phi1 phix phiy phiz phixy phixz phiyz phix2 phiy2 phiz2 null null null null null null null null null null null null null null null null null null null null;
         null null null null null null null null null null phi1 phix phiy phiz phixy phixz phiyz phix2 phiy2 phiz2 null null null null null null null null null null;
         null null null null null null null null null null null null null null null null null null null null phi1 phix phiy phiz phixy phixz phiyz phix2 phiy2 phiz2]

    return phi
end
