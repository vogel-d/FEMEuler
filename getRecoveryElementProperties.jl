function getRecoveryElementProperties(type::Symbol)

    null(xyz,m) = 0.0

    #Ansatzfunktionen
    phi1(xyz,m) = 1.0;

    phix(xyz,m) = (xyz[1]-m[1]);
    phiy(xyz,m) = (xyz[2]-m[2]);
    phiz(xyz,m) = (xyz[3]-m[3]);

    phixy(xyz,m) = (xyz[1]-m[1])*(xyz[2]-m[2])
    phixz(xyz,m) = (xyz[1]-m[1])*(xyz[3]-m[3])
    phiyz(xyz,m) = (xyz[2]-m[2])*(xyz[3]-m[3])

    phix2(xyz,m) = (xyz[1]-m[1])^2
    phiy2(xyz,m) = (xyz[2]-m[2])^2
    phiz2(xyz,m) = (xyz[3]-m[3])^2

    #Ableitungen
    phi2x(xyz,m) = 2*(xyz[1]-m[1]);
    phi2y(xyz,m) = 2*(xyz[2]-m[2]);
    phi2z(xyz,m) = 2*(xyz[3]-m[3]);

    if type==:R1
        phi=[phi1, phix, phiy, phiz]

        divphi=[null, null, null, null];

        gradphi=[null  phi1   null    null;
                 null  null   phi1    null;
                 null  null   null    phi1];

        nFace=4;
        nEdge=0;
        nVert=0;

        cm=Dict([1,2]=>[0,0,0,0], [2,3]=>[0,0,0,0], [3,4]=>[0,0,0,0], [1,4]=>[0,0,0,0]);

    elseif type==:R2
        phi=[phi1 phix phiy phiz phixy phixz phiyz phix2 phiy2 phiz2]

        divphi=[null, null, null, null, null, null, null, null, null, null];

        gradphi=[null  phi1   null    null    phiy   phiz   null   phi2x   null     null;
                 null  null   phi1    null    phix   null   phiz   null    phi2y    null;
                 null  null   null    phi1    null   phix   phiy   null    null     phi2z]

        nFace=10;
        nEdge=0;
        nVert=0;

        cm=Dict([1,2]=>[0,0,0,0,0,0,0,0,0,0], [2,3]=>[0,0,0,0,0,0,0,0,0,0],
            [3,4]=>[0,0,0,0,0,0,0,0,0,0], [1,4]=>[0,0,0,0,0,0,0,0,0,0]);

    else
        error("Unzulässiger finite-Elemente-Raum");
    end
    return phi, divphi, gradphi, cm, nFace, nEdge, nVert;
end
