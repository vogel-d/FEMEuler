function getRecoveryElementProperties(type::Symbol)

    null(xyz,m) = 0.0

    #Ansatzfunktionen
    phi1(xyz,m) = 1.0;

    phix(xyz,m) = (xyz[1]-m[1])^2;
    phiy(xyz,m) = (xyz[2]-m[2])^2;
    phiz(xyz,m) = (xyz[3]-m[3])^2;

    phixy(xyz,m) = (xyz[1]-m[1])^2*(xyz[2]-m[2])^2
    phixz(xyz,m) = (xyz[1]-m[1])^2*(xyz[3]-m[3])^2
    phiyz(xyz,m) = (xyz[2]-m[2])^2*(xyz[3]-m[3])^2

    phix2(xyz,m) = (xyz[1]-m[1])^4
    phiy2(xyz,m) = (xyz[2]-m[2])^4
    phiz2(xyz,m) = (xyz[3]-m[3])^4

    #Ableitungen
    dxphix(xyz,m) = 2*(xyz[1]-m[1]);
    dyphiy(xyz,m) = 2*(xyz[2]-m[2]);
    dzphiz(xyz,m) = 2*(xyz[3]-m[3]);

    dxphixy(xyz,m) = 2*(xyz[1]-m[1])*(xyz[2]-m[2])^2
    dyphixy(xyz,m) = 2*(xyz[1]-m[1])^2*(xyz[2]-m[2])
    dxphixz(xyz,m) = 2*(xyz[1]-m[1])*(xyz[3]-m[3])^2
    dzphixz(xyz,m) = 2*(xyz[1]-m[1])^2*(xyz[3]-m[3])
    dyphiyz(xyz,m) = 2*(xyz[2]-m[2])*(xyz[3]-m[3])^2
    dzphiyz(xyz,m) = 2*(xyz[2]-m[2])^2*(xyz[3]-m[3])

    dxphix2(xyz,m) = 4*(xyz[1]-m[1])^3
    dyphiy2(xyz,m) = 4*(xyz[2]-m[2])^3
    dzphiz2(xyz,m) = 4*(xyz[3]-m[3])^3


    if type==:R1
        phi=[phi1, phix, phiy, phiz]

        divphi=[null, null, null, null];

        gradphi=[null  dxphix   null    null;
                 null  null     dyphiy  null;
                 null  null     null    dzphiz];

        nFace=4;
        nEdge=0;
        nVert=0;

        cm=Dict([1,2]=>[0,0,0,0], [2,3]=>[0,0,0,0], [3,4]=>[0,0,0,0], [1,4]=>[0,0,0,0]);

    elseif type==:R2
        phi=[phi1 phix phiy phiz phixy phixz phiyz phix2 phiy2 phiz2]

        divphi=[null, null, null, null, null, null, null, null, null, null];

        gradphi=[null  dxphix   null    null    dxphixy   dxphixz   null      dxphix2   null       null;
                 null  null     dyphiy  null    dyphixy   null      dyphiyz   null      dyphiy1    null;
                 null  null     null    dzphiz  null      dzphixz   dzphiyz   null      null       dzphiz2]

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
