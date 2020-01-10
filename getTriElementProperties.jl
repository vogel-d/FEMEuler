function getTriElementProperties(type::Symbol)
    null(x,y)=0.0;
    f1(x,y)=1.0;
    f1s(x,y)=1.0/sqrt(2);
    f2(x,y)=2.0;
    f2s(x,y)=2.0/sqrt(2);
    fx(x,y)=x;
    fmx(x,y)=-x;
    fsx(x,y)=(0.7071067811865475244)*x;
    fy(x,y)=y;
    fmy(x,y)=-y;
    fsy(x,y)=(0.7071067811865475244)*y;
    fxm1(x,y)=x-1;
    fym1(x,y)=y-1;
    f1mxy(x,y)=1-x*y;


    if type==:DG0
        phi=[f1];
        #c=[1/3 0.0; 1/3 0.0];
        #c=c[:,1:1];
        divphi=[null];
        gradphi=reshape([null,null],2,1);

        nFace=1;
        nEdge=0;
        nVert=0;

        cm=[0.5 0.0 0.0; 0.5 0.5 0.0; 0.0 0.5 0.0];
        #comp=[0];

    elseif type==:RT0

        phi=[fxm1 fx fsx; fy fym1 fsy];
        #c=[0.5 0.0 0.5; 0.0 0.5 0.5];

        divphi=[f2, f2, f2s];
        gradphi=[f1 null f1 null f1s null; null f1 null f1 null f1s];

        nFace=0;
        nEdge=1;
        nVert=0;

        cm=[0.5 0.0 0.0 1.0 0.0;
            0.5 0.5 0.0 0.0 1.0;
            0.0 0.5 1.0 0.0 0.0];
        #comp=[0];


    elseif type==:P1
        phi=[f1mxy fx fy];
        #c=[0.0 1.0 0.0; 0.0 0.0 1.0];

        divphi=[null null null];

        gradphi=[fmy f1 null;
                 fmx null f1];

        nFace=0;
        nEdge=0;
        nVert=1;

        cm=[0.5 0.0 1.0 1.0 0.0;
            0.5 0.5 0.0 1.0 1.0;
            0.0 0.5 1.0 0.0 1.0];

        #comp=[0, 0, 0];
    end
    return phi, divphi, gradphi, cm, nFace, nEdge, nVert;
end
