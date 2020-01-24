function getTriElementProperties(type::Symbol)
    null(x,y)=0.0;
    f1(x,y)=1.0;
    fm1(x,y)=-1.0;
    f2(x,y)=2.0;
    fx(x,y)=x;
    fmx(x,y)=-x;
    
    f1s(x,y)=1.0/sqrt(2);
    f2s(x,y)=2.0/sqrt(2);
    fsx(x,y)=(0.7071067811865475244)*x;
    fsy(x,y)=(0.7071067811865475244)*y;

    #=
    f1s(x,y)=1.4142135623730951;
    f2s(x,y)=2.0*1.4142135623730951;
    fsx(x,y)=1.4142135623730951*x;
    fsy(x,y)=1.4142135623730951*y;
    =#
    fy(x,y)=y;
    fmy(x,y)=-y;
    fxm1(x,y)=x-1;
    fym1(x,y)=y-1;
    f1mxy(x,y)=1-x*y;
    f1mxmy(x,y)=1-x-y;
    f1mx(x,y)=1-x;
    f1my(x,y)=1-y;



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


    elseif type==:RT0

        #phi=[fxm1 fx fsx;
        #     fy fym1 fsy];
         phi=[fx fsx fxm1;
              fym1 fsy fy];
        #c=[0.5 0.0 0.5; 0.0 0.5 0.5];

        #divphi=[f2, f2, f2s];
        #gradphi=[f1 null f1 null f1s null;
        #         null f1 null f1 null f1s];

        divphi=[f2, f2s, f2];
        gradphi=[f1 null f1s null f1 null;
                 null f1 null f1s null f1];

        nFace=0;
        nEdge=1;
        nVert=0;

        cm=[0.5 0.0 1.0 0.0 0.0;
            0.5 0.5 0.0 1.0 0.0;
            0.0 0.5 0.0 0.0 1.0];



    elseif type==:P1
        phi=[f1mxmy, fx, fy];
        #c=[0.0 1.0 0.0; 0.0 0.0 1.0];

        divphi=[null, null, null];

        gradphi=[fm1 f1 null;
                 fm1 null f1];

        nFace=0;
        nEdge=0;
        nVert=1;

        cm=[0.5 0.0 1.0 1.0 0.0;
            0.5 0.5 0.0 1.0 1.0;
            0.0 0.5 1.0 0.0 1.0];


    end
    return phi, divphi, gradphi, cm, nFace, nEdge, nVert;
end
