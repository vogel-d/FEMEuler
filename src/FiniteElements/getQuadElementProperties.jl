function getQuadElementProperties(type::Symbol)

    null(x,y)=0.0;

    g0(x)=1.0;
    g10(x)=(x-1.0)/(0.0-1.0);
    g11(x)=(x-0.0)/(1.0-0.0);
    g20(x)=(x-0.5)*(x-1.0)/((0.0-0.5)*(0.0-1.0));
    g21(x)=(x-0.0)*(x-1.0)/((0.5-0.0)*(0.5-1.0));
    g22(x)=(x-0.0)*(x-0.5)/((1.0-0.0)*(1.0-0.5));

    Dg0(x)=0.0;
    Dg10(x)=1.0/(0.0-1.0);
    Dg11(x)=1.0/(1.0-0.0);
    Dg20(x)=((x-0.5)+(x-1.0))/((0.0-0.5)*(0.0-1.0));
    Dg21(x)=((x-0.0)+(x-1.0))/((0.5-0.0)*(0.5-1.0));
    Dg22(x)=((x-0.0)+(x-0.5))/((1.0-0.0)*(1.0-0.5));

    h0_0(x,y)=g0(x)*g0(y);
    Dxh0_0(x,y)=Dg0(x)*g0(y);
    Dyh0_0(x,y)=g0(x)*Dg0(y);

    h10_0(x,y)=g10(x)*g0(y);
    h11_0(x,y)=g11(x)*g0(y);
    h0_10(x,y)=g0(x)*g10(y);
    h0_11(x,y)=g0(x)*g11(y);
    Dxh10_0(x,y)=Dg10(x)*g0(y);
    Dxh11_0(x,y)=Dg11(x)*g0(y);
    Dxh0_10(x,y)=Dg0(x)*g10(y);
    Dxh0_11(x,y)=Dg0(x)*g11(y);
    Dyh10_0(x,y)=g10(x)*Dg0(y);
    Dyh11_0(x,y)=g11(x)*Dg0(y);
    Dyh0_10(x,y)=g0(x)*Dg10(y);
    Dyh0_11(x,y)=g0(x)*Dg11(y);

    h10_0neg(x,y)=(-1)*g10(x)*g0(y);
    h11_0neg(x,y)=(-1)*g11(x)*g0(y);
    h0_10neg(x,y)=(-1)*g0(x)*g10(y);
    h0_11neg(x,y)=(-1)*g0(x)*g11(y);
    Dxh10_0neg(x,y)=(-1)*Dg10(x)*g0(y);
    Dxh11_0neg(x,y)=(-1)*Dg11(x)*g0(y);
    Dxh0_10neg(x,y)=(-1)*Dg0(x)*g10(y);
    Dxh0_11neg(x,y)=(-1)*Dg0(x)*g11(y);
    Dyh10_0neg(x,y)=(-1)*g10(x)*Dg0(y);
    Dyh11_0neg(x,y)=(-1)*g11(x)*Dg0(y);
    Dyh0_10neg(x,y)=(-1)*g0(x)*Dg10(y);
    Dyh0_11neg(x,y)=(-1)*g0(x)*Dg11(y);

    h10_10(x,y)=g10(x)*g10(y);
    h10_11(x,y)=g10(x)*g11(y);
    h11_10(x,y)=g11(x)*g10(y);
    h11_11(x,y)=g11(x)*g11(y);
    Dxh10_10(x,y)=Dg10(x)*g10(y);
    Dxh10_11(x,y)=Dg10(x)*g11(y);
    Dxh11_10(x,y)=Dg11(x)*g10(y);
    Dxh11_11(x,y)=Dg11(x)*g11(y);
    Dyh10_10(x,y)=g10(x)*Dg10(y);
    Dyh10_11(x,y)=g10(x)*Dg11(y);
    Dyh11_10(x,y)=g11(x)*Dg10(y);
    Dyh11_11(x,y)=g11(x)*Dg11(y);

    h20_20(x,y)=g20(x)*g20(y);
    h21_20(x,y)=g21(x)*g20(y);
    h22_20(x,y)=g22(x)*g20(y);
    h20_21(x,y)=g20(x)*g21(y);
    h21_21(x,y)=g21(x)*g21(y);
    h22_21(x,y)=g22(x)*g21(y);
    h20_22(x,y)=g20(x)*g22(y);
    h21_22(x,y)=g21(x)*g22(y);
    h22_22(x,y)=g22(x)*g22(y);
    Dxh20_20(x,y)=Dg20(x)*g20(y);
    Dxh21_20(x,y)=Dg21(x)*g20(y);
    Dxh22_20(x,y)=Dg22(x)*g20(y);
    Dxh20_21(x,y)=Dg20(x)*g21(y);
    Dxh21_21(x,y)=Dg21(x)*g21(y);
    Dxh22_21(x,y)=Dg22(x)*g21(y);
    Dxh20_22(x,y)=Dg20(x)*g22(y);
    Dxh21_22(x,y)=Dg21(x)*g22(y);
    Dxh22_22(x,y)=Dg22(x)*g22(y);
    Dyh20_20(x,y)=g20(x)*Dg20(y);
    Dyh21_20(x,y)=g21(x)*Dg20(y);
    Dyh22_20(x,y)=g22(x)*Dg20(y);
    Dyh20_21(x,y)=g20(x)*Dg21(y);
    Dyh21_21(x,y)=g21(x)*Dg21(y);
    Dyh22_21(x,y)=g22(x)*Dg21(y);
    Dyh20_22(x,y)=g20(x)*Dg22(y);
    Dyh21_22(x,y)=g21(x)*Dg22(y);
    Dyh22_22(x,y)=g22(x)*Dg22(y);

    h20_10(x,y)=g20(x)*g10(y);
    h21_10(x,y)=g21(x)*g10(y);
    h22_10(x,y)=g22(x)*g10(y);
    h20_11(x,y)=g20(x)*g11(y);
    h21_11(x,y)=g21(x)*g11(y);
    h22_11(x,y)=g22(x)*g11(y);

    h10_20(x,y)=g10(x)*g20(y);
    h11_20(x,y)=g11(x)*g20(y);
    h10_21(x,y)=g10(x)*g21(y);
    h11_21(x,y)=g11(x)*g21(y);
    h10_22(x,y)=g10(x)*g22(y);
    h11_22(x,y)=g11(x)*g22(y);

    Dxh20_10(x,y)=Dg20(x)*g10(y);
    Dxh21_10(x,y)=Dg21(x)*g10(y);
    Dxh22_10(x,y)=Dg22(x)*g10(y);
    Dxh20_11(x,y)=Dg20(x)*g11(y);
    Dxh21_11(x,y)=Dg21(x)*g11(y);
    Dxh22_11(x,y)=Dg22(x)*g11(y);

    Dyh20_10(x,y)=g20(x)*Dg10(y);
    Dyh21_10(x,y)=g21(x)*Dg10(y);
    Dyh22_10(x,y)=g22(x)*Dg10(y);
    Dyh20_11(x,y)=g20(x)*Dg11(y);
    Dyh21_11(x,y)=g21(x)*Dg11(y);
    Dyh22_11(x,y)=g22(x)*Dg11(y);

    Dxh10_20(x,y)=Dg10(x)*g20(y);
    Dxh11_20(x,y)=Dg11(x)*g20(y);
    Dxh10_21(x,y)=Dg10(x)*g21(y);
    Dxh11_21(x,y)=Dg11(x)*g21(y);
    Dxh10_22(x,y)=Dg10(x)*g22(y);
    Dxh11_22(x,y)=Dg11(x)*g22(y);

    Dyh10_20(x,y)=g10(x)*Dg20(y);
    Dyh11_20(x,y)=g11(x)*Dg20(y);
    Dyh10_21(x,y)=g10(x)*Dg21(y);
    Dyh11_21(x,y)=g11(x)*Dg21(y);
    Dyh10_22(x,y)=g10(x)*Dg22(y);
    Dyh11_22(x,y)=g11(x)*Dg22(y);

    f1(x,y)=1.0
    fx(x,y)=x-0.5
    fy(x,y)=y-0.5
    fxy(x,y)=(x-0.5)*(y-0.5)
    fx2(x,y)=(x-0.5)^2
    fy2(x,y)=(y-0.5)^2
    fx2y(x,y)=(x-0.5)^2*(y-0.5)
    fxy2(x,y)=(x-0.5)*(y-0.5)^2
    fx3(x,y)=(x-0.5)^3
    fy3(x,y)=(y-0.5)^3
    dxfx2(x,y)=2*(x-0.5)
    dyfy2(x,y)=2*(y-0.5)
    dxfx2y(x,y)=2*(x-0.5)*(y-0.5)
    dyfxy2(x,y)=(x-0.5)*2*(y-0.5)
    dxfx3(x,y)=3*(x-0.5)^2
    dyfy3(x,y)=3*(y-0.5)^2

    if type==:DG0
        phi=[h0_0];
        divphi=[null];
        gradphi=reshape([Dxh0_0,Dxh0_0],2,1);
        nFace=1;
        nEdge=0;
        nVert=0;

        cm=Dict([1,2]=>[0], [2,3]=>[0], [3,4]=>[0], [1,4]=>[0]);

    elseif type==:P1
        phi=[h10_10,h11_10,h11_11,h10_11];

        divphi=[null, null, null, null];
        gradphi=[Dxh10_10 Dxh11_10 Dxh11_11 Dxh10_11;
                 Dyh10_10 Dyh11_10 Dyh11_11 Dyh10_11];
        nFace=0;
        nEdge=0;
        nVert=1;

        cm=Dict([1,2]=>[1,1,0,0], [2,3]=>[0,1,1,0], [3,4]=>[0,0,1,1], [1,4]=>[1,0,0,1]);


    elseif type==:DG1
        phi=[h10_10,h11_10,h11_11,h10_11];

        divphi=[null, null, null, null];
        gradphi=[Dxh10_10 Dxh11_10 Dxh11_11 Dxh10_11;
                 Dyh10_10 Dyh11_10 Dyh11_11 Dyh10_11];

        nFace=4;
        nEdge=0;
        nVert=0;

        cm=Dict([1,2]=>[0,0,0,0], [2,3]=>[0,0,0,0], [3,4]=>[0,0,0,0], [1,4]=>[0,0,0,0]);

    elseif type==:DGLin
        phi=[f1,fx,fy];
        divphi=[null, null, null];
        gradphi=[null f1 null;
                 null null f1];

        nFace=3;
        nEdge=0;
        nVert=0;

        cm=Dict([1,2]=>[0,0,0], [2,3]=>[0,0,0], [3,4]=>[0,0,0], [1,4]=>[0,0,0]);

    elseif type==:DGQuad
        phi=[f1,fx,fy,fxy,fx2,fy2];
        divphi=[null, null, null, null, null, null];
        gradphi=[null f1 null fy dxfx2 null;
                 null null f1 fx null dyfy2];

        nFace=6;
        nEdge=0;
        nVert=0;

        cm=Dict([1,2]=>[0,0,0,0,0,0], [2,3]=>[0,0,0,0,0,0], [3,4]=>[0,0,0,0,0,0], [1,4]=>[0,0,0,0,0,0]);

    elseif type==:DGTri
        phi=[f1,fx,fy,fxy,fx2,fy2,fx2y,fxy2,fx3,fy3];
        divphi=[null, null, null, null, null, null, null, null, null, null];
        gradphi=[null f1 null fy dxfx2 null dxfx2y fy2    dxfx3 null;
                 null null f1 fx null dyfy2 fx2    dyfxy2 null  dyfy3];

        nFace=10;
        nEdge=0;
        nVert=0;

        cm=Dict([1,2]=>[0,0,0,0,0,0,0,0,0,0], [2,3]=>[0,0,0,0,0,0,0,0,0,0], [3,4]=>[0,0,0,0,0,0,0,0,0,0], [1,4]=>[0,0,0,0,0,0,0,0,0,0]);

    elseif type==:P2
        phi=[h21_21, h21_20, h22_21, h21_22, h20_21, h20_20, h22_20, h22_22, h20_22]

        divphi=[null, null, null, null, null, null, null, null, null];
        gradphi=[Dxh21_21 Dxh21_20 Dxh22_21 Dxh21_22 Dxh20_21 Dxh20_20 Dxh22_20 Dxh22_22 Dxh20_22;
                 Dyh21_21 Dyh21_20 Dyh22_21 Dyh21_22 Dyh20_21 Dyh20_20 Dyh22_20 Dyh22_22 Dyh20_22]

        nFace=1;
        nEdge=1;
        nVert=1;

        cm=Dict([1,2]=>[0,1,0,0,0,1,1,0,0], [2,3]=>[0,0,1,0,0,0,1,1,0], [3,4]=>[0,0,0,1,0,0,0,1,1], [1,4]=>[0,0,0,0,1,1,0,0,1]);

    elseif type==:DG2
        phi=[h20_20, h22_20, h22_22, h20_22, h21_20, h22_21, h21_22, h20_21, h21_21]

        divphi=[null, null, null, null, null, null, null, null, null];
        gradphi=[Dxh20_20 Dxh22_20 Dxh22_22 Dxh20_22 Dxh21_20 Dxh22_21 Dxh21_22 Dxh20_21 Dxh21_21;
                 Dyh20_20 Dyh22_20 Dyh22_22 Dyh20_22 Dyh21_20 Dyh22_21 Dyh21_22 Dyh20_21 Dyh21_21]

        nFace=9;
        nEdge=0;
        nVert=0;

        cm=Dict([1,2]=>[0,0,0,0,0,0,0,0,0], [2,3]=>[0,0,0,0,0,0,0,0,0], [3,4]=>[0,0,0,0,0,0,0,0,0], [1,4]=>[0,0,0,0,0,0,0,0,0]);

    elseif type==:RT0
        #=
        # Compund Finite Elements
        # TO DO: Fix bug to use these normals for regular finite elements in the plane
        phi=[null h11_0 null h10_0;
             h0_10neg null h0_11neg null];

        divphi=[Dyh0_10neg, Dxh11_0, Dyh0_11neg, Dxh10_0];
        gradphi=[null null    Dxh11_0 null null null    Dxh10_0 null;
                 null Dyh0_10neg null    null null Dyh0_11neg null    null];
        =#
        phi=[null h11_0 null h10_0;
             h0_10 null h0_11 null];

        divphi=[Dyh0_10, Dxh11_0, Dyh0_11, Dxh10_0];
        gradphi=[null null    Dxh11_0 null null null    Dxh10_0 null;
                 null Dyh0_10 null    null null Dyh0_11 null    null];

        nFace=0;
        nEdge=1;
        nVert=0;

        cm=Dict([1,2]=>[1,0,0,0], [2,3]=>[0,1,0,0], [3,4]=>[0,0,1,0], [1,4]=>[0,0,0,1]);

    elseif type==:RT0B #Broken RT0
        #=
        # Compund Finite Elements
        phi=[null h11_0 null h10_0;
             h0_10neg null h0_11neg null];

        divphi=[Dyh0_10neg, Dxh11_0, Dyh0_11neg, Dxh10_0];
        gradphi=[null null    Dxh11_0 null null null    Dxh10_0 null;
                 null Dyh0_10neg null    null null Dyh0_11neg null    null];
        =#
        phi=[null h11_0 null h10_0;
             h0_10 null h0_11 null];

        divphi=[Dyh0_10, Dxh11_0, Dyh0_11, Dxh10_0];
        gradphi=[null null    Dxh11_0 null null null    Dxh10_0 null;
                 null Dyh0_10 null    null null Dyh0_11 null    null];

        nFace=4;
        nEdge=0;
        nVert=0;

        cm=Dict([1,2]=>[0,0,0,0], [2,3]=>[0,0,0,0], [3,4]=>[0,0,0,0], [1,4]=>[0,0,0,0]);

    elseif type==:RT1

        phi=[h21_10 null   h21_11   null   null   null   h22_10 h22_11 null   null   h20_10 h20_11;
             null   h11_21 null     h10_21 h10_20 h11_20 null   null   h10_22 h11_22 null   null  ]
        divphi=[Dxh21_10 Dyh11_21 Dxh21_11 Dyh10_21 Dyh10_20 Dyh11_20 Dxh22_10 Dxh22_11 Dyh10_22 Dyh11_22 Dxh20_10 Dxh20_11]
        gradphi=Matrix(undef, 2, 24)
        gradphi[1:2, 1:2]=[Dxh21_10 Dyh21_10;
                          null     null]
        gradphi[1:2, 3: 4]=[null     null;
                            Dxh11_21 Dyh11_21]
        gradphi[1:2, 5: 6]=[Dxh21_11 Dyh21_11;
                            null     null]
        gradphi[1:2, 7: 8]=[null     null;
                            Dxh10_21 Dyh10_21]
        gradphi[1:2, 9:10]=[null     null;
                            Dxh10_20 Dyh10_20]
        gradphi[1:2,11:12]=[null     null;
                            Dxh11_20 Dyh11_20]
        gradphi[1:2,13:14]=[Dxh22_10 Dyh22_10;
                            null     null]
        gradphi[1:2,15:16]=[Dxh22_11 Dyh22_11;
                            null     null]
        gradphi[1:2,17:18]=[null     null;
                            Dxh10_22 Dyh10_22]
        gradphi[1:2,19:20]=[null     null;
                            Dxh11_22 Dyh11_22]
        gradphi[1:2,21:22]=[Dxh20_10 Dyh20_10;
                            null     null]
        gradphi[1:2,23:24]=[Dxh20_11 Dyh20_11;
                            null     null]

        nFace=4;
        nEdge=2;
        nVert=0;

        cm=Dict([1,2]=>[0,0,0,0,1,1,0,0,0,0,0,0], [2,3]=>[0,0,0,0,0,0,1,1,0,0,0,0], [3,4]=>[0,0,0,0,0,0,0,0,1,1,0,0], [1,4]=>[0,0,0,0,0,0,0,0,0,0,1,1]);

    elseif type==:RT1B

        phi=[h21_10 null   h21_11   null   null   null   h22_10 h22_11 null   null   h20_10 h20_11;
             null   h11_21 null     h10_21 h10_20 h11_20 null   null   h10_22 h11_22 null   null  ]
        divphi=[Dxh21_10 Dyh11_21 Dxh21_11 Dyh10_21 Dyh10_20 Dyh11_20 Dxh22_10 Dxh22_11 Dyh10_22 Dyh11_22 Dxh20_10 Dxh20_11]
        gradphi=Matrix(undef, 2, 24)
        gradphi[1:2, 1:2]=[Dxh21_10 Dyh21_10;
                          null     null]
        gradphi[1:2, 3: 4]=[null     null;
                            Dxh11_21 Dyh11_21]
        gradphi[1:2, 5: 6]=[Dxh21_11 Dyh21_11;
                            null     null]
        gradphi[1:2, 7: 8]=[null     null;
                            Dxh10_21 Dyh10_21]
        gradphi[1:2, 9:10]=[null     null;
                            Dxh10_20 Dyh10_20]
        gradphi[1:2,11:12]=[null     null;
                            Dxh11_20 Dyh11_20]
        gradphi[1:2,13:14]=[Dxh22_10 Dyh22_10;
                            null     null]
        gradphi[1:2,15:16]=[Dxh22_11 Dyh22_11;
                            null     null]
        gradphi[1:2,17:18]=[null     null;
                            Dxh10_22 Dyh10_22]
        gradphi[1:2,19:20]=[null     null;
                            Dxh11_22 Dyh11_22]
        gradphi[1:2,21:22]=[Dxh20_10 Dyh20_10;
                            null     null]
        gradphi[1:2,23:24]=[Dxh20_11 Dyh20_11;
                            null     null]
        nFace=12;
        nEdge=0;
        nVert=0;

        cm=Dict([1,2]=>[0,0,0,0,0,0,0,0,0,0,0,0], [2,3]=>[0,0,0,0,0,0,0,0,0,0,0,0], [3,4]=>[0,0,0,0,0,0,0,0,0,0,0,0], [1,4]=>[0,0,0,0,0,0,0,0,0,0,0,0]);

    elseif type==:VecP1

        phi=[h10_10 null h11_10 null h11_11 null h10_11 null;
             null h10_10 null h11_10 null h11_11 null h10_11]
        divphi=[Dxh10_10,Dyh10_10,Dxh11_10,Dyh11_10,Dxh11_11,Dyh11_11,Dxh10_11,Dyh10_11];
        gradphi=Matrix(undef, 2, 16)
        gradphi[1:2, 1: 2]=[Dxh10_10 Dyh10_10;
                           null     null]
        gradphi[1:2, 3: 4]=[null     null
                           Dxh10_10 Dyh10_10]
        gradphi[1:2, 5: 6]=[Dxh11_10 Dyh11_10;
                           null     null]
        gradphi[1:2, 7: 8]=[null     null
                           Dxh11_10 Dyh11_10]
        gradphi[1:2, 9:10]=[Dxh11_11 Dyh11_11;
                           null     null]
        gradphi[1:2,11:12]=[null     null
                           Dxh11_11 Dyh11_11]
        gradphi[1:2,13:14]=[Dxh10_11 Dyh10_11;
                           null     null]
        gradphi[1:2,15:16]=[null     null
                           Dxh10_11 Dyh10_11]

        nFace=0;
        nEdge=0;
        nVert=2;

        cm=Dict([1,2]=>[0,1,0,1,0,0,0,0], [2,3]=>[0,0,1,0,1,0,0,0], [3,4]=>[0,0,0,0,0,1,0,1], [1,4]=>[1,0,0,0,0,0,1,0]);

    elseif type==:VecDG1

        phi=[h10_10 null h11_10 null h11_11 null h10_11 null;
             null h10_10 null h11_10 null h11_11 null h10_11]
        divphi=[Dxh10_10,Dyh10_10,Dxh11_10,Dyh11_10,Dxh11_11,Dyh11_11,Dxh10_11,Dyh10_11];
        gradphi=Matrix(undef, 2, 16)
        gradphi[1:2, 1: 2]=[Dxh10_10 Dyh10_10;
                           null     null]
        gradphi[1:2, 3: 4]=[null     null
                           Dxh10_10 Dyh10_10]
        gradphi[1:2, 5: 6]=[Dxh11_10 Dyh11_10;
                           null     null]
        gradphi[1:2, 7: 8]=[null     null
                           Dxh11_10 Dyh11_10]
        gradphi[1:2, 9:10]=[Dxh11_11 Dyh11_11;
                           null     null]
        gradphi[1:2,11:12]=[null     null
                           Dxh11_11 Dyh11_11]
        gradphi[1:2,13:14]=[Dxh10_11 Dyh10_11;
                           null     null]
        gradphi[1:2,15:16]=[null     null
                           Dxh10_11 Dyh10_11]

        nFace=8;
        nEdge=0;
        nVert=0;

        cm=Dict([1,2]=>[0,0,0,0,0,0,0,0], [2,3]=>[0,0,0,0,0,0,0,0], [3,4]=>[0,0,0,0,0,0,0,0], [1,4]=>[0,0,0,0,0,0,0,0]);

    elseif type==:VecDGLin
        phi=[f1     fx      fy    null   null    null;
             null   null    null  f1     fx      fy];
        divphi=[null, f1, null,null,null,f1];
        gradphi=[null f1 null;
                 null null f1];
        gradphi=Matrix(undef, 2, 12)
        gradphi[1:2, 1: 2]=[null     null;
                            null     null]
        gradphi[1:2, 3: 4]=[f1       null;
                            null     null]
        gradphi[1:2, 5: 6]=[null     f1;
                            null     null]
        gradphi[1:2, 7: 8]=[null     null;
                            null     null]
        gradphi[1:2, 9:10]=[null     null;
                            f1       null]
        gradphi[1:2,11:12]=[null     null;
                            null     f1]

        nFace=6;
        nEdge=0;
        nVert=0;

        cm=Dict([1,2]=>[0,0,0,0,0,0], [2,3]=>[0,0,0,0,0,0], [3,4]=>[0,0,0,0,0,0], [1,4]=>[0,0,0,0,0,0]);

    elseif type==:VecDGQuad
        phi=[f1     fx     fy     fxy     fx2    fy2    null   null   null   null    null   null;
             null   null   null   null    null   null   f1     fx     fy     fxy    fx2     fy2];
        divphi=[null, f1, null, fy, dxfx2, null, null,null,f1,fx,null,dyfy2];
        gradphi=Matrix(undef, 2, 24)
        gradphi[1:2, 1: 2]=[null     null;
                            null     null]
        gradphi[1:2, 3: 4]=[f1       null;
                            null     null]
        gradphi[1:2, 5: 6]=[null     f1;
                            null     null]
        gradphi[1:2, 7: 8]=[fy       fx;
                            null     null]
        gradphi[1:2, 9:10]=[dxfx2    null;
                            null     null]
        gradphi[1:2,11:12]=[null     dyfy2;
                            null     null]
        gradphi[1:2,13:14]=[null     null;
                            null     null]
        gradphi[1:2,15:16]=[null     null;
                            f1       null]
        gradphi[1:2,17:18]=[null     null;
                            null     f1]
        gradphi[1:2,19:20]=[null     null;
                            fy       fx]
        gradphi[1:2,21:22]=[null     null;
                            dxfx2    null]
        gradphi[1:2,23:24]=[null     null;
                            null     dyfy2]
        nFace=12;
        nEdge=0;
        nVert=0;

        cm=Dict([1,2]=>[0,0,0,0,0,0,0,0,0,0,0,0], [2,3]=>[0,0,0,0,0,0,0,0,0,0,0,0], [3,4]=>[0,0,0,0,0,0,0,0,0,0,0,0], [1,4]=>[0,0,0,0,0,0,0,0,0,0,0,0]);

    elseif type==:VecDGTri
        phi=[f1     fx     fy     fxy     fx2    fy2    fx2y   fxy2   fx3    fy3    null   null   null   null    null   null   null   null   null   null;
             null   null   null   null    null   null   null   null   null   null   f1     fx     fy     fxy    fx2     fy2    fx2y   fxy2   fx3    fy3];
        divphi=[null, f1, null, fy, dxfx2, null, dxfx2y, fy2, dxfx3, null, null, null, f1, fx, null, dyfy2, fx2, dyfxy2, null, dyfy3];
        gradphi=Matrix(undef, 2, 40)
        gradphi[1:2, 1: 2]=[null     null;
                            null     null]
        gradphi[1:2, 3: 4]=[f1       null;
                            null     null]
        gradphi[1:2, 5: 6]=[null     f1;
                            null     null]
        gradphi[1:2, 7: 8]=[fy       fx;
                            null     null]
        gradphi[1:2, 9:10]=[dxfx2    null;
                            null     null]
        gradphi[1:2,11:12]=[null     dyfy2;
                            null     null]
        gradphi[1:2,13:14]=[dxfx2y   fx2;
                            null     null]
        gradphi[1:2,15:16]=[fy2      dyfxy2;
                            null     null]
        gradphi[1:2,17:18]=[dxfx3    null;
                            null     null]
        gradphi[1:2,19:20]=[null     dyfy3;
                            null     null]
        gradphi[1:2,21:22]=[null     null;
                            null     null]
        gradphi[1:2,23:24]=[null     null;
                            f1       null]
        gradphi[1:2,25:26]=[null     null;
                            null     f1]
        gradphi[1:2,27:28]=[null     null;
                            fy       fx]
        gradphi[1:2,29:30]=[null     null;
                            dxfx2    null]
        gradphi[1:2,31:32]=[null     null;
                            null     dyfy2]
        gradphi[1:2,33:34]=[null     null;
                            dxfx2y   fx2]
        gradphi[1:2,35:36]=[null     null;
                            fy2      dyfxy2]
        gradphi[1:2,37:38]=[null     null;
                            dxfx3    null]
        gradphi[1:2,39:40]=[null     null;
                            null     dyfy3]

        nFace=20;
        nEdge=0;
        nVert=0;

        cm=Dict([1,2]=>[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [2,3]=>[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [3,4]=>[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [1,4]=>[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]);
    elseif type==:VecP2

        phi=[h21_21 null   h21_20 null   h22_21 null   h21_22 null   h20_21 null   h20_20 null  h22_20 null   h22_22 null   h20_22 null  ;
             null   h21_21 null   h21_20 null   h22_21 null   h21_22 null   h20_21 null  h20_20 null   h22_20 null   h22_22 null   h20_22]
        divphi=[Dxh21_21,Dyh21_21,Dxh21_20,Dyh21_20,Dxh22_21,Dyh22_21,Dxh21_22,Dyh21_22,Dxh20_21,Dyh20_21,
               Dxh20_20,Dyh20_20,Dxh22_20,Dyh22_20,Dxh22_22,Dyh22_22,Dxh20_22,Dyh20_22];
        gradphi=Matrix(undef, 2, 36)
        gradphi[1:2, 1: 2]=[Dxh21_21 Dyh21_21;
                           null     null]
        gradphi[1:2, 3: 4]=[null     null
                           Dxh21_21 Dyh21_21]
        gradphi[1:2, 5: 6]=[Dxh21_20 Dyh21_20;
                           null     null]
        gradphi[1:2, 7: 8]=[null     null
                           Dxh21_20 Dyh21_20]
        gradphi[1:2, 9:10]=[Dxh22_21 Dyh22_21;
                           null     null]
        gradphi[1:2,11:12]=[null     null
                           Dxh22_21 Dyh22_21]
        gradphi[1:2,13:14]=[Dxh21_22 Dyh21_22;
                           null     null]
        gradphi[1:2,15:16]=[null     null
                           Dxh21_22 Dyh21_22]
        gradphi[1:2,17:18]=[Dxh20_21 Dyh20_21;
                           null     null]
        gradphi[1:2,19:20]=[null     null
                           Dxh20_21 Dyh20_21]
        gradphi[1:2,21:22]=[Dxh20_20 Dyh20_20;
                           null     null]
        gradphi[1:2,23:24]=[null     null
                           Dxh20_20 Dyh20_20]
        gradphi[1:2,25:26]=[Dxh22_20 Dyh22_20;
                           null     null]
        gradphi[1:2,27:28]=[null     null
                           Dxh22_20 Dyh22_20]
        gradphi[1:2,29:30]=[Dxh22_22 Dyh22_22;
                           null     null]
        gradphi[1:2,31:32]=[null     null
                           Dxh22_22 Dyh22_22]
        gradphi[1:2,33:34]=[Dxh20_22 Dyh20_22;
                           null     null]
        gradphi[1:2,35:36]=[null     null
                           Dxh20_22 Dyh20_22]

        nFace=2;
        nEdge=2;
        nVert=2;

        cm=Dict([1,2]=>[0,0,1,1,0,0,0,0,0,0,1,1,1,1,0,0,0,0], [2,3]=>[0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1,0,0], [3,4]=>[0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1], [1,4]=>[0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1]);

    elseif type==:VecDG2

        phi=[h21_21 null   h21_20 null   h22_21 null   h21_22 null   h20_21 null   h20_20 null  h22_20 null   h22_22 null   h20_22 null  ;
             null   h21_21 null   h21_20 null   h22_21 null   h21_22 null   h20_21 null  h20_20 null   h22_20 null   h22_22 null   h20_22]
        divphi=[Dxh21_21,Dyh21_21,Dxh21_20,Dyh21_20,Dxh22_21,Dyh22_21,Dxh21_22,Dyh21_22,Dxh20_21,Dyh20_21,
               Dxh20_20,Dyh20_20,Dxh22_20,Dyh22_20,Dxh22_22,Dyh22_22,Dxh20_22,Dyh20_22];
        gradphi=Matrix(undef, 2, 36)
        gradphi[1:2, 1: 2]=[Dxh21_21 Dyh21_21;
                           null     null]
        gradphi[1:2, 3: 4]=[null     null
                           Dxh21_21 Dyh21_21]
        gradphi[1:2, 5: 6]=[Dxh21_20 Dyh21_20;
                           null     null]
        gradphi[1:2, 7: 8]=[null     null
                           Dxh21_20 Dyh21_20]
        gradphi[1:2, 9:10]=[Dxh22_21 Dyh22_21;
                           null     null]
        gradphi[1:2,11:12]=[null     null
                           Dxh22_21 Dyh22_21]
        gradphi[1:2,13:14]=[Dxh21_22 Dyh21_22;
                           null     null]
        gradphi[1:2,15:16]=[null     null
                           Dxh21_22 Dyh21_22]
        gradphi[1:2,17:18]=[Dxh20_21 Dyh20_21;
                           null     null]
        gradphi[1:2,19:20]=[null     null
                           Dxh20_21 Dyh20_21]
        gradphi[1:2,21:22]=[Dxh20_20 Dyh20_20;
                           null     null]
        gradphi[1:2,23:24]=[null     null
                           Dxh20_20 Dyh20_20]
        gradphi[1:2,25:26]=[Dxh22_20 Dyh22_20;
                           null     null]
        gradphi[1:2,27:28]=[null     null
                           Dxh22_20 Dyh22_20]
        gradphi[1:2,29:30]=[Dxh22_22 Dyh22_22;
                           null     null]
        gradphi[1:2,31:32]=[null     null
                           Dxh22_22 Dyh22_22]
        gradphi[1:2,33:34]=[Dxh20_22 Dyh20_22;
                           null     null]
        gradphi[1:2,35:36]=[null     null
                           Dxh20_22 Dyh20_22]

        nFace=18;
        nEdge=0;
        nVert=0;

        cm=Dict([1,2]=>[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [2,3]=>[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [3,4]=>[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [1,4]=>[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]);

    else
        error("Unzulässiger finite-Elemente-Raum");
    end
    return phi, divphi, gradphi, cm, nFace, nEdge, nVert;
end
