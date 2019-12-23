function getQuadElementProperties(type::Symbol, kubPoints::Array{Float64,2})
    null(x,y)=0.0;
    f1(x,y)=1.0;
    fm1(x,y)=-1.0;
    fx(x,y)=x;
    fmx(x,y)=-x;
    fy(x,y)=y;
    fmy(x,y)=-y;
    f1mx(x,y)=1-x;
    fxm1(x,y)=x-1;
    f1my(x,y)=1-y;
    fym1(x,y)=y-1;
    fymxy(x,y)=y-x*y;
    fxmxy(x,y)=x-x*y;
    fxymy(x,y)=x*y-y;
    fxymx(x,y)=x*y-x;
    fxy(x,y)=x*y;
    fl(x,y)=1-x-y+x*y;
    fml(x,y)=-(1-x-y+x*y);

    f1P2(x,y)=((2*x-1)*(x-1))*((2*y-1)*(y-1));
    f2P2(x,y)=(4*(x-x^2))*((2*y-1)*(y-1));
    f3P2(x,y)=(2*x^2-x)*((2*y-1)*(y-1));
    f4P2(x,y)=((2*x-1)*(x-1))*(4*(y-y^2));
    f5P2(x,y)=(4*(x-x^2))*(4*(y-y^2));
    f6P2(x,y)=(2*x^2-x)*(4*(y-y^2));
    f7P2(x,y)=((2*x-1)*(x-1))*(2*y^2-y);
    f8P2(x,y)=(4*(x-x^2))*(2*y^2-y);
    f9P2(x,y)=(2*x^2-x)*(2*y^2-y);

    fgrad11P2(x,y)=(-3+4*x)*(1-3*y+2*y^2);
    fgrad12P2(x,y)=(1-3*x+2*x^2)*(-3+4*y);
    fgrad21P2(x,y)=4*(1-2*x)*(-1+y)*(-1+2*y);
    fgrad22P2(x,y)=-4*(-1+x)*x*(-3+4*y);
    fgrad31P2(x,y)=(-1+4*x)*(-1+y)*(-1+2*y);
    fgrad32P2(x,y)=x*(-1+2*x)*(-3+4*y);
    fgrad41P2(x,y)=-4*(-3+4*x)*(-1+y)*y;
    fgrad42P2(x,y)=4*(-1+x)*(-1+2*x)*(1-2*y);
    fgrad51P2(x,y)=16*(-1+2*x)*(-1+y)*y;
    fgrad52P2(x,y)=16*(-1+x)*x*(-1+2*y);
    fgrad61P2(x,y)=-4*(-1+4*x)*(-1+y)*y;
    fgrad62P2(x,y)=4*x*(-1+2*x)*(1-2*y);
    fgrad71P2(x,y)=(-3+4*x)*y*(-1+2*y);
    fgrad72P2(x,y)=(-1+x)*(-1+2*x)*(-1+4*y);
    fgrad81P2(x,y)=4*(1-2*x)*y*(-1+2*y);
    fgrad82P2(x,y)=-4*(-1+x)*x*(-1+4*y);
    fgrad91P2(x,y)=(-1+4*x)*y*(-1+2*y);
    fgrad92P2(x,y)=x*(-1+2*x)*(-1+4*y);


    f1RT1(x,y)=(1-x)*(2*y-1)*(y-1);
    f2RT1(x,y)=4*(1-y)*(x-x^2);
    f3RT1(x,y)=x*(2*y-1)*(y-1);
    f4RT1(x,y)=(1-y)*(2*x^2-x);
    f5RT1(x,y)=4*x*(y-y^2);
    f6RT1(x,y)=y*(2*x^2-x);
    f7RT1(x,y)=x*(2*y^2-y);
    f8RT1(x,y)=4*y*(x-x^2);
    f9RT1(x,y)=(1-x)*(2*y^2-y);
    f10RT1(x,y)=y*(2*x-1)*(x-1);
    f11RT1(x,y)=4*(1-x)*(y-y^2);
    f12RT1(x,y)=(1-y)*(2*x-1)*(x-1);

    fdiv1RT1(x,y)=(1-x)*(4*y-3);
    fdiv2RT1(x,y)=4*(1-y)*(1-2*x);
    fdiv3RT1(x,y)=x*(4*y-3);
    fdiv4RT1(x,y)=(1-y)*(4*x-1);
    fdiv5RT1(x,y)=4*x*(1-2*y);
    fdiv6RT1(x,y)=y*(4*x-1);
    fdiv7RT1(x,y)=x*(4*y-1);
    fdiv8RT1(x,y)=4*y*(1-2*x);
    fdiv9RT1(x,y)=(1-x)*(4*y-1);
    fdiv10RT1(x,y)=y*(4*x-3);
    fdiv11RT1(x,y)=4*(1-x)*(1-2*y);
    fdiv12RT1(x,y)=(1-y)*(4*x-3);


    #Namen: fgrad-indexphi-indexingradphimatrix-RT1
    fgrad121RT1(x,y)=-(2*y-1)*(y-1);
    fgrad122RT1(x,y)=(1-x)*(4*y-3);
    fgrad211RT1(x,y)=4*(1-y)*(1-2*x);
    fgrad212RT1(x,y)=-4*(x-x^2);
    fgrad321RT1(x,y)=(2*y-1)*(y-1);
    fgrad322RT1(x,y)=x*(4*y-3);
    fgrad411RT1(x,y)=(1-y)*(4*x-1);
    fgrad412RT1(x,y)=-(2*x^2-x);
    fgrad521RT1(x,y)=4*(y-y^2);
    fgrad522RT1(x,y)=4*x*(1-2*y);
    fgrad611RT1(x,y)=y*(4*x-1);
    fgrad612RT1(x,y)=2*x^2-x;
    fgrad721RT1(x,y)=2*y^2-y;
    fgrad722RT1(x,y)=x*(4*y-1);
    fgrad811RT1(x,y)=4*y*(1-2*x);
    fgrad812RT1(x,y)=4*(x-x^2);
    fgrad921RT1(x,y)=-(2*y^2-y);
    fgrad922RT1(x,y)=(1-x)*(4*y-1);
    fgrad1011RT1(x,y)=y*(4*x-3);
    fgrad1012RT1(x,y)=(2*x-1)*(x-1);
    fgrad1121RT1(x,y)=-4*(y-y^2);
    fgrad1122RT1(x,y)=4*(1-x)*(1-2*y);
    fgrad1211RT1(x,y)=(1-y)*(4*x-3);
    fgrad1212RT1(x,y)=-(2*x-1)*(x-1);

    if type==:DG0
        phi=[f1];
        #c=reshape([0.5; 0.5],2,1);
        divphi=[null];
        gradphi=reshape([null,null],2,1);

        nFace=1;
        nEdge=0;
        nVert=0;

        #cm=[0.5 0.0 0.0; 1.0 0.5 0.0; 0.5 1.0 0.0; 0.0 0.5 0.0];
        comp=[0];

    elseif type==:P1
        phi=[fl,fxmxy,fxy,fymxy];
        #c=[0.0 1.0 1.0 0.0 ;
        #    0.0 0.0 1.0 1.0];

        divphi=[null null null null];
        gradphi=[fym1 f1my fy fmy;
                 fxm1 fmx fx f1mx];

        nFace=0;
        nEdge=0;
        nVert=1;

        #cm=[0.5 0.0 1.0 1.0 0.0 0.0;
        #   1.0 0.5 0.0 1.0 1.0 0.0;
        #   0.5 1.0 0.0 0.0 1.0 1.0;
        #   0.0 0.5 1.0 0.0 0.0 1.0];
        comp=[0, 0, 0, 0];

    #=
    elseif type==:P1x
        phi=[f1mx,fx];
        #c=[0.0 1.0;0.5 0.5];

        divphi=[null null];
        gradphi=[fm1 f1;
                 null null];

        nFace=0;
        nEdge=1;
        nVert=0;

        #cm=[0.5 0.0 0.0 0.0; 1.0 0.5 0.0 1.0; 0.5 1.0 0.0 0.0; 0.0 0.5 1.0 0.0];
        comp=[0, 0];

    elseif type==:P1y
        phi=[f1my,fy];
        #c=[0.5 0.5;0.0 1.0];

        divphi=[null null];
        gradphi=[null null;
                 fm1 f1];

        nFace=0;
        nEdge=1;
        nVert=0;

        #cm=[0.5 0.0 1.0 0.0; 1.0 0.5 0.0 0.0; 0.5 1.0 0.0 1.0; 0.0 0.5 0.0 0.0];
        comp=[0, 0];
    =#
    elseif type==:DG1
        phi=[fl,fxmxy,fxy,fymxy];
        #c=[0.0 1.0 1.0 0.0 ;
        #   0.0 0.0 1.0 1.0];

        divphi=[null null null null];
        gradphi=[fym1 f1my fy fmy;
                 fxm1 fmx fx f1mx];


        nFace=4;
        nEdge=0;
        nVert=0;

        #cm=[0.5 0.0 0.0 0.0 0.0 0.0;
        #   1.0 0.5 0.0 0.0 0.0 0.0;
        #   0.5 1.0 0.0 0.0 0.0 0.0;
        #   0.0 0.5 0.0 0.0 0.0 0.0];
        comp=[0, 0, 0, 0];
    #=
    elseif type==:DG1x
        phi=[f1mx,fx];
        #c=[0.0 1.0;0.5 0.5];

        divphi=[null null];
        gradphi=[fm1 f1;
                 null null];

        nFace=2;
        nEdge=0;
        nVert=0;

        #cm=[0.5 0.0 0.0 0.0; 1.0 0.5 0.0 0.0; 0.5 1.0 0.0 0.0; 0.0 0.5 0.0 0.0];
        comp=[0, 0];

    elseif type==:DG1y
        phi=[f1mx,fx];
        #c=[0.0 1.0;0.5 0.5];

        divphi=[null null];
        gradphi=[fm1 f1;
                 null null];

        nFace=2;
        nEdge=0;
        nVert=0;

        #cm=[0.5 0.0 0.0 0.0; 1.0 0.5 0.0 0.0; 0.5 1.0 0.0 0.0; 0.0 0.5 0.0 0.0];
        comp=[0, 0];
    =#

    elseif type==:P2
        #phi=[f1P2, f2P2, f3P2, f4P2, f5P2, f6P2, f7P2, f8P2, f9P2];
        phi=[f5P2, f2P2, f6P2, f8P2, f4P2, f1P2, f3P2, f9P2, f7P2];

        #c=[0.5 0.5 1.0 0.5 0.0 0.0 1.0 1.0 0.0;
        #   0.5 0.0 0.5 1.0 0.5 0.0 0.0 1.0 1.0];

        divphi=[null null null null null null null null null];
        gradphi=[fgrad51P2 fgrad21P2 fgrad61P2 fgrad81P2 fgrad41P2 fgrad11P2 fgrad31P2 fgrad91P2 fgrad71P2;
                 fgrad52P2 fgrad22P2 fgrad62P2 fgrad82P2 fgrad42P2 fgrad12P2 fgrad32P2 fgrad92P2 fgrad72P2];

        nFace=1;
        nEdge=1;
        nVert=1;

        comp=[0, 0, 0, 0, 0, 0, 0, 0, 0];

    elseif type==:DG2
        phi=[f5P2, f2P2, f6P2, f8P2, f4P2, f1P2, f3P2, f9P2, f7P2];
        #c=[0.5 0.5 1.0 0.5 0.0 0.0 1.0 1.0 0.0;
        #   0.5 0.0 0.5 1.0 0.5 0.0 0.0 1.0 1.0];

        divphi=[null null null null null null null null null];
        gradphi=[fgrad51P2 fgrad21P2 fgrad61P2 fgrad81P2 fgrad41P2 fgrad11P2 fgrad31P2 fgrad91P2 fgrad71P2;
                 fgrad52P2 fgrad22P2 fgrad62P2 fgrad82P2 fgrad42P2 fgrad12P2 fgrad32P2 fgrad92P2 fgrad72P2];

        comp=[0, 0, 0, 0, 0, 0, 0, 0, 0];

    elseif type==:RT0

        phi=[null fx null f1mx;
             f1my null fy null];
        #c=[0.5 1.0 0.5 0.0;
        #   0.0 0.5 1.0 0.5];

        divphi=[fm1, f1, f1, fm1];
        gradphi=[null null f1 null null null fm1 null;
                 null fm1 null null null f1 null null];

        nFace=0;
        nEdge=1;
        nVert=0;

        #cm=[0.5 0.0 1.0 0.0 0.0 0.0; 1.0 0.5 0.0 1.0 0.0 0.0; 0.5 1.0 0.0 0.0 1.0 0.0; 0.0 0.5 0.0 0.0 0.0 1.0];
        comp=[2,1,2,1];

    elseif type==:RT0B #Broken RT0

        phi=[null fx null f1mx; f1my null fy null];
        #c=[0.5 1.0 0.5 0.0; 0.0 0.5 1.0 0.5];

        divphi=[fm1, f1, f1, fm1];
        gradphi=[null null f1 null null null fm1 null;
                 null fm1 null null null f1 null null];

        nFace=4;
        nEdge=0;
        nVert=0;

        #cm=[0.5 0.0 0.0 0.0 0.0 0.0; 1.0 0.5 0.0 0.0 0.0 0.0; 0.5 1.0 0.0 0.0 0.0 0.0; 0.0 0.5 0.0 0.0 0.0 0.0];
        comp=[2,1,2,1];



    elseif type==:RT1

        phi=[f2RT1 null f8RT1 null null null f4RT1 f6RT1 null null f12RT1 f10RT1;
             null f5RT1 null f11RT1 f1RT1 f3RT1 null null f9RT1 f7RT1 null null];
         #c=[0.5 1.0 0.5 0.0 0.0 1.0 1.0 1.0 1.0 0.0 0.0 0.0;
         #   0.0 0.5 1.0 0.5 0.0 0.0 0.0 1.0 1.0 1.0 1.0 0.0];

        divphi=[fdiv2RT1, fdiv5RT1, fdiv8RT1, fdiv11RT1, fdiv1RT1, fdiv3RT1, fdiv4RT1, fdiv6RT1, fdiv9RT1, fdiv7RT1, fdiv12RT1, fdiv10RT1];

        gradphi=[fgrad211RT1 fgrad212RT1 null null fgrad811RT1 fgrad812RT1 null null null null null null fgrad411RT1 fgrad412RT1 fgrad611RT1 fgrad612RT1 null null null null fgrad1211RT1 fgrad1212RT1 fgrad1011RT1 fgrad1012RT1;
                 null null fgrad521RT1 fgrad522RT1 null null fgrad1121RT1 fgrad1122RT1 fgrad121RT1 fgrad122RT1 fgrad321RT1 fgrad322RT1 null null null null fgrad921RT1 fgrad922RT1 fgrad721RT1 fgrad722RT1 null null null null];


        nFace=4;
        nEdge=2;
        nVert=0;

        comp=[1,2,1,2,2,2,1,1,2,2,1,1];

    elseif type==:RT1B

        phi=[f2RT1 null f8RT1 null null null f4RT1 f6RT1 null null f10RT1 f12RT1;
             null f5RT1 null f11RT1 f1RT1 f3RT1 null null f7RT1 f9RT1 null null];
         #c=[0.5 1.0 0.5 0.0 0.0 1.0 1.0 1.0 1.0 0.0 0.0 0.0;
         #   0.0 0.5 1.0 0.5 0.0 0.0 0.0 1.0 1.0 1.0 1.0 0.0];

        divphi=[fdiv2RT1, fdiv5RT1, fdiv8RT1, fdiv11RT1, fdiv1RT1, fdiv3RT1, fdiv4RT1, fdiv6RT1, fdiv7RT1, fdiv9RT1, fdiv10RT1, fdiv12RT1];

        gradphi=[fgrad211RT1 fgrad212RT1 null null fgrad811RT1 fgrad812RT1 null null null null null null fgrad411RT1 fgrad412RT1 fgrad611RT1 fgrad612RT1 null null null null fgrad1011RT1 fgrad1012RT1 fgrad1211RT1 fgrad1212RT1;
                 null null fgrad521RT1 fgrad522RT1 null null fgrad1121RT1 fgrad1122RT1 fgrad121RT1 fgrad122RT1 fgrad321RT1 fgrad322RT1 null null null null fgrad721RT1 fgrad722RT1 fgrad921RT1 fgrad922RT1 null null null null];


        nFace=12;
        nEdge=0;
        nVert=0;

        comp=[1,2,1,2,1,2,1,2,1,2,1,2];

    elseif type==:VecP1

        phi=[fl null fxmxy null fxy null fymxy null;
             null fl null fxmxy null fxy null fymxy];

        divphi=[fym1, fxm1, f1my, fmx, fy, fx, fmy, f1mx];
        gradphi=[fym1 fxm1 null null f1my fmx null null fy fx null null fmy f1mx null null ;
                null null fym1 fxm1 null null f1my fmx null null fy fx null null fmy f1mx ];

        #c=[0.0 0.0 1.0 1.0 1.0 1.0 0.0 0.0;
        #   0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0];

        #cm=[0.5 0.0 0.0 1.0 0.0 1.0 0.0 0.0 0.0 0.0;
        #    1.0 0.5 0.0 0.0 1.0 0.0 1.0 0.0 0.0 0.0;
        #    0.5 1.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 1.0;
        #    0.0 0.5 1.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0];

        nFace=0;
        nEdge=0;
        nVert=2;

        comp=[1, 2, 1, 2, 1, 2, 1, 2];

    elseif type==:VecDG1

        phi=[fl null fxmxy null fxy null fymxy null;
             null fl null fxmxy null fxy null fymxy];

        divphi=[fym1, fxm1, f1my, fmx, fy, fx, fmy, f1mx];
        gradphi=[fym1 fxm1 null null f1my fmx null null fy fx null null fmy f1mx null null ;
                null null fym1 fxm1 null null f1my fmx null null fy fx null null fmy f1mx ];


        #c=[0.0 0.0 1.0 1.0 1.0 1.0 0.0 0.0;
        #   0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0];

        #cm=[0.5 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
        #    1.0 0.5 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
        #    0.5 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
        #    0.0 0.5 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0];

        comp=[1, 2, 1, 2, 1, 2, 1, 2];

        nFace=8;
        nEdge=0;
        nVert=0;

    else
        error("Unzulässiger finite-Elemente-Raum");
    end

    sk=size(kubPoints,2);
    kubPhi=Array{Array{Float64,2},ndims(phi)}(undef,size(phi));
    for k=1:length(phi)
        kubVal=Array{Float64,2}(undef,sk,sk);
        for i=1:sk, j=1:sk
            kubVal[i,j]=phi[k](kubPoints[1,i], kubPoints[2,j]);
        end
        kubPhi[k]=kubVal;
    end

    kubDiv=Array{Array{Float64,2},1}(undef,length(divphi));
    for k in 1:length(divphi)
        kubVal=Array{Float64,2}(undef,sk,sk);
        for i=1:sk, j=1:sk
            kubVal[i,j]=divphi[k](kubPoints[1,i], kubPoints[2,j]);
        end
        kubDiv[k]=kubVal;
    end

    kubGrad=Array{Array{Float64,2},2}(undef,size(gradphi));
    for ki=1:size(gradphi,1), kj=1:size(gradphi,2)
        kubVal=Array{Float64,2}(undef,sk,sk);
        for i=1:sk, j=1:sk
            kubVal[i,j]=gradphi[ki,kj](kubPoints[1,i], kubPoints[2,j]);
        end
        kubGrad[ki,kj]=kubVal;
    end

    return kubPhi, kubDiv,  kubGrad, comp, nFace, nEdge, nVert
end


function getTriElementProperties(type::Symbol, kubPoints::Array{Float64,2})
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

        #cm=[0.5 0.0 0.0; 0.5 0.5 0.0; 0.0 0.5 0.0];
        comp=[0];

    elseif type==:RT0

        phi=[fxm1 fx fsx; fy fym1 fsy];
        #c=[0.5 0.0 0.5; 0.0 0.5 0.5];

        divphi=[f2, f2, f2s];
        gradphi=[f1 null f1 null f1s null; null f1 null f1 null f1s];

        nFace=0;
        nEdge=1;
        nVert=0;

        #cm=[0.5 0.0 0.0 1.0 0.0;
        #    0.5 0.5 0.0 0.0 1.0;
        #    0.0 0.5 1.0 0.0 0.0];
        comp=[0];


    elseif type==:P1
        phi=[f1mxy fx fy];
        #c=[0.0 1.0 0.0; 0.0 0.0 1.0];

        divphi=[null null null];

        gradphi=[fmy f1 null;
                 fmx null f1];

        nFace=0;
        nEdge=0;
        nVert=1;

        #cm=[0.5 0.0 1.0 1.0 0.0;
        #    0.5 0.5 0.0 1.0 1.0;
        #    0.0 0.5 1.0 0.0 1.0];

        comp=[0, 0, 0];
    end

    sk=size(kubPoints,2);
    kubPhi=Array{Array{Float64,2},ndims(phi)}(undef,size(phi));
    for k=1:length(phi)
        kubVal=Array{Float64,2}(undef,sk,sk);
        for i=1:sk, j=1:sk
            kubVal[i,j]=phi[k](kubPoints[1,i], kubPoints[2,j]);
        end
        kubPhi[k]=kubVal;
    end

    kubDiv=Array{Array{Float64,2},1}(undef,length(divphi));
    for k in 1:length(divphi)
        kubVal=Array{Float64,2}(undef,sk,sk);
        for i=1:sk, j=1:sk
            kubVal[i,j]=divphi[k](kubPoints[1,i], kubPoints[2,j]);
        end
        kubDiv[k]=kubVal;
    end

    kubGrad=Array{Array{Float64,2},2}(undef,size(gradphi));
    for ki=1:size(gradphi,1), kj=1:size(gradphi,2)
        kubVal=Array{Float64,2}(undef,sk,sk);
        for i=1:sk, j=1:sk
            kubVal[i,j]=gradphi[ki,kj](kubPoints[1,i], kubPoints[2,j]);
        end
        kubGrad[ki,kj]=kubVal;
    end

    return kubPhi, kubDiv,  comp, nFace, nEdge, nVert
end
