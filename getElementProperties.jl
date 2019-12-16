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
    fxy(x,y)=x*y;
    fl(x,y)=1-x-y+x*y;

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
