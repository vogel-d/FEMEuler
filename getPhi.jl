function getPhi(type::Symbol)

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


    if type==:DG0
        phi=[f1];
        psize=(1,1);

    elseif type==:P1 || type==:DG1
        phi=[fl,fxmxy,fxy,fymxy];
        psize=(1,4);
    #=
    elseif type==:P1x || type==:DG1x
        phi=[f1mx,fx];
        psize=(1,2);

    elseif type==:P1y || type==:DG1y
        phi=[f1my,fy];
        psize=(1,2);
    =#
    elseif type==:P2 || type==:DG2
        phi=[f1P2, f2P2, f3P2, f4P2, f5P2, f6P2, f7P2, f8P2, f9P2];
        psize=(1,9)

    elseif type==:RT0 || type==:RT0B
        phi=[null fx null f1mx;
             f1my null fy null];

        psize=(2,4);

    elseif type==:RT1 || type==:RT1B
        phi=[null f2RT1 null f4RT1 null f6RT1 null f8RT1 null f10RT1 null f12RT1;
             f1RT1 null f3RT1 null f5RT1 null f7RT1 null f9RT1 null f11RT1 null];
        psize=(2,12);

    elseif type==:VecP1 || type==:VecDG1
        phi=[fl null fxmxy null fxy null fymxy null;
             null fl null fxmxy null fxy null fymxy];

        psize=(2,8);
    else
        error("Ansatzfunktionen zu finitem Element $type nicht gefunden.")
    end

    return phi, psize;
end
