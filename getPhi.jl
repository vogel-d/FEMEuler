function getPhi(type::Symbol)
    if type==:DG0
        phi1(x,y)=1.0;
        phi=[phi1];
        psize=(1,1);

    elseif type==:P1
        phip1(x,y)=(1-x-y+x*y); #(0,0)
        phip2(x,y)=(x-x*y); #(1,0)
        phip3(x,y)=(x*y); #(1,1)
        phip4(x,y)=(y-x*y); #(0,1)

        phi=[phip1,phip2,phip3,phip4];
        psize=(1,4);
    #=
    elseif type==:P1x || type==:DG1x
        phipx1(x,y)=(1-x);
        phipx2(x,y)=(x);
        phi=[phipx1,phipx2];
        psize=(1,2);

    elseif type==:P1y || type==:DG1y
        phipy1(x,y)=(1-y);
        phipy2(x,y)=(y);
        phi=[phipy1,phipy2];
        psize=(1,2);
    =#
    elseif type==:DG1
        phid1(x,y)=(1-x-y+x*y); #(0,0)
        phid2(x,y)=(x-x*y); #(1,0)
        phid3(x,y)=(x*y); #(1,1)
        phid4(x,y)=(y-x*y); #(0,1)

        phi=[phid1,phid2,phid3,phid4];
        psize=(1,4);

    elseif type==:RT0
        phiv1(x,y)=x;
        phiv0(x,y)=0;
        phiv2(x,y)=(1-x);
        phiv3(x,y)=y;
        phiv4(x,y)=(1-y);

        phi=[phiv0 phiv1 phiv0 phiv2;
            phiv4 phiv0 phiv3 phiv0];
        #phi=[phiv0 phiv4; phiv1 phiv0; phiv0 phiv3; phiv2 phiv0];
        psize=(2,4);

    elseif type==:RT0B
        phiv1d(x,y)=x;
        phiv0d(x,y)=0;
        phiv2d(x,y)=(1-x);
        phiv3d(x,y)=y;
        phiv4d(x,y)=(1-y);

        phi=[phiv0d phiv1d phiv0d phiv2d;
             phiv4d phiv0d phiv3d phiv0d];
        #phi=[phiv0d phi4d; phiv1d phiv0d; phiv0d phiv3d; phiv2d phiv0d]';
        psize=(2,4);

    elseif type==:VecP1
        phiw0(x,y)=0;
        phiw1(x,y)=(1-x)*(1-y);
        phiw2(x,y)=x*(1-y);
        phiw3(x,y)=x*y;
        phiw4(x,y)=(1-x)*y;

        phi=[phiw0 phiw0 phiw2 phiw3 phiw0 phiw0 phiw4 phiw1;
             phiw1 phiw2 phiw0 phiw0 phiw3 phiw4 phiw0 phiw0];
        #phi=[phiw0 phiw1; phiw0 phiw2; phiw2 phiw0; phiw3 phiw0; phiw0 phiw3; phiw0  phiw4; phiw4 phiw0;phiw1 phiw0]

        psize=(2,8);

    elseif type==:VecDG1
        phiw0d(x,y)=0;
        phiw1d(x,y)=(1-x)*(1-y);
        phiw2d(x,y)=x*(1-y);
        phiw3d(x,y)=x*y;
        phiw4d(x,y)=(1-x)*y;

        phi=[phiw0d phiw0d phiw2d phiw3d phiw0d phiw0d phiw4d phiw1d;
             phiw1d phiw2d phiw0d phiw0d phiw3d phiw4d phiw0d phiw0d];
        #phi=[phiw0d phiw1d; phiw0d phiw2d; phiw2d phiw0d; phiw3d phiw0d; phiw0d phiw3d; phiw0d  phiw4d; phiw4d phiw0d;phiw1d phiw0d]

        psize=(2,8);
    end

    return phi, psize;
end
