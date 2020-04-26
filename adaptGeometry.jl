function adaptGeometry!(m::mesh,hm::Float64,a::Float64)
    H=m.geometry.r[2];
    xM=0.5(m.geometry.l[1]+m.geometry.r[1])
    h(x)=(hm*a^2)/((x-xM)^2+a^2); #Witch of Agnesi
    coord=m.geometry.coordinates;
    for i in 1:size(coord,2)
        z0=h(coord[1,i])
        coord[2,i]=H*(coord[2,i]+z0)/(H+z0); #Gal-Chen & Sommerville transformation
    end
    return nothing;
end

function adaptGeometry!(m::mesh,hm::Float64,a::Float64,lambda::Float64)
    H=m.geometry.r[2];
    xM=0.5(m.geometry.l[1]+m.geometry.r[1])
    function h(x) #Schaer
        h2=hm*exp(-((x-xM)/a)^2);
        h1=h2*cos(pi*(x-xM)/lambda)^2;
        return h1, h2;
    end
    coord=m.geometry.coordinates;
    sC=4000.0
    sF=1000.0
    nC=1.35
    nF=1.35
    for i in 1:size(coord,2)
        zSF, zSC=h(coord[1,i])
        coord[2,i]+=sinh((H/sC)^nC-(coord[2,i]/sC)^nC)/sinh((H/sC)^nC)*zSC+sinh((H/sF)^nF-(coord[2,i]/sF)^nF)/sinh((H/sF)^nF)*(zSF-zSC); #SLEVE
    end
    return nothing;
end

function adaptGeometry!(m::mesh,vid::Array{Int64,1},c::Array{Float64,2})
    m.geometry.coordinates[:,vid]=c;
    return nothing;
end

function adaptGeometry!(m::mesh,vid::Int64,c::Array{Float64,1})
    m.geometry.coordinates[:,vid]=c;
    return nothing;
end

function adaptGeometry!(m::mesh,r::Float64,adaptBoundary::Bool=false)
    coord=m.geometry.coordinates;
    xR=m.geometry.r[1]; xL=m.geometry.l[1]; yR=m.geometry.r[2]; yL=m.geometry.l[2]
    if adaptBoundary
        for i in 1:size(coord,2)
            x=r*rand([-1,1])*rand(1)[1];
            y=rand([-1,1])*sqrt(r*r-x*x)*rand(1)[1];
            coord[:,i]=[coord[1,i]+x, coord[2,i]+y];
        end
    else
        for i in 1:size(coord,2)
            if coord[1,i]!=xL && coord[1,i]!=xR && coord[2,i]!=yL && coord[2,i]!=yR
                x=r*rand([-1,1])*rand(1)[1];
                y=rand([-1,1])*sqrt(r*r-x*x)*rand(1)[1];
                coord[:,i]=[coord[1,i]+x, coord[2,i]+y];
            end
        end
    end

    return nothing;
end

function adaptGeometry!(m::mesh,pert::Tuple{Float64,Float64},adaptBoundary::Bool=false)
    coord=m.geometry.coordinates;
    xR=m.geometry.r[1]; xL=m.geometry.l[1]; yR=m.geometry.r[2]; yL=m.geometry.l[2]
    dx=(xR-xL)/m.topology.n[1];
    dy=(yR-yL)/m.topology.n[2];

    if adaptBoundary
        for i in 1:size(coord,2)
            x=(coord[1,i]);
            y=(coord[2,i]);
            hx=coord[1,i]+pert[1]*sin(2*pi*(y-yL)/(yR-yL))*dx;
            hy=coord[2,i]+pert[2]*sin(2*pi*(x-xL)/(xR-xL))*dy;
            coord[:,i]=[hx, hy];
        end
    else
        for i in 1:size(coord,2)
            if coord[1,i]!=xL && coord[1,i]!=xR && coord[2,i]!=yL && coord[2,i]!=yR
                x=(coord[1,i]);
                y=(coord[2,i]);
                hx=coord[1,i]+pert[1]*sin(2*pi*(y-yL)/(yR-yL))*dx;
                hy=coord[2,i]+pert[2]*sin(2*pi*(x-xL)/(xR-xL))*dy;
                coord[:,i]=[hx, hy];
            end
        end
    end

    return nothing;
end
