#mt=3 fehlt
#=
    einige jacobi-Routinen um je nach Argumenten die gewünschten transformierten Größen zurückzugeben
    bspw. Jacobi-Matrix, transformierte Ansatzfunktionen, Determinante der Jacobi-Matrix, ...
=#

function jacobi!(J::Array{Array{Float64,2},2},dJ::Array{Float64,2},m::mesh, fid::Int64, kubPoints::Array{Float64,2}, coord::Array{Float64,2})
    key="20";
    mt=m.meshType;
    rstart=m.topology.offset[key][fid];
    z=1;
    for j in rstart:rstart+mt-1
        for i in 1:3
            coord[i,z]=m.geometry.coordinates[i,m.topology.incidence[key][j]];
        end
        z+=1;
    end
    sk=size(kubPoints,2);
    sign=m.orientation[fid];
    if mt==3
        a=coord[:,2]-coord[:,1];
        b=coord[:,3]-coord[:,1];
        for i in 1:sk
            dJ[i]=abs(a[1]*b[2]-b[1]*a[2]);

            J[1,1][i]=a[1];
            J[1,2][i]=b[1];
            J[2,1][i]=a[2];
            J[2,2][i]=b[2];
        end
    elseif mt==4
        r=m.geometry.r[1];
         for i=1:sk, j=1:sk
             x=kubPoints[1,i]
             y=kubPoints[2,j]
             f1=(coord[1,1]+(coord[1,2]-coord[1,1])*x+(coord[1,4]-coord[1,1])*y+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*x*y)^2;
             f2=(coord[2,1]+(coord[2,2]-coord[2,1])*x+(coord[2,4]-coord[2,1])*y+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*x*y)^2;
             f3=(coord[3,1]+(coord[3,2]-coord[3,1])*x+(coord[3,4]-coord[3,1])*y+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*x*y)^2;
             f=r*(f1+f2+f3)^(-3/2)
             J[1,1][i,j]=f*(f2+f3)*((coord[1,2]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*y);
             J[2,1][i,j]=f*(f1+f3)*((coord[2,2]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*y);
             J[3,1][i,j]=f*(f1+f2)*((coord[3,2]-coord[3,1])+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*y);
             J[1,2][i,j]=f*(f2+f3)*((coord[1,4]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*x);
             J[2,2][i,j]=f*(f1+f3)*((coord[2,4]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*x);
             J[3,2][i,j]=f*(f1+f2)*((coord[3,4]-coord[3,1])+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*x);
             dJ[i,j]=sign*sqrt((J[2,1][i,j]*J[3,2][i,j]-J[3,1][i,j]*J[2,2][i,j])^2+(J[3,1][i,j]*J[1,2][i,j]-J[1,1][i,j]*J[3,2][i,j])^2+(J[1,1][i,j]*J[2,2][i,j]-J[2,1][i,j]*J[1,2][i,j])^2);
         end
    end
    return nothing;
end

function jacobi!(J::Array{Float64,2},m::mesh, fid::Int64, x::Float64, y::Float64, coord::Array{Float64,2})
    key="20";
    mt=m.meshType;
    rstart=m.topology.offset[key][fid];
    z=1;
    for j in rstart:rstart+mt-1
        for i in 1:3
            coord[i,z]=m.geometry.coordinates[i,m.topology.incidence[key][j]];
        end
        z+=1;
    end
    sign=m.orientation[fid];
    if mt==3
        a=coord[:,2]-coord[:,1];
        b=coord[:,3]-coord[:,1];
        dJ=abs(a[1]*b[2]-b[1]*a[2]);
        J[1,1]=a[1];
        J[1,2]=b[1];
        J[2,1]=a[2];
        J[2,2]=b[2];
    elseif mt==4
        f1=(coord[1,1]+(coord[1,2]-coord[1,1])*x+(coord[1,4]-coord[1,1])*y+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*x*y)^2;
        f2=(coord[2,1]+(coord[2,2]-coord[2,1])*x+(coord[2,4]-coord[2,1])*y+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*x*y)^2;
        f3=(coord[3,1]+(coord[3,2]-coord[3,1])*x+(coord[3,4]-coord[3,1])*y+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*x*y)^2;
        f=m.geometry.r[1];*(f1+f2+f3)^(-3/2)
        J[1,1]=f*(f2+f3)*((coord[1,2]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*y);
        J[2,1]=f*(f1+f3)*((coord[2,2]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*y);
        J[3,1]=f*(f1+f2)*((coord[3,2]-coord[3,1])+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*y);
        J[1,2]=f*(f2+f3)*((coord[1,4]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*x);
        J[2,2]=f*(f1+f3)*((coord[2,4]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*x);
        J[3,2]=f*(f1+f2)*((coord[3,4]-coord[3,1])+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*x);
        dJ=sign*norm(cross(J[:,1],J[:,2]),2)
    end
    return dJ;
end

function jacobi!(J::Array{Array{Float64,2},2},ddJ::Array{Float64,2},jphi::Array{Array{Float64,2},2},m::mesh, fid::Int64, kubPoints::Array{Float64,2}, phi::Array{Array{Float64,2},2}, coord::Array{Float64,2})
    key="20";
    mt=m.meshType;
    rstart=m.topology.offset[key][fid];
    z=1;
    for j in rstart:rstart+mt-1
        for i in 1:3
            coord[i,z]=m.geometry.coordinates[i,m.topology.incidence[key][j]];
        end
        z+=1;
    end

    sk=size(kubPoints,2);
    sign=m.orientation[fid];
    if mt==3
        a=coord[:,2]-coord[:,1];
        b=coord[:,3]-coord[:,1];
        for i in 1:sk
            ddJ[i]=1/abs(a[1]*b[2]-b[1]*a[2]);

            J[1,1][1,i]=a[1];
            J[1,2][1,i]=b[1];
            J[2,1][1,i]=a[2];
            J[2,2][1,i]=b[2];

            for k in 1:size(phi,2)
                jphi[1,k][1,i]=(J[1,1][1,i]*phi[1,k][1,i]+J[1,2][1,i]*phi[2,k][1,i]);
                jphi[2,k][1,i]=(J[2,1][1,i]*phi[1,k][1,i]+J[2,2][1,i]*phi[2,k][1,i]);
            end
        end
        println("coord: $coord")
        println([a[1] b[1];a[2] b[2]])
    elseif mt==4
        r=m.geometry.r[1];
        for i=1:sk, j=1:sk
            x=kubPoints[1,i]
            y=kubPoints[2,j]
            f1=(coord[1,1]+(coord[1,2]-coord[1,1])*x+(coord[1,4]-coord[1,1])*y+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*x*y)^2;
            f2=(coord[2,1]+(coord[2,2]-coord[2,1])*x+(coord[2,4]-coord[2,1])*y+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*x*y)^2;
            f3=(coord[3,1]+(coord[3,2]-coord[3,1])*x+(coord[3,4]-coord[3,1])*y+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*x*y)^2;
            f=r*(f1+f2+f3)^(-3/2)
            J[1,1][i,j]=f*(f2+f3)*((coord[1,2]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*y);
            J[2,1][i,j]=f*(f1+f3)*((coord[2,2]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*y);
            J[3,1][i,j]=f*(f1+f2)*((coord[3,2]-coord[3,1])+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*y);
            J[1,2][i,j]=f*(f2+f3)*((coord[1,4]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*x);
            J[2,2][i,j]=f*(f1+f3)*((coord[2,4]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*x);
            J[3,2][i,j]=f*(f1+f2)*((coord[3,4]-coord[3,1])+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*x);
            ddJ[i,j]=sign*1/sqrt((J[2,1][i,j]*J[3,2][i,j]-J[3,1][i,j]*J[2,2][i,j])^2+(J[3,1][i,j]*J[1,2][i,j]-J[1,1][i,j]*J[3,2][i,j])^2+(J[1,1][i,j]*J[2,2][i,j]-J[2,1][i,j]*J[1,2][i,j])^2);
            for k in 1:size(phi,2)
                jphi[1,k][i,j]=(J[1,1][i,j]*phi[1,k][i,j]+J[1,2][i,j]*phi[2,k][i,j]);
                jphi[2,k][i,j]=(J[2,1][i,j]*phi[1,k][i,j]+J[2,2][i,j]*phi[2,k][i,j]);
                jphi[3,k][i,j]=(J[3,1][i,j]*phi[1,k][i,j]+J[3,2][i,j]*phi[2,k][i,j]);
            end
        end
    end
    return nothing;
end

function jacobi!(ddJ::Array{Float64,2}, jphi::Array{Array{Float64,2},2}, m::mesh, fid::Int64, kubPoints::Array{Float64,2}, phi::Array{Array{Float64,2},2}, coord::Array{Float64,2})
    key="20";
    mt=m.meshType;
    rstart=m.topology.offset[key][fid];

    z=1;
    for j in rstart:rstart+mt-1
        for i in 1:3
            coord[i,z]=m.geometry.coordinates[i,m.topology.incidence[key][j]];
        end
        z+=1;
    end

    sk=size(kubPoints,2);
    sign=m.orientation[fid];
    if mt==3
        a=coord[:,2]-coord[:,1];
        b=coord[:,3]-coord[:,1];
        for i in 1:sk
            ddJ[i]=1/abs(a[1]*b[2]-b[1]*a[2]);
            JtJ11[i]=a[1]^2+a[2]^2;
            JtJ12[i]=a[1]*b[1]+a[2]*b[2];
            JtJ21[i]=a[1]*b[1]+a[2]*b[2];
            JtJ22[i]=b[1]^2+b[2]^2;
            for k in 1:size(phi,2)
                jphi[1,k][i]=(JtJ11*phi[1,k][i]+JtJ12*phi[2,k][i]);
                jphi[2,k][i]=(JtJ21*phi[1,k][i]+JtJ22*phi[2,k][i]);
            end
        end
    elseif mt==4
        r=m.geometry.r[1];
        for i=1:sk, j=1:sk
            x=kubPoints[1,i]
            y=kubPoints[2,j]
            f1=(coord[1,1]+(coord[1,2]-coord[1,1])*x+(coord[1,4]-coord[1,1])*y+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*x*y)^2;
            f2=(coord[2,1]+(coord[2,2]-coord[2,1])*x+(coord[2,4]-coord[2,1])*y+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*x*y)^2;
            f3=(coord[3,1]+(coord[3,2]-coord[3,1])*x+(coord[3,4]-coord[3,1])*y+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*x*y)^2;
            f=r*(f1+f2+f3)^(-3/2)
            j11=f*(f2+f3)*((coord[1,2]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*y);
            j21=f*(f1+f3)*((coord[2,2]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*y);
            j31=f*(f1+f2)*((coord[3,2]-coord[3,1])+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*y);
            j12=f*(f2+f3)*((coord[1,4]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*x);
            j22=f*(f1+f3)*((coord[2,4]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*x);
            j32=f*(f1+f2)*((coord[3,4]-coord[3,1])+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*x);
            ddJ[i,j]=sign*1/sqrt((j21*j32-j31*j22)^2+(j31*j12-j11*j32)^2+(j11*j22-j21*j12)^2);
            JtJ11=j11^2+j21^2+j31^2;
            JtJ12=JtJ21=j11*j12+j21*j22+j31*j32;
            JtJ22=j12^2+j22^2+j32^2;
            for k in 1:size(phi,2)
                jphi[1,k][i,j]=(JtJ11*phi[1,k][i,j]+JtJ12*phi[2,k][i,j]);
                jphi[2,k][i,j]=(JtJ21*phi[1,k][i,j]+JtJ22*phi[2,k][i,j]);
            end
        end
    end
    return nothing;
end

function jacobi!(J::Array{Array{Float64,2},2},ddJ::Array{Float64,2},jphi::Array{Array{Float64,2},2},jpsi::Array{Array{Float64,2},2},m::mesh, fid::Int64, kubPoints::Array{Float64,2}, phi::Array{Array{Float64,2},2}, psi::Array{Array{Float64,2},2}, coord::Array{Float64,2})
    key="20";
    mt=m.meshType;
    rstart=m.topology.offset[key][fid];

    z=1;
    for j in rstart:rstart+mt-1
        for i in 1:3
            coord[i,z]=m.geometry.coordinates[i,m.topology.incidence[key][j]];
        end
        z+=1;
    end

    sk=size(kubPoints,2);
    sign=m.orientation[fid];
    if mt==3
        a=coord[:,2]-coord[:,1];
        b=coord[:,3]-coord[:,1];
        for i in 1:sk
            ddJ[i]=1/abs(a[1]*b[2]-b[1]*a[2]);

            J[1,1][i]=a[1];
            J[1,2][i]=b[1];
            J[2,1][i]=a[2];
            J[2,2][i]=b[2];

            for k in 1:size(phi,2)
                jphi[1,k][i]=(J[1,1][i]*phi[1,k][i]+J[1,2][i]*phi[2,k][i]);
                jphi[2,k][i]=(J[2,1][i]*phi[1,k][i]+J[2,2][i]*phi[2,k][i]);
            end
            for k in 1:size(psi,2)
                jpsi[1,k][i]=(J[1,1][i]*psi[1,k][i]+J[1,2][i]*psi[2,k][i]);
                jpsi[2,k][i]=(J[2,1][i]*psi[1,k][i]+J[2,2][i]*psi[2,k][i]);
            end
        end
    elseif mt==4
        r=m.geometry.r[1];
        for i=1:sk, j=1:sk
            x=kubPoints[1,i]
            y=kubPoints[2,j]
            f1=(coord[1,1]+(coord[1,2]-coord[1,1])*x+(coord[1,4]-coord[1,1])*y+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*x*y)^2;
            f2=(coord[2,1]+(coord[2,2]-coord[2,1])*x+(coord[2,4]-coord[2,1])*y+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*x*y)^2;
            f3=(coord[3,1]+(coord[3,2]-coord[3,1])*x+(coord[3,4]-coord[3,1])*y+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*x*y)^2;
            f=r*(f1+f2+f3)^(-3/2)
            J[1,1][i,j]=f*(f2+f3)*((coord[1,2]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*y);
            J[2,1][i,j]=f*(f1+f3)*((coord[2,2]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*y);
            J[3,1][i,j]=f*(f1+f2)*((coord[3,2]-coord[3,1])+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*y);
            J[1,2][i,j]=f*(f2+f3)*((coord[1,4]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*x);
            J[2,2][i,j]=f*(f1+f3)*((coord[2,4]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*x);
            J[3,2][i,j]=f*(f1+f2)*((coord[3,4]-coord[3,1])+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*x);
            ddJ[i,j]=sign*1/sqrt((J[2,1][i,j]*J[3,2][i,j]-J[3,1][i,j]*J[2,2][i,j])^2+(J[3,1][i,j]*J[1,2][i,j]-J[1,1][i,j]*J[3,2][i,j])^2+(J[1,1][i,j]*J[2,2][i,j]-J[2,1][i,j]*J[1,2][i,j])^2);
            for k in 1:size(phi,2)
                jphi[1,k][i,j]=(J[1,1][i,j]*phi[1,k][i,j]+J[1,2][i,j]*phi[2,k][i,j]);
                jphi[2,k][i,j]=(J[2,1][i,j]*phi[1,k][i,j]+J[2,2][i,j]*phi[2,k][i,j]);
                jphi[2,k][i,j]=(J[3,1][i,j]*phi[1,k][i,j]+J[3,2][i,j]*phi[2,k][i,j]);
            end
            for k in 1:size(psi,2)
                jpsi[1,k][i,j]=(J[1,1][i,j]*psi[1,k][i,j]+J[1,2][i,j]*psi[2,k][i,j]);
                jpsi[2,k][i,j]=(J[2,1][i,j]*psi[1,k][i,j]+J[2,2][i,j]*psi[2,k][i,j]);
                jpsi[3,k][i,j]=(J[3,1][i,j]*psi[1,k][i,j]+J[3,2][i,j]*psi[2,k][i,j]);
            end
        end
    end

    return nothing;
end

function jacobi!(J::Array{Array{Float64,1},2},dJ::Array{Float64,1},m::mesh, fid::Int64, kubPoints::Array{Float64,2}, coord::Array{Float64,2})
    key="20";
    mt=m.meshType;
    rstart=m.topology.offset[key][fid];
    z=1;
    for j in rstart:rstart+mt-1
        for i in 1:3
            coord[i,z]=m.geometry.coordinates[i,m.topology.incidence[key][j]];
        end
        z+=1;
    end

    sk=size(kubPoints,2);
    sign=m.orientation[fid];
    if mt==3
        a=coord[:,2]-coord[:,1];
        b=coord[:,3]-coord[:,1];
        for i in 1:sk
            dJ[i]=abs(a[1]*b[2]-b[1]*a[2]);
            J[1,1][i]=a[1];
            J[1,2][i]=b[1];
            J[2,1][i]=a[2];
            J[2,2][i]=b[2];
        end
    elseif mt==4
        r=m.geometry.r[1];
         for i=1:sk
             x=kubPoints[1,i]
             y=kubPoints[2,i]
             f1=(coord[1,1]+(coord[1,2]-coord[1,1])*x+(coord[1,4]-coord[1,1])*y+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*x*y)^2;
             f2=(coord[2,1]+(coord[2,2]-coord[2,1])*x+(coord[2,4]-coord[2,1])*y+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*x*y)^2;
             f3=(coord[3,1]+(coord[3,2]-coord[3,1])*x+(coord[3,4]-coord[3,1])*y+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*x*y)^2;
             f=r*(f1+f2+f3)^(-3/2)
             J[1,1][i]=f*(f2+f3)*((coord[1,2]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*y);
             J[2,1][i]=f*(f1+f3)*((coord[2,2]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*y);
             J[3,1][i]=f*(f1+f2)*((coord[3,2]-coord[3,1])+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*y);
             J[1,2][i]=f*(f2+f3)*((coord[1,4]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*x);
             J[2,2][i]=f*(f1+f3)*((coord[2,4]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*x);
             J[3,2][i]=f*(f1+f2)*((coord[3,4]-coord[3,1])+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*x);
             dJ[i,j]=sign*sqrt((J[2,1][i]*J[3,2][i]-J[3,1][i]*J[2,2][i])^2+(J[3,1][i]*J[1,2][i]-J[1,1][i]*J[3,2][i])^2+(J[1,1][i]*J[2,2][i]-J[2,1][i]*J[1,2][i])^2);
         end
    end
    return nothing;
end

function jacobi!(J::Array{Array{Float64,1},2},ddJ::Array{Float64,1},jphi::Array{Array{Float64,1},2},jpsi::Array{Array{Float64,1},2},m::mesh, fid::Int64, kubPoints::Array{Float64,2}, phi::Array{Array{Float64,1},2}, psi::Array{Array{Float64,1},2}, coord::Array{Float64,2})
    key="20";
    mt=m.meshType;
    rstart=m.topology.offset[key][fid];
    z=1;
    for j in rstart:rstart+mt-1
        for i in 1:3
            coord[i,z]=m.geometry.coordinates[i,m.topology.incidence[key][j]];
        end
        z+=1;
    end

    sk=size(kubPoints,2);
    sign=m.orientation[fid];
    if mt==3
        a=coord[:,2]-coord[:,1];
        b=coord[:,3]-coord[:,1];
        for i in 1:sk
            ddJ[i]=1/abs(a[1]*b[2]-b[1]*a[2]);
            J[1,1][i]=a[1];
            J[1,2][i]=b[1];
            J[2,1][i]=a[2];
            J[2,2][i]=b[2];

            for k in 1:size(phi,2)
                jphi[1,k][i]=(J[1,1][i]*phi[1,k][i]+J[1,2][i]*phi[2,k][i]);
                jphi[2,k][i]=(J[2,1][i]*phi[1,k][i]+J[2,2][i]*phi[2,k][i]);
            end
            for k in 1:size(psi,2)
                jpsi[1,k][i]=(J[1,1][i]*psi[1,k][i]+J[1,2][i]*psi[2,k][i]);
                jpsi[2,k][i]=(J[2,1][i]*psi[1,k][i]+J[2,2][i]*psi[2,k][i]);
            end
        end
    elseif mt==4
        r=m.geometry.r[1];
        for i=1:sk
            x=kubPoints[1,i]
            y=kubPoints[2,i]
            f1=(coord[1,1]+(coord[1,2]-coord[1,1])*x+(coord[1,4]-coord[1,1])*y+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*x*y)^2;
            f2=(coord[2,1]+(coord[2,2]-coord[2,1])*x+(coord[2,4]-coord[2,1])*y+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*x*y)^2;
            f3=(coord[3,1]+(coord[3,2]-coord[3,1])*x+(coord[3,4]-coord[3,1])*y+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*x*y)^2;
            f=r*(f1+f2+f3)^(-3/2)
            J[1,1][i]=f*(f2+f3)*((coord[1,2]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*y);
            J[2,1][i]=f*(f1+f3)*((coord[2,2]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*y);
            J[3,1][i]=f*(f1+f2)*((coord[3,2]-coord[3,1])+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*y);
            J[1,2][i]=f*(f2+f3)*((coord[1,4]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*x);
            J[2,2][i]=f*(f1+f3)*((coord[2,4]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*x);
            J[3,2][i]=f*(f1+f2)*((coord[3,4]-coord[3,1])+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*x);
            ddJ[i]=sign*1/sqrt((J[2,1][i]*J[3,2][i]-J[3,1][i]*J[2,2][i])^2+(J[3,1][i]*J[1,2][i]-J[1,1][i]*J[3,2][i])^2+(J[1,1][i]*J[2,2][i]-J[2,1][i]*J[1,2][i])^2);
            for k in 1:size(phi,2)
                jphi[1,k][i]=(J[1,1][i]*phi[1,k][i]+J[1,2][i]*phi[2,k][i]);
                jphi[2,k][i]=(J[2,1][i]*phi[1,k][i]+J[2,2][i]*phi[2,k][i]);
                jphi[3,k][i]=(J[3,1][i]*phi[1,k][i]+J[3,2][i]*phi[2,k][i]);
            end
            for k in 1:size(psi,2)
                jpsi[1,k][i]=(J[1,1][i]*psi[1,k][i]+J[1,2][i]*psi[2,k][i]);
                jpsi[2,k][i]=(J[2,1][i]*psi[1,k][i]+J[2,2][i]*psi[2,k][i]);
                jpsi[2,k][i]=(J[3,1][i]*psi[1,k][i]+J[3,2][i]*psi[2,k][i]);
            end
        end
    end

    return nothing;
end
