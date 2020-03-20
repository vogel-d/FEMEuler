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
             X1=kubPoints[1,i]
             X2=kubPoints[2,j]
             ξ1=coord[1,1]+(coord[1,2]-coord[1,1])*X1+(coord[1,4]-coord[1,1])*X2+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*X1*X2;
             ξ2=coord[2,1]+(coord[2,2]-coord[2,1])*X1+(coord[2,4]-coord[2,1])*X2+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*X1*X2;
             ξ3=coord[3,1]+(coord[3,2]-coord[3,1])*X1+(coord[3,4]-coord[3,1])*X2+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*X1*X2;
             f=r*(ξ1^2+ξ2^2+ξ3^2)^(-3/2)

             dx1dξ1=f*(ξ2^2+ξ3^2); dx1dξ2=-f*ξ1*ξ2; dx1dξ3=-f*ξ1*ξ3;
             dx2dξ1=dx1dξ2; dx2dξ2=f*(ξ1^2+ξ3^2); dx2dξ3=-f*ξ2*ξ3;
             dx3dξ1=dx1dξ3; dx3dξ2=dx2dξ3; dx3dξ3=f*(ξ1^2+ξ2^2)

             dξ1dX1=(coord[1,2]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*X2;
             dξ1dX2=(coord[1,4]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*X1;
             dξ2dX1=(coord[2,2]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*X2;
             dξ2dX2=(coord[2,4]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*X1;
             dξ3dX1=(coord[3,2]-coord[3,1])+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*X2;
             dξ3dX2=(coord[3,4]-coord[3,1])+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*X1;

             J[1,1][i,j]=dx1dξ1*dξ1dX1+dx1dξ2*dξ2dX1+dx1dξ3*dξ3dX1;
             J[2,1][i,j]=dx2dξ1*dξ1dX1+dx2dξ2*dξ2dX1+dx2dξ3*dξ3dX1;
             J[3,1][i,j]=dx3dξ1*dξ1dX1+dx3dξ2*dξ2dX1+dx3dξ3*dξ3dX1;
             J[1,2][i,j]=dx1dξ1*dξ1dX2+dx1dξ2*dξ2dX2+dx1dξ3*dξ3dX2;
             J[2,2][i,j]=dx2dξ1*dξ1dX2+dx2dξ2*dξ2dX2+dx2dξ3*dξ3dX2;
             J[3,2][i,j]=dx3dξ1*dξ1dX2+dx3dξ2*dξ2dX2+dx3dξ3*dξ3dX2;

             dJ[i,j]=sign*sqrt((J[2,1][i,j]*J[3,2][i,j]-J[3,1][i,j]*J[2,2][i,j])^2+(J[3,1][i,j]*J[1,2][i,j]-J[1,1][i,j]*J[3,2][i,j])^2+(J[1,1][i,j]*J[2,2][i,j]-J[2,1][i,j]*J[1,2][i,j])^2);
         end
    end
    return nothing;
end

function jacobi!(J::Array{Float64,2},m::mesh, fid::Int64, X1::Float64, X2::Float64, coord::Array{Float64,2})
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
        ξ1=coord[1,1]+(coord[1,2]-coord[1,1])*X1+(coord[1,4]-coord[1,1])*X2+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*X1*X2;
        ξ2=coord[2,1]+(coord[2,2]-coord[2,1])*X1+(coord[2,4]-coord[2,1])*X2+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*X1*X2;
        ξ3=coord[3,1]+(coord[3,2]-coord[3,1])*X1+(coord[3,4]-coord[3,1])*X2+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*X1*X2;
        f=m.geometry.r[1]*(ξ1^2+ξ2^2+ξ3^2)^(-3/2)

        dx1dξ1=f*(ξ2^2+ξ3^2); dx1dξ2=-f*ξ1*ξ2; dx1dξ3=-f*ξ1*ξ3;
        dx2dξ1=dx1dξ2; dx2dξ2=f*(ξ1^2+ξ3^2); dx2dξ3=-f*ξ2*ξ3;
        dx3dξ1=dx1dξ3; dx3dξ2=dx2dξ3; dx3dξ3=f*(ξ1^2+ξ2^2)

        dξ1dX1=(coord[1,2]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*X2;
        dξ1dX2=(coord[1,4]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*X1;
        dξ2dX1=(coord[2,2]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*X2;
        dξ2dX2=(coord[2,4]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*X1;
        dξ3dX1=(coord[3,2]-coord[3,1])+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*X2;
        dξ3dX2=(coord[3,4]-coord[3,1])+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*X1;

        J[1,1]=dx1dξ1*dξ1dX1+dx1dξ2*dξ2dX1+dx1dξ3*dξ3dX1;
        J[2,1]=dx2dξ1*dξ1dX1+dx2dξ2*dξ2dX1+dx2dξ3*dξ3dX1;
        J[3,1]=dx3dξ1*dξ1dX1+dx3dξ2*dξ2dX1+dx3dξ3*dξ3dX1;
        J[1,2]=dx1dξ1*dξ1dX2+dx1dξ2*dξ2dX2+dx1dξ3*dξ3dX2;
        J[2,2]=dx2dξ1*dξ1dX2+dx2dξ2*dξ2dX2+dx2dξ3*dξ3dX2;
        J[3,2]=dx3dξ1*dξ1dX2+dx3dξ2*dξ2dX2+dx3dξ3*dξ3dX2;

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
            X1=kubPoints[1,i]
            X2=kubPoints[2,j]
            ξ1=coord[1,1]+(coord[1,2]-coord[1,1])*X1+(coord[1,4]-coord[1,1])*X2+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*X1*X2;
            ξ2=coord[2,1]+(coord[2,2]-coord[2,1])*X1+(coord[2,4]-coord[2,1])*X2+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*X1*X2;
            ξ3=coord[3,1]+(coord[3,2]-coord[3,1])*X1+(coord[3,4]-coord[3,1])*X2+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*X1*X2;
            f=r*(ξ1^2+ξ2^2+ξ3^2)^(-3/2)

            dx1dξ1=f*(ξ2^2+ξ3^2); dx1dξ2=-f*ξ1*ξ2; dx1dξ3=-f*ξ1*ξ3;
            dx2dξ1=dx1dξ2; dx2dξ2=f*(ξ1^2+ξ3^2); dx2dξ3=-f*ξ2*ξ3;
            dx3dξ1=dx1dξ3; dx3dξ2=dx2dξ3; dx3dξ3=f*(ξ1^2+ξ2^2)

            dξ1dX1=(coord[1,2]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*X2;
            dξ1dX2=(coord[1,4]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*X1;
            dξ2dX1=(coord[2,2]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*X2;
            dξ2dX2=(coord[2,4]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*X1;
            dξ3dX1=(coord[3,2]-coord[3,1])+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*X2;
            dξ3dX2=(coord[3,4]-coord[3,1])+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*X1;

            J[1,1][i,j]=dx1dξ1*dξ1dX1+dx1dξ2*dξ2dX1+dx1dξ3*dξ3dX1;
            J[2,1][i,j]=dx2dξ1*dξ1dX1+dx2dξ2*dξ2dX1+dx2dξ3*dξ3dX1;
            J[3,1][i,j]=dx3dξ1*dξ1dX1+dx3dξ2*dξ2dX1+dx3dξ3*dξ3dX1;
            J[1,2][i,j]=dx1dξ1*dξ1dX2+dx1dξ2*dξ2dX2+dx1dξ3*dξ3dX2;
            J[2,2][i,j]=dx2dξ1*dξ1dX2+dx2dξ2*dξ2dX2+dx2dξ3*dξ3dX2;
            J[3,2][i,j]=dx3dξ1*dξ1dX2+dx3dξ2*dξ2dX2+dx3dξ3*dξ3dX2;

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
            X1=kubPoints[1,i]
            X2=kubPoints[2,j]
            ξ1=coord[1,1]+(coord[1,2]-coord[1,1])*X1+(coord[1,4]-coord[1,1])*X2+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*X1*X2;
            ξ2=coord[2,1]+(coord[2,2]-coord[2,1])*X1+(coord[2,4]-coord[2,1])*X2+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*X1*X2;
            ξ3=coord[3,1]+(coord[3,2]-coord[3,1])*X1+(coord[3,4]-coord[3,1])*X2+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*X1*X2;
            f=r*(ξ1^2+ξ2^2+ξ3^2)^(-3/2)

            dx1dξ1=f*(ξ2^2+ξ3^2); dx1dξ2=-f*ξ1*ξ2; dx1dξ3=-f*ξ1*ξ3;
            dx2dξ1=dx1dξ2; dx2dξ2=f*(ξ1^2+ξ3^2); dx2dξ3=-f*ξ2*ξ3;
            dx3dξ1=dx1dξ3; dx3dξ2=dx2dξ3; dx3dξ3=f*(ξ1^2+ξ2^2)

            dξ1dX1=(coord[1,2]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*X2;
            dξ1dX2=(coord[1,4]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*X1;
            dξ2dX1=(coord[2,2]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*X2;
            dξ2dX2=(coord[2,4]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*X1;
            dξ3dX1=(coord[3,2]-coord[3,1])+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*X2;
            dξ3dX2=(coord[3,4]-coord[3,1])+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*X1;

            j11=dx1dξ1*dξ1dX1+dx1dξ2*dξ2dX1+dx1dξ3*dξ3dX1;
            j21=dx2dξ1*dξ1dX1+dx2dξ2*dξ2dX1+dx2dξ3*dξ3dX1;
            j31=dx3dξ1*dξ1dX1+dx3dξ2*dξ2dX1+dx3dξ3*dξ3dX1;
            j12=dx1dξ1*dξ1dX2+dx1dξ2*dξ2dX2+dx1dξ3*dξ3dX2;
            j22=dx2dξ1*dξ1dX2+dx2dξ2*dξ2dX2+dx2dξ3*dξ3dX2;
            j32=dx3dξ1*dξ1dX2+dx3dξ2*dξ2dX2+dx3dξ3*dξ3dX2;

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
            X1=kubPoints[1,i]
            X2=kubPoints[2,j]
            ξ1=coord[1,1]+(coord[1,2]-coord[1,1])*X1+(coord[1,4]-coord[1,1])*X2+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*X1*X2;
            ξ2=coord[2,1]+(coord[2,2]-coord[2,1])*X1+(coord[2,4]-coord[2,1])*X2+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*X1*X2;
            ξ3=coord[3,1]+(coord[3,2]-coord[3,1])*X1+(coord[3,4]-coord[3,1])*X2+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*X1*X2;
            f=r*(ξ1^2+ξ2^2+ξ3^2)^(-3/2)

            dx1dξ1=f*(ξ2^2+ξ3^2); dx1dξ2=-f*ξ1*ξ2; dx1dξ3=-f*ξ1*ξ3;
            dx2dξ1=dx1dξ2; dx2dξ2=f*(ξ1^2+ξ3^2); dx2dξ3=-f*ξ2*ξ3;
            dx3dξ1=dx1dξ3; dx3dξ2=dx2dξ3; dx3dξ3=f*(ξ1^2+ξ2^2)

            dξ1dX1=(coord[1,2]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*X2;
            dξ1dX2=(coord[1,4]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*X1;
            dξ2dX1=(coord[2,2]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*X2;
            dξ2dX2=(coord[2,4]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*X1;
            dξ3dX1=(coord[3,2]-coord[3,1])+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*X2;
            dξ3dX2=(coord[3,4]-coord[3,1])+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*X1;

            J[1,1][i,j]=dx1dξ1*dξ1dX1+dx1dξ2*dξ2dX1+dx1dξ3*dξ3dX1;
            J[2,1][i,j]=dx2dξ1*dξ1dX1+dx2dξ2*dξ2dX1+dx2dξ3*dξ3dX1;
            J[3,1][i,j]=dx3dξ1*dξ1dX1+dx3dξ2*dξ2dX1+dx3dξ3*dξ3dX1;
            J[1,2][i,j]=dx1dξ1*dξ1dX2+dx1dξ2*dξ2dX2+dx1dξ3*dξ3dX2;
            J[2,2][i,j]=dx2dξ1*dξ1dX2+dx2dξ2*dξ2dX2+dx2dξ3*dξ3dX2;
            J[3,2][i,j]=dx3dξ1*dξ1dX2+dx3dξ2*dξ2dX2+dx3dξ3*dξ3dX2;

            ddJ[i,j]=sign*1/sqrt((J[2,1][i,j]*J[3,2][i,j]-J[3,1][i,j]*J[2,2][i,j])^2+(J[3,1][i,j]*J[1,2][i,j]-J[1,1][i,j]*J[3,2][i,j])^2+(J[1,1][i,j]*J[2,2][i,j]-J[2,1][i,j]*J[1,2][i,j])^2);
            for k in 1:size(phi,2)
                jphi[1,k][i,j]=(J[1,1][i,j]*phi[1,k][i,j]+J[1,2][i,j]*phi[2,k][i,j]);
                jphi[2,k][i,j]=(J[2,1][i,j]*phi[1,k][i,j]+J[2,2][i,j]*phi[2,k][i,j]);
                jphi[3,k][i,j]=(J[3,1][i,j]*phi[1,k][i,j]+J[3,2][i,j]*phi[2,k][i,j]);
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

function jacobi!(ddJe::Array{Float64,1},m::mesh, fid::Int64,n,kubPoints::Array{Float64,2}, coord::Array{Float64,2})
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
             X1=kubPoints[1,i]
             X2=kubPoints[2,1]
             ξ1=coord[1,1]+(coord[1,2]-coord[1,1])*X1+(coord[1,4]-coord[1,1])*X2+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*X1*X2;
             ξ2=coord[2,1]+(coord[2,2]-coord[2,1])*X1+(coord[2,4]-coord[2,1])*X2+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*X1*X2;
             ξ3=coord[3,1]+(coord[3,2]-coord[3,1])*X1+(coord[3,4]-coord[3,1])*X2+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*X1*X2;
             f=r*(ξ1^2+ξ2^2+ξ3^2)^(-3/2)

             dx1dξ1=f*(ξ2^2+ξ3^2); dx1dξ2=-f*ξ1*ξ2; dx1dξ3=-f*ξ1*ξ3;
             dx2dξ1=dx1dξ2; dx2dξ2=f*(ξ1^2+ξ3^2); dx2dξ3=-f*ξ2*ξ3;
             dx3dξ1=dx1dξ3; dx3dξ2=dx2dξ3; dx3dξ3=f*(ξ1^2+ξ2^2)

             dξ1dX1=(coord[1,2]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*X2;
             dξ1dX2=(coord[1,4]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*X1;
             dξ2dX1=(coord[2,2]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*X2;
             dξ2dX2=(coord[2,4]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*X1;
             dξ3dX1=(coord[3,2]-coord[3,1])+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*X2;
             dξ3dX2=(coord[3,4]-coord[3,1])+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*X1;

             j11=dx1dξ1*dξ1dX1+dx1dξ2*dξ2dX1+dx1dξ3*dξ3dX1;
             j21=dx2dξ1*dξ1dX1+dx2dξ2*dξ2dX1+dx2dξ3*dξ3dX1;
             j31=dx3dξ1*dξ1dX1+dx3dξ2*dξ2dX1+dx3dξ3*dξ3dX1;
             j12=dx1dξ1*dξ1dX2+dx1dξ2*dξ2dX2+dx1dξ3*dξ3dX2;
             j22=dx2dξ1*dξ1dX2+dx2dξ2*dξ2dX2+dx2dξ3*dξ3dX2;
             j32=dx3dξ1*dξ1dX2+dx3dξ2*dξ2dX2+dx3dξ3*dξ3dX2;

             JtJ11=j11^2+j21^2+j31^2;
             JtJ12=JtJ21=j11*j12+j21*j22+j31*j32;
             JtJ22=j12^2+j22^2+j32^2;

             fj=1/(JtJ11*JtJ22-JtJ12^2)
             jInvT11=fj*(JtJ22*j11-JtJ12*j12)
             jInvT12=fj*(JtJ11*j12-JtJ12*j11)
             jInvT21=fj*(JtJ22*j21-JtJ12*j22)
             jInvT22=fj*(JtJ11*j22-JtJ12*j21)
             jInvT31=fj*(JtJ22*j31-JtJ12*j32)
             jInvT32=fj*(JtJ11*j32-JtJ12*j31)

             dJ=sqrt((j21*j32-j31*j22)^2+(j31*j12-j11*j32)^2+(j11*j22-j21*j12)^2);
             ddJe[i]=1/(dJ*sqrt((jInvT11*n[1]+jInvT12*n[2])^2+(jInvT21*n[1]+jInvT22*n[2])^2+(jInvT31*n[1]+jInvT32*n[2])^2))
         end
    end
    return nothing;
end

function jacobi!(ddJ::Array{Float64,1},ddJe::Array{Float64,1},jphi::Array{Array{Float64,1},2},jpsi::Array{Array{Float64,1},2},m::mesh, fid::Int64, n, kubPoints::Array{Float64,2}, phi::Array{Array{Float64,1},2}, psi::Array{Array{Float64,1},2}, coord::Array{Float64,2})
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
            X1=kubPoints[1,i]
            X2=kubPoints[2,i]
            ξ1=coord[1,1]+(coord[1,2]-coord[1,1])*X1+(coord[1,4]-coord[1,1])*X2+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*X1*X2;
            ξ2=coord[2,1]+(coord[2,2]-coord[2,1])*X1+(coord[2,4]-coord[2,1])*X2+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*X1*X2;
            ξ3=coord[3,1]+(coord[3,2]-coord[3,1])*X1+(coord[3,4]-coord[3,1])*X2+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*X1*X2;
            f=r*(ξ1^2+ξ2^2+ξ3^2)^(-3/2)

            dx1dξ1=f*(ξ2^2+ξ3^2); dx1dξ2=-f*ξ1*ξ2; dx1dξ3=-f*ξ1*ξ3;
            dx2dξ1=dx1dξ2; dx2dξ2=f*(ξ1^2+ξ3^2); dx2dξ3=-f*ξ2*ξ3;
            dx3dξ1=dx1dξ3; dx3dξ2=dx2dξ3; dx3dξ3=f*(ξ1^2+ξ2^2)

            dξ1dX1=(coord[1,2]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*X2;
            dξ1dX2=(coord[1,4]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*X1;
            dξ2dX1=(coord[2,2]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*X2;
            dξ2dX2=(coord[2,4]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*X1;
            dξ3dX1=(coord[3,2]-coord[3,1])+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*X2;
            dξ3dX2=(coord[3,4]-coord[3,1])+(coord[3,3]-coord[3,4]-coord[3,2]+coord[3,1])*X1;

            j11=dx1dξ1*dξ1dX1+dx1dξ2*dξ2dX1+dx1dξ3*dξ3dX1;
            j21=dx2dξ1*dξ1dX1+dx2dξ2*dξ2dX1+dx2dξ3*dξ3dX1;
            j31=dx3dξ1*dξ1dX1+dx3dξ2*dξ2dX1+dx3dξ3*dξ3dX1;
            j12=dx1dξ1*dξ1dX2+dx1dξ2*dξ2dX2+dx1dξ3*dξ3dX2;
            j22=dx2dξ1*dξ1dX2+dx2dξ2*dξ2dX2+dx2dξ3*dξ3dX2;
            j32=dx3dξ1*dξ1dX2+dx3dξ2*dξ2dX2+dx3dξ3*dξ3dX2;

            JtJ11=j11^2+j21^2+j31^2;
            JtJ12=JtJ21=j11*j12+j21*j22+j31*j32;
            JtJ22=j12^2+j22^2+j32^2;

            fj=1/(JtJ11*JtJ22-JtJ12^2)
            jInvT11=fj*(JtJ22*j11-JtJ12*j12)
            jInvT12=fj*(JtJ11*j12-JtJ12*j11)
            jInvT21=fj*(JtJ22*j21-JtJ12*j22)
            jInvT22=fj*(JtJ11*j22-JtJ12*j21)
            jInvT31=fj*(JtJ22*j31-JtJ12*j32)
            jInvT32=fj*(JtJ11*j32-JtJ12*j31)

            dJ=sqrt((j21*j32-j31*j22)^2+(j31*j12-j11*j32)^2+(j11*j22-j21*j12)^2);
            ddJe[i]=1/(dJ*sqrt((jInvT11*n[1]+jInvT12*n[2])^2+(jInvT21*n[1]+jInvT22*n[2])^2+(jInvT31*n[1]+jInvT32*n[2])^2))
            ddJ[i]=sign*1/dJ;
            for k in 1:size(phi,2)
                jphi[1,k][i]=(j11*phi[1,k][i]+j12*phi[2,k][i]);
                jphi[2,k][i]=(j21*phi[1,k][i]+j22*phi[2,k][i]);
                jphi[3,k][i]=(j31*phi[1,k][i]+j32*phi[2,k][i]);
            end
            for k in 1:size(psi,2)
                jpsi[1,k][i]=(j11*psi[1,k][i]+j12*psi[2,k][i]);
                jpsi[2,k][i]=(j21*psi[1,k][i]+j22*psi[2,k][i]);
                jpsi[3,k][i]=(j31*psi[1,k][i]+j32*psi[2,k][i]);
            end
        end
    end

    return nothing;
end
