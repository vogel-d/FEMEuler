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
        for i in 1:2
            coord[i,z]=m.geometry.coordinates[i,m.topology.incidence[key][j]];
        end
        z+=1;
    end
    sk=size(kubPoints,2);

    if mt==3
        a=coord[:,2]-coord[:,1];
        b=coord[:,3]-coord[:,1];
        ones=ones(sk,sk);
        dJ=abs(a[1]*b[2]-b[1]*a[2])*ones;
        J[1,1]=a[1]*ones;
        J[1,2]=b[1]*ones;
        J[2,1]=a[2]*ones;
        J[2,2]=b[2]*ones;
    elseif mt==4
         #J[1,1]=Array{Float64,2}(undef,sk,sk);
         #J[1,2]=Array{Float64,2}(undef,sk,sk);
         #J[2,1]=Array{Float64,2}(undef,sk,sk);
         #J[2,2]=Array{Float64,2}(undef,sk,sk);
         for i=1:sk, j=1:sk
             J[1,1][i,j]=(coord[1,2]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*kubPoints[2,j];
             J[2,1][i,j]=(coord[2,2]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*kubPoints[2,j];
             J[1,2][i,j]=(coord[1,4]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*kubPoints[1,i];
             J[2,2][i,j]=(coord[2,4]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*kubPoints[1,i];
             dJ[i,j]=abs(J[1,1][i,j]*J[2,2][i,j]-J[2,1][i,j]*J[1,2][i,j]);
         end
    end
    return nothing;
end

function jacobi!(J::Array{Array{Float64,2},2},ddJ::Array{Float64,2},jphi::Array{Float64,4},m::mesh, fid::Int64, kubPoints::Array{Float64,2}, phi::Array{Float64,4}, coord::Array{Float64,2})
    key="20";
    mt=m.meshType;
    rstart=m.topology.offset[key][fid];
    z=1;
    for j in rstart:rstart+mt-1
        for i in 1:2
            coord[i,z]=m.geometry.coordinates[i,m.topology.incidence[key][j]];
        end
        z+=1;
    end

    sk=size(kubPoints,2);

    if mt==3
        a=coord[:,2]-coord[:,1];
        b=coord[:,3]-coord[:,1];
        ones=ones(sk,sk);
        ddJ=1/abs(a[1]*b[2]-b[1]*a[2])*ones;
        J[1,1]=a[1]*ones;
        J[1,2]=b[1]*ones;
        J[2,1]=a[2]*ones;
        J[2,2]=b[2]*ones;
    elseif mt==4
        #J[1,1]=Array{Float64,2}(undef,sk,sk);
        #J[1,2]=Array{Float64,2}(undef,sk,sk);
        #J[2,1]=Array{Float64,2}(undef,sk,sk);
        #J[2,2]=Array{Float64,2}(undef,sk,sk);
        for i=1:sk, j=1:sk
            J[1,1][i,j]=(coord[1,2]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*kubPoints[2,j];
            J[2,1][i,j]=(coord[2,2]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*kubPoints[2,j];
            J[1,2][i,j]=(coord[1,4]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*kubPoints[1,i];
            J[2,2][i,j]=(coord[2,4]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*kubPoints[1,i];
            ddJ[i,j]=1/abs(J[1,1][i,j]*J[2,2][i,j]-J[2,1][i,j]*J[1,2][i,j]);
            for k in 1:size(phi,4)
                jphi[1,i,j,k]=(J[1,1][i,j]*phi[1,i,j,k]+J[1,2][i,j]*phi[2,i,j,k]);
                jphi[2,i,j,k]=(J[2,1][i,j]*phi[1,i,j,k]+J[2,2][i,j]*phi[2,i,j,k]);
            end
        end
    end

    return nothing;
end

function jacobi!(ddJ::Array{Float64,2}, jphi::Array{Float64,4}, m::mesh, fid::Int64, kubPoints::Array{Float64,2}, phi::Array{Float64,4}, coord::Array{Float64,2})
    key="20";
    mt=m.meshType;
    rstart=m.topology.offset[key][fid];

    z=1;
    for j in rstart:rstart+mt-1
        for i in 1:2
            coord[i,z]=m.geometry.coordinates[i,m.topology.incidence[key][j]];
        end
        z+=1;
    end
    sk=size(kubPoints,2);

    if mt==3
        a=coord[:,2]-coord[:,1];
        b=coord[:,3]-coord[:,1];
        ones=ones(sk,sk)
        ddJ=1/abs(a[1]*b[2]-b[1]*a[2])*ones;
        JtJ11=(a[1]^2+a[2]^2)*ones;
        JtJ12=(a[1]*b[1]+a[2]*b[2])*ones;
        JtJ2=(a[1]*b[1]+a[2]*b[2])*ones;
        JtJ22=(b[1]^2+b[2]^2)*ones;
        for k in 1:size(phi,2)
            jphi[1,k][i,j]=(JtJ11*phi[1,k][i,j]+JtJ12*phi[2,k][i,j]);
            jphi[2,k][i,j]=(JtJ21*phi[1,k][i,j]+JtJ22*phi[2,k][i,j]);
        end
    elseif mt==4
        for i=1:sk, j=1:sk
            a1=coord[1,2]-coord[1,1]+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*kubPoints[2,j];
            b1=coord[1,4]-coord[1,1]+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*kubPoints[1,i];
            a2=coord[2,2]-coord[2,1]+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*kubPoints[2,j];
            b2=coord[2,4]-coord[2,1]+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*kubPoints[1,i];
            ddJ[i,j]=1/abs(a1*b2-b1*a2);
            JtJ11=a1^2+a2^2;
            JtJ12=a1*b1+a2*b2;
            JtJ21=a1*b1+a2*b2;
            JtJ22=b1^2+b2^2;
            for k in 1:size(phi,4)
                jphi[1,i,j,k]=(JtJ11*phi[1,i,j,k]+JtJ12*phi[2,i,j,k]);
                jphi[2,i,j,k]=(JtJ21*phi[1,i,j,k]+JtJ22*phi[2,i,j,k]);
            end
        end
    end
    return nothing;
end

function jacobi!(J::Array{Array{Float64,2},2},ddJ::Array{Float64,2},jphi::Array{Float64,4},jpsi::Array{Float64,4},m::mesh, fid::Int64, kubPoints::Array{Float64,2}, phi::Array{Float64,4}, psi::Array{Float64,4}, coord::Array{Float64,2})
    key="20";
    mt=m.meshType;
    rstart=m.topology.offset[key][fid];

    z=1;
    for j in rstart:rstart+mt-1
        for i in 1:2
            coord[i,z]=m.geometry.coordinates[i,m.topology.incidence[key][j]];
        end
        z+=1;
    end
    sk=size(kubPoints,2);

    sk=size(kubPoints,2);

    if mt==3
        a=coord[:,2]-coord[:,1];
        b=coord[:,3]-coord[:,1];
        ones=ones(sk,sk);
        ddJ=1/abs(a[1]*b[2]-b[1]*a[2])*ones;
        J[1,1]=a[1]*ones;
        J[1,2]=b[1]*ones;
        J[2,1]=a[2]*ones;
        J[2,2]=b[2]*ones;
    elseif mt==4
        #J[1,1]=Array{Float64,2}(undef,sk,sk);
        #J[1,2]=Array{Float64,2}(undef,sk,sk);
        #J[2,1]=Array{Float64,2}(undef,sk,sk);
        #J[2,2]=Array{Float64,2}(undef,sk,sk);
        for i=1:sk, j=1:sk
            J[1,1][i,j]=(coord[1,2]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*kubPoints[2,j];
            J[2,1][i,j]=(coord[2,2]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*kubPoints[2,j];
            J[1,2][i,j]=(coord[1,4]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*kubPoints[1,i];
            J[2,2][i,j]=(coord[2,4]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*kubPoints[1,i];
            ddJ[i,j]=1/abs(J[1,1][i,j]*J[2,2][i,j]-J[2,1][i,j]*J[1,2][i,j]);
            for k in 1:size(phi,4)
                jphi[1,i,j,k]=(J[1,1][i,j]*phi[1,i,j,k]+J[1,2][i,j]*phi[2,i,j,k]);
                jphi[2,i,j,k]=(J[2,1][i,j]*phi[1,i,j,k]+J[2,2][i,j]*phi[2,i,j,k]);
            end
            for k in 1:size(psi,4)
                jpsi[1,i,j,k]=(J[1,1][i,j]*psi[1,i,j,k]+J[1,2][i,j]*psi[2,i,j,k]);
                jpsi[2,i,j,k]=(J[2,1][i,j]*psi[1,i,j,k]+J[2,2][i,j]*psi[2,i,j,k]);
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
        for i in 1:2
            coord[i,z]=m.geometry.coordinates[i,m.topology.incidence[key][j]];
        end
        z+=1;
    end

    sk=size(kubPoints,2);

    if mt==3
        a=coord[:,2]-coord[:,1];
        b=coord[:,3]-coord[:,1];
        ones=ones(sk,sk);
        dJ=abs(a[1]*b[2]-b[1]*a[2])*ones;
        J[1,1]=a[1]*ones;
        J[1,2]=b[1]*ones;
        J[2,1]=a[2]*ones;
        J[2,2]=b[2]*ones;
    elseif mt==4
         #J[1,1]=Array{Float64,1}(undef,sk);
         #J[1,2]=Array{Float64,1}(undef,sk);
         #J[2,1]=Array{Float64,1}(undef,sk);
         #J[2,2]=Array{Float64,1}(undef,sk);
         for i=1:sk
             J[1,1][i]=(coord[1,2]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*kubPoints[2,i];
             J[2,1][i]=(coord[2,2]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*kubPoints[2,i];
             J[1,2][i]=(coord[1,4]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*kubPoints[1,i];
             J[2,2][i]=(coord[2,4]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*kubPoints[1,i];
             dJ[i]=abs(J[1,1][i]*J[2,2][i]-J[2,1][i]*J[1,2][i]);
         end
    end
    return nothing;
end

function jacobi!(J::Array{Array{Float64,1},2},ddJ::Array{Float64,1},jphi::Array{Float64,3},jpsi::Array{Float64,3},m::mesh, fid::Int64, kubPoints::Array{Float64,2}, phi, psi, coord::Array{Float64,2})
    key="20";
    mt=m.meshType;
    rstart=m.topology.offset[key][fid];
    z=1;
    for j in rstart:rstart+mt-1
        for i in 1:2
            coord[i,z]=m.geometry.coordinates[i,m.topology.incidence[key][j]];
        end
        z+=1;
    end

    sk=size(kubPoints,2);

    if mt==3
        a=coord[:,2]-coord[:,1];
        b=coord[:,3]-coord[:,1];
        ones=ones(sk);
        ddJ=1/abs(a[1]*b[2]-b[1]*a[2])*ones;
        J[1,1]=a[1]*ones;
        J[1,2]=b[1]*ones;
        J[2,1]=a[2]*ones;
        J[2,2]=b[2]*ones;
    elseif mt==4
        #J[1,1]=Array{Float64,1}(undef,sk);
        #J[1,2]=Array{Float64,1}(undef,sk);
        #J[2,1]=Array{Float64,1}(undef,sk);
        #J[2,2]=Array{Float64,1}(undef,sk);
        for i=1:sk
            J[1,1][i]=(coord[1,2]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*kubPoints[2,i];
            J[2,1][i]=(coord[2,2]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*kubPoints[2,i];
            J[1,2][i]=(coord[1,4]-coord[1,1])+(coord[1,3]-coord[1,4]-coord[1,2]+coord[1,1])*kubPoints[1,i];
            J[2,2][i]=(coord[2,4]-coord[2,1])+(coord[2,3]-coord[2,4]-coord[2,2]+coord[2,1])*kubPoints[1,i];
            ddJ[i]=1/abs(J[1,1][i]*J[2,2][i]-J[2,1][i]*J[1,2][i]);
            for k in 1:size(phi,3)
                jphi[1,i,k]=(J[1,1][i]*phi[1,i,k]+J[1,2][i]*phi[2,i,k]);
                jphi[2,i,k]=(J[2,1][i]*phi[1,i,k]+J[2,2][i]*phi[2,i,k]);
            end
            for k in 1:size(psi,3)
                jpsi[1,i,k]=(J[1,1][i]*psi[1,i,k]+J[1,2][i]*psi[2,i,k]);
                jpsi[2,i,k]=(J[2,1][i]*psi[1,i,k]+J[2,2][i]*psi[2,i,k]);
            end
        end
    end

    return nothing;
end
