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

function jacobi!(J::Array{Array{Float64,1},2},ddJ::Array{Float64,1},jphi::Array{Float64,3},jpsi::Array{Float64,3},m::mesh, fid::Int64, kubPoints::Array{Float64,2}, phi::Array{Float64,3}, psi::Array{Float64,3}, coord::Array{Float64,2})
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


#=
function initPhi(size::Tuple{Int64,Int64}, sk::Tuple{Int64,Int64})
    phi=Array{Float64,4}(undef,sk[1],sk[2],size[1],size[2]);
    return phi;
end

function initPhi(size::Tuple{Int64,Int64}, sk::Int64)
    phi=Array{Float64,3}(undef,sk,size[1],size[2]);
    return phi;
end
=#

function coordTrans(mt::Int64, normals::Array{Float64,2}, type::Array{Symbol,1}, n::Int64)
    quadPoints, quadWeights=getQuad(2*n-1);
    sk=length(quadPoints);

    nquadPhi=Dict{Symbol, Array{Float64,4}}();
    nquadPoints=Array{Array{Float64,2},1}(undef, size(normals,2));

    for i in 1:size(normals,2)
        n=normals[:,i];
        if n==[0.0,-1.0]
            newQuadPoints=zeros(2,sk);
            newQuadPoints[1,:]=quadPoints;
        elseif n==[-1.0,0.0]
            newQuadPoints=zeros(2,sk);
            newQuadPoints[2,:]=quadPoints;
        elseif n==[0.0, 1.0]
            newQuadPoints=ones(2,sk);
            newQuadPoints[1,:]=quadPoints;
        elseif n==[1.0,0.0]
            newQuadPoints=ones(2,sk);
            newQuadPoints[2,:]=quadPoints;
        elseif n==[0.7071067811865475244, 0.7071067811865475244]
            newQuadPoints=zeros(2,sk);
            newQuadPoints[2,:]=quadPoints;
            newQuadPoints[1,:]=1 .- quadPoints;
        end
        nquadPoints[i]=newQuadPoints;
    end

    for k in type
        phi, psize =getPhi(k);
        nquadPhi[k]=Array{Float64,4}(undef,psize[1],sk,psize[2],size(normals,2));
        for m in 1:size(normals,2)
            for l in 1:psize[2]
                for n in 1:psize[1]
                    for i=1:sk
                        nquadPhi[k][n,i,l,m]=phi[l,n](nquadPoints[m][1,i], nquadPoints[m][2,i]);
                    end
                end
            end
        end
    end
    return nquadPhi, nquadPoints;
end

function splitExplicit(y0::solution,Y::Array{solution,1},FY::Array{solution,1},SthY::Array{SparseMatrixCSC{Float64,Int64},1},
                      p::femProblem, gamma::Float64,
                      nquadPhi::Dict{Symbol, Array{Float64,4}},
                      nquadPoints::Array{Array{Float64,2},1},MrT::SparseMatrixCSC{Float64,Int64}, MrV::SparseMatrixCSC{Float64,Int64},
                      MIS::MIS,t0::Float64,dt::Float64,ns::Int64)
  stage=MIS.nStage;
  #Y=Array{solution,1}(undef,stage+1);
  #FY=Array{solution,1}(undef,stage);
  #SthY=Array{SparseMatrixCSC{Float64,Int64},1}(undef,stage);

  numRho=p.degFBoundary[p.femType[:rho][1]].num
  numRhoV=p.degFBoundary[p.femType[:rhoV][1]].num
  numRhoTheta=p.degFBoundary[p.femType[:rhoTheta][1]].num

  for i in 1:(stage+1)
    Y[i]=y0;
  end

  for i in 1:stage
    for j in 1:i
      Y[i+1]=MIS.alpha[i+1,j]*(Y[j]-y0)+Y[i+1];
    end

    FY[i], SthY[i]=advection(p,gamma,Y[i],nquadPoints,nquadPhi,MrT,MrV);


    fSlow=createSolution(length(y0.rho),length(y0.rhoV),length(y0.rhoTheta),length(y0.v),length(y0.theta));
    SthSlow=spzeros(size(SthY[i],1),size(SthY[i],2))

    for j in 1:i
      fSlow=fSlow+MIS.beta[i+1,j]*FY[j]+(MIS.gamma[i+1,j]/dt)*(Y[j]-y0);
      SthSlow=SthSlow+MIS.beta[i+1,j]*SthY[j];
    end

    nsLoc=ceil(Int64,ns*MIS.d[i+1]);
    dtLoc=dt*MIS.d[i+1];
    dtauLoc=dtLoc/nsLoc;

    Y[i+1].rhoV[numRhoV+1:end]=y0.rhoV[numRhoV+1:end];
    Y[i+1].rho[numRho+1:end]=y0.rho[numRho+1:end];
    Y[i+1].rhoTheta[numRhoTheta+1:end]=y0.rhoTheta[numRhoTheta+1:end];

    symplektischerEuler!(Y[i+1],p,fSlow,SthSlow,nsLoc,dtauLoc);
  end
  return Y[stage+1];
end

function advection(p::femProblem, gamma::Float64, y::solution,
                  nquadPoints::Array{Array{Float64,2},1}, nquadPhi::Dict{Symbol, Array{Float64,4}},
                  MrT::SparseMatrixCSC{Float64,Int64}, MrV::SparseMatrixCSC{Float64,Int64})

  fTtheta=p.femType[:rhoTheta]; fTv=p.femType[:rhoV];
  Fv=@views p.massM[fTv[1]];
  nRhoV=p.degFBoundary[p.femType[:rhoV][1]].num;
  nRhoTheta=p.degFBoundary[p.femType[:rhoTheta][1]].num;

  if p.advection
    cR=projectRhoChi(p,y.rho,y.rhoV,:rho,:rhoV,MrV);
    if p.taskRecovery
      cR=recovery(p,fTv,cR);
      Sv=advectionStiff(p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],
                        p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],sparse(y.rhoV),
                        p.degFBoundary[fTv[3]],nquadPhi[fTv[3]],cR,
                        gamma,p.mesh,p.kubPoints,p.kubWeights,
                        nquadPoints,p.edgeData);
      rCv=Fv\Sv;
    else
      Sv=advectionStiff(p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],
                        p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],sparse(y.rhoV),
                        p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],cR,
                        gamma,p.mesh,p.kubPoints,p.kubWeights,
                        nquadPoints,p.edgeData);
      rCv=Fv\Sv;
    end
    cR=projectRhoChi(p,y.rho,y.rhoTheta,:rho,:rhoTheta,MrT);
    if p.taskRecovery
      cR=recovery(p,fTtheta,cR,:theta);
      Sth=advectionStiffMatrix(p.degFBoundary[fTtheta[1]],nquadPhi[fTtheta[1]],
                         p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],y.rhoV,
                         p.degFBoundary[fTtheta[3]],nquadPhi[fTtheta[3]],cR,
                         gamma,p.mesh,p.kubPoints,p.kubWeights,
                         nquadPoints, p.edgeData);

    else
      Sth=advectionStiffMatrix(p.degFBoundary[fTtheta[1]],nquadPhi[fTtheta[1]],
                         p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],y.rhoV,
                         p.degFBoundary[fTtheta[1]],nquadPhi[fTtheta[1]],cR,
                         gamma,p.mesh,p.kubPoints,p.kubWeights,
                         nquadPoints, p.edgeData);
    end
  else
    rCv=zeros(nRhoV);
    Sth=spzeros(nRhoTheta,length(y.rhoV));
  end

  f=createSolution(length(y.rho),length(y.rhoV),length(y.rhoTheta),length(y.v),length(y.theta));
  f.rhoV[1:nRhoV]=rCv[:,1];

  return f, Sth;
end
