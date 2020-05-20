function generateCubedSphere(n::Int,r::Float64,nz::Int=0,case::Symbol=:purser1) #nz rausnehmen?

  piFourth=atan(1.0)
  d=2.0*piFourth/n
  dd=2.0/n;
  NumberOfNodesPlane=(6*(n-1)*(n-1)+12*(n-1)+8)
  NumberOfNodes=(nz+1)*(6*(n-1)*(n-1)+12*(n-1)+8)


  #Definieren der Koordinaten
  #coord=Array{Float64,2}(undef,3,NumberOfNodes);
  coord=Array{Float64,2}(undef,3,NumberOfNodesPlane);
  NodeNumber=1
  #Faces
  # West
  NodeNumberW=NodeNumber
  NodeNumber=getCubePoints!(coord,[-1.0,0.0,0.0],3,2,NodeNumber,n,case)
  # East
  NodeNumberE=NodeNumber
  NodeNumber=getCubePoints!(coord,[1.0,0.0,0.0],3,2,NodeNumber,n,case)
  # South
  NodeNumberS=NodeNumber
  NodeNumber=getCubePoints!(coord,[0.0,-1.0,0.0],3,1,NodeNumber,n,case)
  # North
  NodeNumberN=NodeNumber
  NodeNumber=getCubePoints!(coord,[0.0,1.0,0.0],3,1,NodeNumber,n,case)
  # Bottom
  NodeNumberB=NodeNumber
  NodeNumber=getCubePoints!(coord,[0.0,0.0,-1.0],2,1,NodeNumber,n,case)
  # Top
  NodeNumberT=NodeNumber
  NodeNumber=getCubePoints!(coord,[0.0,0.0,1.0],2,1,NodeNumber,n,case)

  #Edges
  #West East
  NodeNumberWEmm=NodeNumber
  NodeNumber=getCubePoints!(coord,[0.0,-1.0,-1.0],1,NodeNumber,n,case)
  NodeNumberWEpm=NodeNumber
  NodeNumber=getCubePoints!(coord,[0.0,1.0,-1.0],1,NodeNumber,n,case)
  NodeNumberWEmp=NodeNumber
  NodeNumber=getCubePoints!(coord,[0.0,-1.0,1.0],1,NodeNumber,n,case)
  NodeNumberWEpp=NodeNumber
  NodeNumber=getCubePoints!(coord,[0.0,1.0,1.0],1,NodeNumber,n,case)
  #South North
  NodeNumberSNmm=NodeNumber
  NodeNumber=getCubePoints!(coord,[-1.0,0.0,-1.0],2,NodeNumber,n,case)
  NodeNumberSNpm=NodeNumber
  NodeNumber=getCubePoints!(coord,[1.0,0.0,-1.0],2,NodeNumber,n,case)
  NodeNumberSNmp=NodeNumber
  NodeNumber=getCubePoints!(coord,[-1.0,0.0,1.0],2,NodeNumber,n,case)
  NodeNumberSNpp=NodeNumber
  NodeNumber=getCubePoints!(coord,[1.0,0.0,1.0],2,NodeNumber,n,case)
  #Bottom Top
  NodeNumberBTmm=NodeNumber
  NodeNumber=getCubePoints!(coord,[-1.0,-1.0,0.0],3,NodeNumber,n,case)
  NodeNumberBTpm=NodeNumber
  NodeNumber=getCubePoints!(coord,[1.0,-1.0,0.0],3,NodeNumber,n,case)
  NodeNumberBTmp=NodeNumber
  NodeNumber=getCubePoints!(coord,[-1.0,1.0,0.0],3,NodeNumber,n,case)
  NodeNumberBTpp=NodeNumber
  NodeNumber=getCubePoints!(coord,[1.0,1.0,0.0],3,NodeNumber,n,case)

  #Nodes
  NodeNumbermmm=NodeNumber
  coord[:,NodeNumber]=cubePoint([-1.0,-1.0,-1.0],[0,0,0],n,case)
  NodeNumber+=1
  NodeNumberpmm=NodeNumber
  coord[:,NodeNumber]=cubePoint([1.0,-1.0,-1.0],[0,0,0],n,case)
  NodeNumber+=1
  NodeNumbermpm=NodeNumber
  coord[:,NodeNumber]=cubePoint([-1.0,1.0,-1.0],[0,0,0],n,case)
  NodeNumber+=1
  NodeNumberppm=NodeNumber
  coord[:,NodeNumber]=cubePoint([1.0,1.0,-1.0],[0,0,0],n,case)
  NodeNumber+=1
  NodeNumbermmp=NodeNumber
  coord[:,NodeNumber]=cubePoint([-1.0,-1.0,1.0],[0,0,0],n,case)
  NodeNumber+=1
  NodeNumberpmp=NodeNumber
  coord[:,NodeNumber]=cubePoint([1.0,-1.0,1.0],[0,0,0],n,case)
  NodeNumber+=1
  NodeNumbermpp=NodeNumber
  coord[:,NodeNumber]=cubePoint([-1.0,1.0,1.0],[0,0,0],n,case)
  NodeNumber+=1
  NodeNumberppp=NodeNumber
  coord[:,NodeNumber]=cubePoint([1.0,1.0,1.0],[0,0,0],n,case)
  NodeNumber+=1

  NumberOfEdges=(nz+1)*(12*(n-1)*n+12*n)+nz*(6*(n-1)*(n-1)+12*(n-1)+8)
  NumberOfEdgesPlane=12*(n-1)*n+12*n
  #Inzidenz 1->0
  ince=Int[];
  EdgeNumber=1
  # Faces
  # West
  EdgeNumber, EdgeNumberW1,EdgeNumberW2=insertFaceEdge!(ince,NodeNumberW,NodeNumberBTmm,NodeNumberBTmp,NodeNumberSNmm,NodeNumberSNmp,EdgeNumber, n)
  # East
  EdgeNumber, EdgeNumberE1,EdgeNumberE2=insertFaceEdge!(ince,NodeNumberE,NodeNumberBTpm,NodeNumberBTpp,NodeNumberSNpm,NodeNumberSNpp,EdgeNumber, n)
  # South
  EdgeNumber, EdgeNumberS1,EdgeNumberS2=insertFaceEdge!(ince,NodeNumberS,NodeNumberBTmm,NodeNumberBTpm,NodeNumberWEmm,NodeNumberWEmp,EdgeNumber, n)
  # North
  EdgeNumber, EdgeNumberN1,EdgeNumberN2=insertFaceEdge!(ince,NodeNumberN,NodeNumberBTmp,NodeNumberBTpp,NodeNumberWEpm,NodeNumberWEpp,EdgeNumber, n)
  # Bottom
  EdgeNumber, EdgeNumberB1,EdgeNumberB2=insertFaceEdge!(ince,NodeNumberB,NodeNumberSNmm,NodeNumberSNpm,NodeNumberWEmm,NodeNumberWEpm,EdgeNumber, n)
  # Top
  EdgeNumber, EdgeNumberT1,EdgeNumberT2=insertFaceEdge!(ince,NodeNumberT,NodeNumberSNmp,NodeNumberSNpp,NodeNumberWEmp,NodeNumberWEpp,EdgeNumber, n)
# Edges
# West East
  EdgeNumber, EdgeNumberWEmm=insertEdgeEdge!(ince,NodeNumberWEmm,NodeNumbermmm,NodeNumberpmm,EdgeNumber, n)
  EdgeNumber, EdgeNumberWEpm=insertEdgeEdge!(ince,NodeNumberWEpm,NodeNumbermpm,NodeNumberppm,EdgeNumber, n)
  EdgeNumber, EdgeNumberWEmp=insertEdgeEdge!(ince,NodeNumberWEmp,NodeNumbermmp,NodeNumberpmp,EdgeNumber, n)
  EdgeNumber, EdgeNumberWEpp=insertEdgeEdge!(ince,NodeNumberWEpp,NodeNumbermpp,NodeNumberppp,EdgeNumber, n)
# South North
  EdgeNumber, EdgeNumberSNmm=insertEdgeEdge!(ince,NodeNumberSNmm,NodeNumbermmm,NodeNumbermpm,EdgeNumber, n)
  EdgeNumber, EdgeNumberSNpm=insertEdgeEdge!(ince,NodeNumberSNpm,NodeNumberpmm,NodeNumberppm,EdgeNumber, n)
  EdgeNumber, EdgeNumberSNmp=insertEdgeEdge!(ince,NodeNumberSNmp,NodeNumbermmp,NodeNumbermpp,EdgeNumber, n)
  EdgeNumber, EdgeNumberSNpp=insertEdgeEdge!(ince,NodeNumberSNpp,NodeNumberpmp,NodeNumberppp,EdgeNumber, n)
# Bottom Top
  EdgeNumber, EdgeNumberBTmm=insertEdgeEdge!(ince,NodeNumberBTmm,NodeNumbermmm,NodeNumbermmp,EdgeNumber, n)
  EdgeNumber, EdgeNumberBTpm=insertEdgeEdge!(ince,NodeNumberBTpm,NodeNumberpmm,NodeNumberpmp,EdgeNumber, n)
  EdgeNumber, EdgeNumberBTmp=insertEdgeEdge!(ince,NodeNumberBTmp,NodeNumbermpm,NodeNumbermpp,EdgeNumber, n)
  EdgeNumber, EdgeNumberBTpp=insertEdgeEdge!(ince,NodeNumberBTpp,NodeNumberppm,NodeNumberppp,EdgeNumber, n)

  NumberOfFaces=(nz+1)*6*n*n+nz*NumberOfEdgesPlane
  NumberOfFacesPlane=6*n*n
  #Inzidenz 2->1
  incfe=Int[];

  # Faces
  # West
  insertFaceFace!(incfe,EdgeNumberW1,EdgeNumberW2,EdgeNumberSNmm,EdgeNumberSNmp,EdgeNumberBTmm,EdgeNumberBTmp,n)
  # East
  insertFaceFace!(incfe,EdgeNumberE1,EdgeNumberE2,EdgeNumberSNpm,EdgeNumberSNpp,EdgeNumberBTpm,EdgeNumberBTpp,n)
  # South
  insertFaceFace!(incfe,EdgeNumberS1,EdgeNumberS2,EdgeNumberWEmm,EdgeNumberWEmp,EdgeNumberBTmm,EdgeNumberBTpm,n)
  # North
  insertFaceFace!(incfe,EdgeNumberN1,EdgeNumberN2,EdgeNumberWEpm,EdgeNumberWEpp,EdgeNumberBTmp,EdgeNumberBTpp,n)
  # Bottom
  insertFaceFace!(incfe,EdgeNumberB1,EdgeNumberB2,EdgeNumberWEmm,EdgeNumberWEpm,EdgeNumberSNmm,EdgeNumberSNpm,n)
  # Top
  insertFaceFace!(incfe,EdgeNumberT1,EdgeNumberT2,EdgeNumberWEmp,EdgeNumberWEpp,EdgeNumberSNmp,EdgeNumberSNpp,n)
  #Initialisieren des Offsets mit den Einträgen "21" und "10"
  off=Dict("21"=>collect(1:4:(4*NumberOfFacesPlane+1)),"10"=>collect(1:2:(2*NumberOfEdgesPlane+1)));

  #Initialisieren der Inzidenz mit den Einträgen "21" und "10"
  inc=Dict("21"=>incfe,"10"=>ince);

  size=Int[NumberOfNodesPlane, NumberOfEdgesPlane, NumberOfFacesPlane]

  bE=spzeros(Int, size[2]);
  bV=spzeros(Int, size[1]);

  coord=(r/sqrt(3)).*coord;

  #Initialisieren der Topologie, Geometrie und damit des Meshes
  n=Int[n,n,nz];
  l=Float64[-r,-r,-r];
  r=Float64[r,r,r];
  mT=meshTopology(inc,off,size,2,n);
  mG=meshGeometry(coord,l,r);
  m=mesh(mT,mG,bE,bV,4);

  inc,off=faceVertices(m);
  m.topology.incidence["20"]=inc;
  m.topology.offset["20"]=off;

  setOrientation!(m)

  # Extension to 3 d
  #Extension3DPolyGrid(nz,PolyGrid,NumberOfNodesPlane,NumberOfEdgesPlane,NumberOfFacesPlane)
  return m;
end

function getCubePoints!(coord::Array{Float64,2},P::Array{Float64,1},pk::Int,pj::Int,NodeNumber::Int,n::Int,case::Symbol)
  ek=zeros(Int,3); ek[pk]=1;
  ej=zeros(Int,3); ej[pj]=1;
  for k in 1:n-1
    for j in 1:n-1
      coord[:,NodeNumber]=cubePoint(P,j.*ej.+k.*ek,n,case)
      NodeNumber+=1
    end
  end
  return NodeNumber;
end

function getCubePoints!(coord::Array{Float64,2},P::Array{Float64,1},pi::Int,NodeNumber::Int,n::Int,case::Symbol)
  ei=zeros(Int,3); ei[pi]=1;
  for i in 1:n-1
    coord[:,NodeNumber]=cubePoint(P,i.*ei,n,case)
    NodeNumber+=1
  end
  return NodeNumber;
end

function cubePoint(p::Array{Float64,1},i::Array{Int,1},n::Int,case::Symbol)
  if case==:cube1
    N=copy(p)
    if i[1]>0
      N[1]=tan(i[1]*pi/(2*n)-0.25*pi);
    end
    if i[2]>0
      N[2]=tan(i[2]*pi/(2*n)-0.25*pi);
    end
    if i[3]>0
      N[3]=tan(i[3]*pi/(2*n)-0.25*pi);
    end
    return N./norm(N).*sqrt(3);

  elseif case==:purser1
    xm=zeros(2);
    dd=2.0/n;
    if (iszero(i[1]))
        if p[1]==1.0
          panel=1; fac=1.0;
        else
          panel=3; fac=-1.0;
        end
        !iszero(i[2]) ? xm[1]=fac*(-1.0+i[2]*dd) : xm[1]=fac*p[2];
        !iszero(i[3]) ? xm[2]=-1.0+i[3]*dd : xm[2]=p[3];
    elseif (iszero(i[2]))
        if p[2]==1.0
          panel=2; fac=-1.0;
        else
          panel=4; fac=1.0;
        end
        !iszero(i[1]) ? xm[1]=fac*(-1.0+i[1]*dd) : xm[1]=p[1];
        !iszero(i[3]) ? xm[2]=-1.0+i[3]*dd :  xm[2]=p[3];
    elseif (iszero(i[3]))
        if p[3]==1.0
          panel=5; fac=-1.0;
        else
          panel=6; fac=1.0;
        end
        !iszero(i[1]) ? xm[2]=fac*(-1.0+i[1]*dd) : xm[2]=p[1];
        !iszero(i[2]) ? xm[1]=-1.0+i[2]*dd : xm[1]=p[2];
      end
      return sqrt(3.0).*xmtoxc(xm,panel)
  end
end

function faceVertices(m::mesh)
  inc=m.topology.incidence["21"];
  off=m.topology.offset["21"];
  ince=m.topology.incidence["10"];
  offe=m.topology.offset["10"];
  coord=m.geometry.coordinates;
  nf=m.topology.size[3];
  incs=Int[];
  for k in 1:nf
    edges=inc[off[k]:off[k+1]-1];
    inck=Int[];
    ind=Set{Int}(edges[2:4])
    append!(inck,ince[offe[edges[1]]:offe[edges[1]+1]-1])
    while length(inck)<4
      vm=inck[end];
      for j in ind
        vh=ince[offe[j]:offe[j+1]-1];
        if vh[1]==vm
          push!(inck,vh[2]);
          setdiff!(ind,j);
          break;
        elseif vh[2]==vm
          push!(inck,vh[1]);
          setdiff!(ind,j);
          break;
        else
          continue;
        end
      end
    end
    append!(incs,inck);
  end
  offs=collect(1:4:(4*m.topology.size[3]+1))
  return incs, offs
end
