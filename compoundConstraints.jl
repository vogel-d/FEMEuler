include("testCompoundAcoustic.jl")
p=testCompoundAcoustic();


#third constraint
#compute integrals of every ansatzfunction over its specific cell
#obviously with compound Mesh
#integratedAnsatzfcn[k][j][i] : integral value from i-th ansatzfcn in subCell j from Cell k
function integrateAnsatzfcnOverCells(p::femProblem)
    m=p.mesh;
    nSubCells=p.compoundData.nSubCells;
    kubWeights=p.kubWeights;
    kubPoints=p.kubPoints;
    sk=size(kubWeights);
    degF=p.degFBoundary[:RT0];
    phi=p.degFBoundary[:RT0].phi;
    nAnsatzfcn=size(phi,2);
    jphi=initJacobi((m.geometry.dim,nAnsatzfcn),sk);
    integratedAnsatzfcn=Array{Array{Array{Float64,1},1},1}();
    integratedAnsatzfcnCell=Array{Array{Float64,1},1}(undef,nSubCells);
    integratedAnsatzfcnSubCell=Array{Float64,1}(undef,nAnsatzfcn);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ=Array{Float64,2}(undef,sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,diff(m.topology.offset["20"][1:2])[1]);
    key="20";
    mt=m.meshType;
    subcoord=Array{Array{Float64,2},1}(undef,nSubCells);


    for Cell in 1:p.mesh.topology.size[p.mesh.geometry.dim+1]
        subcoord=p.compoundData.getSubCells[Cell];
        gvertices=l2g(degF,Cell);
        for subCell in 1:nSubCells
            integratedAnsatzfcnSubCell=zeros(nAnsatzfcn);
            tangent=(subcoord[subCell][:,2].-subcoord[subCell][:,1]);
            tangent=tangent./norm(tangent,2);
            jacobi!(J,ddJ,jphi,kubPoints,phi,subcoord[subCell],mt);
            for i in 1:nAnsatzfcn
                for r in 1:sk[2]
                    for l in 1:sk[1]
                        dotp=0.0;
                        for d in 1:m.geometry.dim
                            dotp+=jphi[d,i][l,r]*tangent[d];
                        end
                        integratedAnsatzfcnSubCell[i]+=kubWeights[l,r]*(ddJ[l,r]/abs(ddJ[l,r]))*dotp;
                    end
                end
            end
            integratedAnsatzfcnCell[subCell]=deepcopy(integratedAnsatzfcnSubCell);
        end
        push!(integratedAnsatzfcn,deepcopy(integratedAnsatzfcnCell));
    end
    return integratedAnsatzfcn;
end

#constraint matrix
#βi,j = j-th beta from sub-element i

#1.
#β1,1=β2,1=1 ,  βi,1=0 i∈{3,..,12}
#βi,2-β1+mod(i,12),3=0 i∈{1,..,12}  (could be βi,2+β1+mod(i,12),3=0)

#2.
#(β1,1+β1,2-β1,3)-(βi,1+βi,2-βi,3)=0 i∈{2,..,12}

#3.
#have a look in integratedAnsatzfcn, its full of [0,1,1] vectors
#∑_(j=1)^n 0*βj,1 + 1*βj,2 + 1*βj,3 = 0

#aufbau: zuerst alle 12 beta1, dann alle 12 beta2, dann alle 12 beta3
b=zeros(36);
A=zeros(36,36);

#first constraint
#outer edges
for i in 1:12
    A[i,i]=1.0;
    b[i]=0.0;
end
b[7]=1.0
b[8]=1.0

#inner edges
for i in 1:12
    A[12+i,12+i]=1.0;
    A[12+i,24+mod(i,12)+1]=-1.0;
#    A[12+i,24+mod(i,12)+1]=1.0;
    b[12+i]=0.0;
end

#second constraint
for i in 1:11
    #FEMEuler divergences
    A[24+i,[1,12+1,24+1]]=[1.0,1.0,-1.0];
    A[24+i,[1+i,12+1+i,24+1+i]]=[-1.0,-1.0,1.0];
    #MelvinThuburn divergences (common divergences)
#    A[24+i,[1,12+1,24+1]]=[1.0,1.0,1.0];
#    A[24+i,[1+i,12+1+i,24+1+i]]=[-1.0,-1.0,-1.0];
    b[24+i]=0.0;
end

#third constraint
coeff=integrateAnsatzfcnOverCells(p);
for subCell in 1:12
    for i in 1:3
        A[36,(i-1)*12+subCell]=coeff[1][subCell][i];
    end
end
#A[36,13:24].=1.0;
#MelvinThuburn ansatzfunctions (common ansatzfunctions)
#A[36,13:24].=1.0;
#A[36,25:36].=-1.0;

#compute beta-values
betas=A\b;


#KITES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#=
bkite=zeros(24);
Akite=zeros(24,24);

#first constraint
#outer edges (β1,i ans β2,i for i∈{1,..,6})
for i in 1:12
    Akite[i,i]=1.0;
    bkite[i]=0.0;
end
bkite[7]=1.0 #first element β2,1
bkite[2]=1.0 #second element β1,2

#inner edges
for i in 1:6
    Akite[12+i,12+i]=1.0;
    Akite[12+i,18+mod(i,6)+1]=-1.0;
#    A[12+i,24+mod(i,6)+1]=1.0;
    bkite[12+i]=0.0;
end

#second constraint
for i in 1:5
    #FEMEuler divergences
    Akite[18+i,[1,6+1,12+1,18+1]]=[1.0,1.0,-1.0,-1.0];
    Akite[18+i,[1+i,6+1+i,12+1+i,18+1+i]]=[-1.0,-1.0,1.0,1.0];
    #MelvinThuburn divergences (common divergences)
#    A[24+i,[1,12+1,24+1]]=[1.0,1.0,1.0];
#    A[24+i,[1+i,12+1+i,24+1+i]]=[-1.0,-1.0,-1.0];
    bkite[18+i]=0.0;
end

#third constraint
coeff=integrateAnsatzfcnOverCells(p);
for subCell in 1:6
    for i in 1:4
        Akite[24,(i-1)*6+subCell]=coeff[1][subCell][i];
    end
end
#A[36,13:24].=1.0;
#MelvinThuburn ansatzfunctions (common ansatzfunctions)
#A[36,13:24].=1.0;
#A[36,25:36].=-1.0;

#compute beta-values
betaskite=Akite\bkite;
=#
