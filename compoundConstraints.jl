include("testCompoundAcoustic.jl")
p=testCompoundAcoustic();

#constraint matrix
#βi,j = j-th beta from sub-element i

#1.
#β1,1=β2,1=1 ,  βi,1=0 i∈{3,..,12}
#βi,2-β1+mod(i,12),3=0 i∈{1,..,12}  (could be βi,2+β1+mod(i,12),3=0)

#2.
#(β1,1+β1,2-β1,3)-(βi,1+βi,2-βi,3)=0 i∈{2,..,12}

#3.
#have a look in integratedAnsatzfct, its full of [0,1,1] vectors
#∑_(j=1)^n 0*βj,1 + 1*βj,2 + 1*βj,3 = 0

#aufbau: zuerst alle 12 beta1, dann alle 12 beta2, dann alle 12 beta3
b=zeros(36);
A=zeros(36,36);

#first constraint
#outer edges
A[1,1]=1.0; b[1]=1.0;
A[2,2]=1.0; b[2]=1.0;
for i in 3:12
    A[i,i]=1.0;
    b[i]=0.0;
end

#inner edges
for i in 1:12
    A[12+i,12+i]=1.0;
#    A[12+i,24+mod(12,i)+1]=-1.0;
    A[12+i,24+mod(12,i)+1]=1.0;
    b[12+i]=0.0;
end

#second constraint
for i in 2:12
    #FEMEuler divergences
#    A[24+i,[1,12+1,24+1]]=[1.0,1.0,-1.0];
#    A[24+i,[i,12+i,24+i]]=[-1.0,-1.0,1.0];
    #MelvinThuburn divergences (common divergences)
    A[24+i,[1,12+1,24+1]]=[1.0,1.0,1.0];
    A[24+i,[i,12+i,24+i]]=[-1.0,-1.0,-1.0];
    b[24+i]=0.0;
end

#third constraint (with a look in integratedAnsatzfct the third constraint
#results in β1,x=0, β2,x=1, β3,x=1 for a HexMesh with one Hexagon in Ω=[0,1]x[0,1])
A[36,13:36].=1.0;

#compute beta-values
betas=A\b;

#third constraint
#compute integrals of every ansatzfunction over its specific cell
#obviously with compound Mesh
#integratedAnsatzfct[k][j][i] : integral value from i-th ansatzfct in subCell j from Cell k
function integrateAnsatzfctOverCells(p::femProblem)
    m=p.mesh;
    nSubCells=p.compoundData.nSubCells;
    kubWeights=p.kubWeights;
    kubPoints=p.kubPoints;
    sk=size(kubWeights);
    degF=p.degFBoundary[:RT0];
    phi=p.degFBoundary[:RT0].phi;
    nAnsatzfct=size(phi,2);
    jphi=initJacobi((m.geometry.dim,nAnsatzfct),sk);
    integratedAnsatzfct=Array{Array{Array{Float64,1},1},1}();
    integratedAnsatzfctCell=Array{Array{Float64,1},1}(undef,nSubCells);
    integratedAnsatzfctSubCell=Array{Float64,1}(undef,nAnsatzfct);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ=Array{Float64,2}(undef,sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,diff(m.topology.offset["20"][1:2])[1]);
    key="20";
    mt=m.meshType;
    subcoord=Array{Array{Float64,2},1}(undef,nSubCells);

    for Cell in 1:p.mesh.topology.size[p.mesh.geometry.dim+1]
        #get coordinates of subCells
        rstart=m.topology.offset[key][Cell];
        rend=m.topology.offset[key][Cell+1]-1;
        z=1;
        for j in rstart:rend
            for i in 1:2
                coord[i,z]=m.geometry.coordinates[i,m.topology.incidence[key][j]];
            end
            z+=1;
        end

        subcoord=p.compoundData.getSubCells[Cell];
        gvertices=l2g(degF,Cell);
        for subCell in 1:nSubCells
            tangent=(subcoord[subCell][:,2].-subcoord[subCell][:,1]);
            tangent=tangent./norm(tangent,2);
            jacobi!(J,ddJ,jphi,kubPoints,phi,subcoord[subCell],mt);
            for i in 1:nAnsatzfct
                for r in 1:sk[2]
                    for l in 1:sk[1]
                        dotp=0.0;
                        for d in 1:m.geometry.dim
                            dotp+=jphi[d,i][l,r]*tangent[d];
                        end
                        integratedAnsatzfctSubCell[i]+=kubWeights[l,r]*(ddJ[l,r]/abs(ddJ[l,r]))*dotp;
                    end
                end
            end
            integratedAnsatzfctCell[subCell]=integratedAnsatzfctSubCell;
        end
        push!(integratedAnsatzfct,integratedAnsatzfctCell);
    end
    return integratedAnsatzfct;
end

integratedAnsatzfct=integrateAnsatzfctOverCells(p);
