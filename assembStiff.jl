function assembStiff!(p::femProblem)
    if p.type==:boussinesq
        degF=p.degFBoundary;
        pkey=p.femType[:p][1];
        bkey=p.femType[:b][1];
        vkey=p.femType[:v][1];

        z=[0.0, 1.0];

        Spv=assembStiff(degF[pkey], degF[vkey], p.mesh, p.kubWeights, p.kubPoints);
        Svp = copy(-Spv');
        Sbv=assembStiff(degF[bkey], degF[vkey], z, p.mesh, p.kubWeights, p.kubPoints)
        Svb = copy(Sbv');

        p.stiffM[:pv]=Spv[1:degF[pkey].num,:];
        p.stiffM[:bv]=Sbv[1:degF[bkey].num,:];
        p.stiffM[:vp]=Svp[1:degF[vkey].num,:];
        p.stiffM[:vb]=Svb[1:degF[vkey].num,:];

    elseif p.type==:compressible
        degF=p.degFBoundary;
        rhokey=p.femType[:rho][1];
        vkey=p.femType[:rhoV][1];
        pkey=p.femType[:p][1];

        z=[0.0,1.0]

        Spv=assembStiff(degF[pkey], degF[vkey], p.mesh, p.kubWeights, p.kubPoints);
        Svp = copy(-Spv');
        Srhov=assembStiff(degF[rhokey], degF[vkey], z, p.mesh, p.kubWeights, p.kubPoints)
        Svrho = copy(Srhov');

        p.stiffM[:rho]=Spv[1:degF[rhokey].num,:];
        p.stiffM[:vp]=Svp[1:degF[vkey].num,:];
        p.stiffM[:vrho]=Svrho[1:degF[vkey].num,:];

    elseif p.type==:shallow
        degF=p.degFBoundary;
        hkey=p.femType[:h][1];
        vkey=p.femType[:hV][1];

        Spv=assembStiff(degF[pkey], degF[vkey], p.mesh, p.kubWeights, p.kubPoints);
        Svp = copy(-Spv');
        Sfv=assembStiff(degF[vkey],degF[vkey],p.mesh, p.kubWeights, p.kubPoints)

        p.stiffM[:h]=Spv[1:degF[hkey].num,:];
        p.stiffM[:vp]=Svp[1:degF[vkey].num,:];
        p.stiffM[:fv]=Sfv[1:degF[vkey].num,:];
    elseif p.type==:linshallow
        degF=p.degFBoundary;
        hkey=p.femType[:h][1];
        vkey=p.femType[:v][1];

        Spv=assembStiff(degF[hkey], degF[vkey], p.mesh, p.kubWeights, p.kubPoints);
        Svp = copy(-Spv');

        p.stiffM[:div]=Spv[1:degF[hkey].num,:];
        p.stiffM[:grad]=Svp[1:degF[vkey].num,:];
    else
        comp=Set{Symbol}()
        for i in collect(keys(p.femType))
            push!(comp,p.femType[i][1])
        end
        for i in comp
            assembStiff!(p,i);
        end
    end
end

#skalare Größe mit Divergenz von vektorieller Größe
function assembStiff(degFs::degF{1,:H1}, degFv::degF{2,:H1div}, m::mesh, kubWeights::Array{Float64,2}, kubPoints::Array{Float64,2})

    nf=m.topology.size[m.topology.dim+1]
    nT=degFs.numB;
    nF=degFv.numB;
    phiT=degFs.phi;
    divphiF=degFv.divphi;

    rows=Int64[];
    cols=Int64[];
    vals=Float64[];

    # hier ist Assemblieren nur einmal nötig, da sich die Determinanten der Jacobi-Matrix jeweils
    # wegkürzt, weshalb die lokale Steifigkeitsmatrix unabhängig von den Koordinaten der jeweiligen Fläche ist.

    lS=zeros(length(phiT), length(divphiF));
    for i in 1:length(phiT)
        for j in 1:length(divphiF)
            currentval=0.0;
            for k in 1:size(kubWeights,2)
                for l in 1:size(kubWeights,1)
                    currentval+=kubWeights[l,k]*phiT[i][l,k]*divphiF[j][l,k];
                end
            end
            lS[i,j] = currentval;
        end
    end

    sk=size(kubWeights)
    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);
    for k in 1:nf
        globalNumT=l2g(degFs,k);
        globalNumF=l2g(degFv,k);
        jacobi!(J,dJ,m,k,kubPoints,coord);
        for j in 1:length(globalNumF)
            for i in 1:length(globalNumT)
                if !isequal(lS[i,j],0.0) || (globalNumT[i]==nT && globalNumF[j]==nF)
                    push!(rows,globalNumT[i]);
                    push!(cols,globalNumF[j]);
                    push!(vals,(abs(dJ[1])/dJ[1])*lS[i,j]);
                    #ddJ[1] reicht um Vorzeichen der Determinante zu identifizieren
                end
            end
        end
    end
    return sparse(rows,cols,vals);
end

#skalare Größe mit vektorieller Größe
function assembStiff(degFs::degF{1,:H1}, degFv::degF{2,:H1div}, z::Array{Float64,1}, m::mesh, kubWeights::Array{Float64,2}, kubPoints::Array{Float64,2})
    nT=degFs.numB;
    nF=degFv.numB;
    phiT=degFs.phi;
    phiF=degFv.phi;

    rows=Int64[];
    cols=Int64[];
    vals=Float64[];

    sk=size(kubWeights)

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphiF=initJacobi((m.geometry.dim,size(phiF,2)),sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    lS=zeros(length(phiT), size(phiF,2));

    for k in 1:m.topology.size[m.topology.dim+1]
        jacobi!(J,ddJ,jphiF,m, k, kubPoints, phiF,coord);
        for i in 1:length(phiT)
            for j in 1:size(phiF,2)
                currentval=0.0;
                for r in 1:sk[2]
                    for l in 1:sk[1]
                        zjphi=0.0
                        for d in 1:m.geometry.dim
                            zjphi+=z[d]*jphiF[d,j][l,r]
                        end
                        currentval+=kubWeights[l,r]*(ddJ[l,r]/abs(ddJ[l,r]))*phiT[i][l,r]*zjphi;
                    end
                end
                lS[i,j] = currentval;
            end
        end

        globalNumT=l2g(degFs,k);
        globalNumF=l2g(degFv,k);

        for j in 1:length(globalNumF)
            for i in 1:length(globalNumT)
                gi=globalNumT[i];
                gj=globalNumF[j];
                if !isequal(lS[i,j],0.0) || (gi==nT && gj==nF)
                    push!(rows,gi);
                    push!(cols,gj);
                    push!(vals,lS[i,j]);
                end
            end
        end
    end
    return sparse(rows,cols,vals);
end

#vektorielle Größe mit Vektorprodukt von Normale und vektorieller Größe
function assembStiff(degFT::degF{2}, degFF::degF{2}, m::mesh, kubWeights::Array{Float64,2}, kubPoints::Array{Float64,2})
    nT=degFT.numB;
    nF=degFF.numB;
    phiT=degFT.phi;
    phiF=degFF.phi;

    rows=Int64[];
    cols=Int64[];
    vals=Float64[];

    sk=size(kubWeights)

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphiT=initJacobi((m.geometry.dim,size(phiT,2)),sk);
    jphiF=initJacobi((m.geometry.dim,size(phiF,2)),sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    #f=zeros(sk); Omega=7.292e-5;
    #fac=[-1.0, 1.0, 0.0]; perm=[3,1,2]; #Kreuzprodukt mit Normalenvektor der Ebene

    lS=zeros(size(phiT,2), size(phiF,2));

    for k in 1:m.topology.size[m.topology.dim+1]
        jacobi!(J,ddJ,jphiT,jphiF,m, k, kubPoints,phiT,phiF,coord);
        #fill!(f,0.0);
        #for i=1:sk[1], j=1:sk[2]
        #    xyz=transformation(m,coord,kubPoints[1,i],kubPoints[2,j])
        #    f[i,j]=2*Omega*xyz[3]/norm(xyz,2);
        #end
        for i in 1:size(phiT,2)
            for j in 1:size(phiF,2)
                currentval=0.0;
                for r in 1:sk[2]
                    for l in 1:sk[1]
                        #=
                        njphi=0.0;
                        for d in 1:m.geometry.dim
                            njphi+=jphiT[d,i][l,r]*fac[d]*jphiF[perm[d],j][l,r]
                        end
                        currentval+=kubWeights[l,r]*f[l,r]*(ddJ[i,j]^3/abs(ddJ[i,j]))*njphi;
                        =#
                        currentval+=kubWeights[l,r]*2*Omega*abs(ddJ[l,r])*(jphiF[2,j][l,r]*jphiT[1,i][l,r]-jphiF[1,j][l,r]*jphiT[2,i][l,r]);
                    end
                end
                lS[i,j] = currentval;
            end
        end

        globalNumT=l2g(degFT,k);
        globalNumF=l2g(degFF,k);

        for j in 1:length(globalNumF)
            for i in 1:length(globalNumT)
                gi=globalNumT[i];
                gj=globalNumF[j];
                if !isequal(lS[i,j],0.0) || (gi==nT && gj==nF)
                    push!(rows,gi);
                    push!(cols,gj);
                    push!(vals,lS[i,j]);
                end
            end
        end
    end
    return sparse(rows,cols,vals);
end



#Gradient einer skalaren Größe mit dem Gradienten einer skalaren Größe (Poisson-Probleme)
function assembStiff!(p::femProblem, comp::Symbol)
    degF=p.degFBoundary;
    m=p.mesh;
    dphiRef=degF[comp].gradphi;

    kubPoints=p.kubPoints;
    kubWeights=p.kubWeights;
    sk=size(kubWeights)

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    nDF=degF[comp].numB;
    S=spzeros(nDF,nDF);
    lS=zeros(size(dphiRef,2), size(dphiRef,2));

    for k in 1:m.topology.size[m.topology.dim+1]
        jacobi!(J,dJ,m,k,kubPoints,coord);

        for j in 1:size(dphiRef,2)
            for i in 1:size(dphiRef,2)
                for k in 1:sk[2]
                    for l in 1:sk[1]
                        lS[i,j]+=kubWeights[l,k]*dJ[l,k]*(dphiRef[1,i][l,k]*dphiRef[1,j][l,k]+dphiRef[2,i][l,k]*dphiRef[2,j][l,k]);
                    end
                end
            end
        end

        globalNum=l2g(degF[comp],k);
        for j in 1:length(globalNum)
            for i in 1:length(globalNum)
                S[globalNum[i],globalNum[j]]+=lS[i,j];
            end
        end
    end
    p.stiffM[comp]=S[1:degF[comp].num,:];
    return nothing;
end
