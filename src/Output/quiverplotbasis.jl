function quiverPlotCompound(method::Symbol)
    if method==:HexToKites
        m=generateHexMesh(0.0,1.0,0.0,1.0,1,:periodic,:constant,meshType=4); #(east/west, top/bottom)
        compoundPhi=3
        #cart=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
        cart=[i/8 for i in 0:8]
        ref_plot_coordx=Float64[];
        ref_plot_coordy=Float64[];
        for i in 1:length(cart)
            for j in 1:length(cart)
                cart[i]==0.0 && continue;
                push!(ref_plot_coordx,cart[i])
                push!(ref_plot_coordy,cart[j])
            end
        end
        ref_plot_coord=zeros(2,length(ref_plot_coordx));
        ref_plot_coord[1,:]=ref_plot_coordx;
        ref_plot_coord[2,:]=ref_plot_coordy;
    elseif method==:HexToTris
        m=generateHexMesh(0.0,1.0,0.0,1.0,1,:periodic,:constant,meshType=3); #(east/west, top/bottom)
        compoundPhi=3
        #ref_plot_coord=[0.1 0.5 0.8 0.1 0.6 0.1 0.4 0.1;
        #                0.1 0.1 0.1 0.3 0.3 0.5 0.5 0.8];
        #cart=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
        cart=[i/8 for i in 0:8]
        ref_plot_coordx=Float64[];
        ref_plot_coordy=Float64[];
        for i in 1:length(cart)
            for j in 1:(length(cart)-i)
                #cart[i]==0.0 && continue;
                push!(ref_plot_coordx,cart[i])
                push!(ref_plot_coordy,cart[j])
            end
        end
        ref_plot_coord=zeros(2,length(ref_plot_coordx));
        ref_plot_coord[1,:]=ref_plot_coordx;
        ref_plot_coord[2,:]=ref_plot_coordy;
    elseif method==:RectToKites
        m=generateRectMesh(1,1,:periodic,:constant,0.0,1.0,0.0,1.0);
        compoundPhi=2
        cart=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
        ref_plot_coordx=Float64[];
        ref_plot_coordy=Float64[];
        for i in 1:length(cart)
            for j in 1:length(cart)
                push!(ref_plot_coordx,cart[i])
                push!(ref_plot_coordy,cart[j])
            end
        end
        ref_plot_coord=zeros(2,length(ref_plot_coordx));
        ref_plot_coord[1,:]=ref_plot_coordx;
        ref_plot_coord[2,:]=ref_plot_coordy;
    elseif method==:RectToTris
        m=generateRectMesh(1,1,:periodic,:constant,0.0,1.0,0.0,1.0,meshType=3);
        compoundPhi=2
        ref_plot_coord=zeros(2,121);
        cart=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
        for i in 1:length(cart)
            for j in 1:(length(cart)-i)
                ref_plot_coord[:,(i-1)*length(cart)+j]=[cart[i],cart[j]]
            end
        end
    else
        error("unknown method");
    end

    femType=Dict(:p=>[:DG0, :P1, :DG1, :DG0], :v=>[:RT0, :VecP1, :VecDG1, :RT0B], :b=>[:DG0, :P1, :DG1, :DG0]);
    p=femProblem(m, femType, compoundMethod=method);
    nSubCells=p.data.compoundData.nSubCells;
    meshType=m.meshType;

    assemblePhiPre!(p);

    subcoord=Array{Array{Float64,2},1}(undef,nSubCells);
    for i in 1:nSubCells
        #fill! causes mutating all entries of subcoord when changing a single entry
        subcoord[i]=Array{Float64,2}(undef,2,p.data.compoundData.nVerticesSubElement);
    end
    coord= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][1]:m.topology.offset["20"][1+1]-1]]
    center=Array{Float64,1}(undef,2)
    getSubCells!(subcoord, coord, center, p.data.compoundData);

    if meshType==4
        phi, divphi, gradphi, cm, nFace, nEdge, nVert=getQuadElementProperties(:RT0);
        #ref_plot_coord=[0.1 0.3 0.5 0.7 0.9 0.1 0.3 0.5 0.7 0.9 0.1 0.3 0.5 0.7 0.9;
        #                0.1 0.1 0.1 0.3 0.3 0.3 0.5 0.5 0.5 0.7 0.7 0.7 0.9 0.9 0.9];

    elseif meshType==3
        phi, divphi, gradphi, cm, nFace, nEdge, nVert=getTriElementProperties(:RT0);

    else
        error("unknown meshType")
    end

    nrefcoords=size(ref_plot_coord,2);

    quiverplot_coord=Array{Float64,2}(undef,2,nrefcoords*nSubCells);
    quiverplotx=Array{Float64,1}(undef,nrefcoords*nSubCells);
    quiverploty=Array{Float64,1}(undef,nrefcoords*nSubCells);

    plot_coord=similar(ref_plot_coord);


    for subCell in 1:nSubCells
        c_now=subcoord[subCell];
        if meshType==3
            for i in 1:nrefcoords
                ref_point=ref_plot_coord[:,i];
                quiverplot_coord[:,(subCell-1)*nrefcoords+i]=c_now[:,1]+
                                    (c_now[:,2]-c_now[:,1])*ref_point[1]+
                                    (c_now[:,3]-c_now[:,1])*ref_point[2];
                J11=c_now[1,2]-c_now[1,1];
                J21=c_now[2,2]-c_now[2,1];
                J12=c_now[1,3]-c_now[1,1];
                J22=c_now[2,3]-c_now[2,1];
                J=[J11 J12;J21 J22];
                dJ=J11*J22-J21*J12;
                vec=zeros(2);
                for subPhi in 1:size(phi,2)
                    vec+=p.data.compoundData.assembledPhiPre[:RT0][1][compoundPhi][subPhi,subCell]*
                        (1/dJ)*J*[phi[1,subPhi](ref_point[1],ref_point[2]);phi[2,subPhi](ref_point[1],ref_point[2])];
                end
                quiverplotx[(subCell-1)*nrefcoords+i]=vec[1];
                quiverploty[(subCell-1)*nrefcoords+i]=vec[2];
            end
        elseif meshType==4
            J=Array{Float64,2}(undef,2,2);
            for i in 1:nrefcoords
                ref_point=ref_plot_coord[:,i];
                if ref_point[1]==0.0 && method==:HexToKites
                    #continue;
                end
                quiverplot_coord[:,(subCell-1)*nrefcoords+i]=c_now[:,1]+
                                    (c_now[:,2]-c_now[:,1])*ref_point[1]+
                                    (c_now[:,4]-c_now[:,1])*ref_point[2]+
                                    (c_now[:,3]-c_now[:,4]-c_now[:,2]+c_now[:,1])*
                                    ref_point[1]*ref_point[2];

                J[1,1]=(c_now[1,2]-c_now[1,1])+(c_now[1,3]-c_now[1,4]-c_now[1,2]+c_now[1,1])*ref_point[2];
                J[2,1]=(c_now[2,2]-c_now[2,1])+(c_now[2,3]-c_now[2,4]-c_now[2,2]+c_now[2,1])*ref_point[2];
                J[1,2]=(c_now[1,4]-c_now[1,1])+(c_now[1,3]-c_now[1,4]-c_now[1,2]+c_now[1,1])*ref_point[1];
                J[2,2]=(c_now[2,4]-c_now[2,1])+(c_now[2,3]-c_now[2,4]-c_now[2,2]+c_now[2,1])*ref_point[1];
                dJ=J[1,1]*J[2,2]-J[2,1]*J[1,2];
                vec=zeros(2);
                for subPhi in 1:size(phi,2)
                    vec+=p.data.compoundData.assembledPhiPre[:RT0][1][compoundPhi][subPhi,subCell]*
                        (1/dJ)*J*[phi[1,subPhi](ref_point[1],ref_point[2]);phi[2,subPhi](ref_point[1],ref_point[2])];
                end
                quiverplotx[(subCell-1)*nrefcoords+i]=vec[1];
                quiverploty[(subCell-1)*nrefcoords+i]=vec[2];
            end
        end
    end

    compress=1/30;
    quiverplotx=compress*quiverplotx;
    quiverploty=compress*quiverploty;

    plot()
    plot!(coord[1,:],coord[2,:],color=:black,legend=:false);
    plot!([coord[1,end], coord[1,1]],[coord[2,end], coord[2,1]],color=:black,legend=:false);
    for subCell in 1:nSubCells
        if meshType==4
            plot!([subcoord[subCell][1,1],subcoord[subCell][1,4]],[subcoord[subCell][2,1],subcoord[subCell][2,4]],color=:black,legend=false);
        elseif meshType==3
            plot!([subcoord[subCell][1,1],subcoord[subCell][1,3]],[subcoord[subCell][2,1],subcoord[subCell][2,3]],color=:black,legend=false);
        end
    end
    quiver!(quiverplot_coord[1,:],quiverplot_coord[2,:],quiver=(quiverplotx,quiverploty),color=:black);
    plot!([0.0],[0.0],legend=false,border=:none)
    display(plot!());
    return nothing;
end


function quiverPlot_simple(method::Symbol)
    if method==:rect
        m=generateRectMesh(1,1,:periodic,:constant,0.0,1.0,0.0,1.0); #(east/west, top/bottom)
        show_phi=2
        cart=[i/20 for i in 0:20];
        ref_plot_coordx=Float64[];
        ref_plot_coordy=Float64[];
        for i in 1:length(cart)
            for j in 1:length(cart)
                cart[i]==0.0 && continue;
                push!(ref_plot_coordx,cart[i])
                push!(ref_plot_coordy,cart[j])
            end
        end
        ref_plot_coord=zeros(2,length(ref_plot_coordx));
        ref_plot_coord[1,:]=ref_plot_coordx;
        ref_plot_coord[2,:]=ref_plot_coordy;
    elseif method==:tri
        m=generateTriMeshEquilateral(0.0,1.0,0.0,1.0,1,:periodic,:constant); #(east/west, top/bottom)
        show_phi=2
        #ref_plot_coord=[0.1 0.5 0.8 0.1 0.6 0.1 0.4 0.1;
        #                0.1 0.1 0.1 0.3 0.3 0.5 0.5 0.8];
        cart=[i/20 for i in 0:20];
        ref_plot_coordx=Float64[];
        ref_plot_coordy=Float64[];
        for i in 1:length(cart)
            for j in 1:(length(cart)+1-i)
                #cart[i]==0.0 && continue;
                push!(ref_plot_coordx,cart[i])
                push!(ref_plot_coordy,cart[j])
            end
        end
        ref_plot_coord=zeros(2,length(ref_plot_coordx));
        ref_plot_coord[1,:]=ref_plot_coordx;
        ref_plot_coord[2,:]=ref_plot_coordy;
    else
        error("unknown method");
    end

    femType=Dict(:p=>[:DG0, :P1, :DG1, :DG0], :v=>[:RT0, :VecP1, :VecDG1, :RT0B], :b=>[:DG0, :P1, :DG1, :DG0]);
    p=femProblem(m, femType);
    meshType=m.meshType;

    coord= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][1]:m.topology.offset["20"][1+1]-1]]

    if meshType==4
        phi, divphi, gradphi, cm, nFace, nEdge, nVert=getQuadElementProperties(:RT0);
        #ref_plot_coord=[0.1 0.3 0.5 0.7 0.9 0.1 0.3 0.5 0.7 0.9 0.1 0.3 0.5 0.7 0.9;
        #                0.1 0.1 0.1 0.3 0.3 0.3 0.5 0.5 0.5 0.7 0.7 0.7 0.9 0.9 0.9];

    elseif meshType==3
        phi, divphi, gradphi, cm, nFace, nEdge, nVert=getTriElementProperties(:RT0);

    else
        error("unknown meshType")
    end

    nrefcoords=size(ref_plot_coord,2);

    quiverplot_coord=Array{Float64,2}(undef,2,nrefcoords);
    quiverplotx=Array{Float64,1}(undef,nrefcoords);
    quiverploty=Array{Float64,1}(undef,nrefcoords);

    plot_coord=similar(ref_plot_coord);


    c_now=coord;
    if meshType==3
        for i in 1:nrefcoords
            ref_point=ref_plot_coord[:,i];
            quiverplot_coord[:,i]=c_now[:,1]+
                                (c_now[:,2]-c_now[:,1])*ref_point[1]+
                                (c_now[:,3]-c_now[:,1])*ref_point[2];
            J11=c_now[1,2]-c_now[1,1];
            J21=c_now[2,2]-c_now[2,1];
            J12=c_now[1,3]-c_now[1,1];
            J22=c_now[2,3]-c_now[2,1];
            J=[J11 J12;J21 J22];
            dJ=J11*J22-J21*J12;

            vec=zeros(2);
            vec+=(1/dJ)*J*[phi[1,show_phi](ref_point[1],ref_point[2]);phi[2,show_phi](ref_point[1],ref_point[2])];

            quiverplotx[i]=vec[1];
            quiverploty[i]=vec[2];
        end
    elseif meshType==4
        J=Array{Float64,2}(undef,2,2);
        for i in 1:nrefcoords
            ref_point=ref_plot_coord[:,i];

            quiverplot_coord[:,i]=c_now[:,1]+
                                (c_now[:,2]-c_now[:,1])*ref_point[1]+
                                (c_now[:,4]-c_now[:,1])*ref_point[2]+
                                (c_now[:,3]-c_now[:,4]-c_now[:,2]+c_now[:,1])*
                                ref_point[1]*ref_point[2];

            J[1,1]=(c_now[1,2]-c_now[1,1])+(c_now[1,3]-c_now[1,4]-c_now[1,2]+c_now[1,1])*ref_point[2];
            J[2,1]=(c_now[2,2]-c_now[2,1])+(c_now[2,3]-c_now[2,4]-c_now[2,2]+c_now[2,1])*ref_point[2];
            J[1,2]=(c_now[1,4]-c_now[1,1])+(c_now[1,3]-c_now[1,4]-c_now[1,2]+c_now[1,1])*ref_point[1];
            J[2,2]=(c_now[2,4]-c_now[2,1])+(c_now[2,3]-c_now[2,4]-c_now[2,2]+c_now[2,1])*ref_point[1];
            dJ=J[1,1]*J[2,2]-J[2,1]*J[1,2];

            vec=zeros(2);
            vec+=(1/dJ)*J*[phi[1,show_phi](ref_point[1],ref_point[2]);phi[2,show_phi](ref_point[1],ref_point[2])];

            quiverplotx[i]=vec[1];
            quiverploty[i]=vec[2];
        end
    end

    compress=1/10;
    quiverplotx=compress*quiverplotx;
    quiverploty=compress*quiverploty;

    plot()
    plot!(coord[1,:],coord[2,:],color=:black,legend=false);
    plot!([coord[1,end], coord[1,1]],[coord[2,end], coord[2,1]],color=:black,legend=false);

    quiver!(quiverplot_coord[1,:],quiverplot_coord[2,:],quiver=(quiverplotx,quiverploty),color=:black);
    plot!([0],[0],legend=false,border=:none)
    display(plot!());
    return nothing;
end
