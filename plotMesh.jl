#Funktion zum Plotten eines zweidimesnionalen Meshes
#Input: m ist das zu plottende Mesh,
#       showann gibt an, ob die Nummerierung angezeigt werden soll.
#Die restlichen Argumente sind keyword-Arguments, d.h. sie müssen, wenn sie
#geändert werden sollen, angegeben werden mit karg=wert (z.B.: schowvertices=false)
#Keyword Arguments: showver gibt an, ob die Knoten markiert werden sollen,
#                   sizevertices bestimmt die Größe der Knotenpunkte,
#                   linecolor bestimmt die Farbe der Linien,
#                   colorann/positionann/sizeann ist jeweils ein Verktor, in dem
#                   die Farbe/Position/Größe der Nummerierung festgelegt wird,
#                   wobei die Anordnung [vertices, edges, facets] gilt.
function plotMesh2D(m::mesh, showann::Bool=true; showvertices::Bool=showann,
    sizevertices::AbstractFloat=2.5, linecolor::Symbol=:blue,
    colorann::Array{Symbol,1}=[:grey, :darkgrey, :black],
    positionann::Array{Symbol,1}=[:bottom, :bottom, :auto],
    sizeann::Array{Int,1}=[12,12,12])

    ince=m.topology.incidence["10"];

    incf=m.topology.incidence["20"];
    off=m.topology.offset["20"];

    coord=m.geometry.coordinates;
    nv, ne, nf=m.topology.size[1:3];
    x=Array{AbstractFloat,2}(undef,2,ne);
    y=Array{AbstractFloat,2}(undef,2,ne);
    coorde=Array{AbstractFloat,2}(undef,2,ne);

    for k in 1:ne
        i=ince[[(2*k-1),2*k]];
        x[:,k]=coord[1,i];
        y[:,k]=coord[2,i];
        coorde[1,k]=0.5*sum(x[:,k]);
        coorde[2,k]=0.5*sum(y[:,k]);
    end

    coordf=Array{AbstractFloat,2}(undef,2,nf);
    ng=off[2]-off[1];

    for k in 1:nf
        i=incf[off[k]:(off[k+1]-1)];
        coordf[1,k]=sum(coord[1,i])/ng;
        coordf[2,k]=sum(coord[2,i])/ng;
    end

    p=plot(x,y, c=linecolor, linewidth=0.5, legend=false, xlabel="x in m", ylabel="z in m");
    if showann
        a=[(coord[1,i], coord[2,i], text("$i", colorann[1], positionann[1], sizeann[1])) for i in 1:nv];
        ae=[(coorde[1,i], coorde[2,i], text("$i", colorann[2], positionann[2], sizeann[2])) for i in 1:ne];
        af=[(coordf[1,i], coordf[2,i], text("$i", colorann[3], positionann[3], sizeann[3])) for i in 1:nf];
        append!(append!(a,ae),af);
        annotate!(a);
    end

    if showvertices
        plot!(coord[1,:], coord[2,:], seriestype=:scatter, c=:black, markersize=sizevertices);
    end

    return p
end


function plotMesh3D(m::mesh, showann::Bool=false; showvertices::Bool=showann,
    sizevertices::AbstractFloat=2.5, linecolor::Symbol=:blue,
    colorann::Array{Symbol,1}=[:grey, :darkgrey, :black],
    positionann::Array{Symbol,1}=[:bottom, :bottom, :auto],
    sizeann::Array{Int,1}=[12,12,12])

    ince=m.topology.incidence["10"];

    #=
    incf=m.topology.incidence["20"];
    offf=m.topology.offset["20"];

    incc=m.topology.incidence["30"];
    offc=m.topology.offset["30"];
    =#

    coord=m.geometry.coordinates;
    #nv, ne, nf, nc=m.topology.size[1:4];
    ne=m.topology.size[2];
    x=Array{AbstractFloat,2}(undef,2,ne);
    y=Array{AbstractFloat,2}(undef,2,ne);
    z=Array{AbstractFloat,2}(undef,2,ne);
    coorde=Array{AbstractFloat,2}(undef,3,ne);
    for k in 1:ne
        i=ince[[(2*k-1),2*k]];
        x[:,k]=coord[1,i];
        y[:,k]=coord[2,i];
        z[:,k]=coord[3,i];
        coorde[1,k]=0.5*sum(x[:,k]);
        coorde[2,k]=0.5*sum(y[:,k]);
        coorde[3,k]=0.5*sum(z[:,k]);
    end
    p=plot(x,y,z, c=linecolor, legend=false);
    #=
    coordf=Array{AbstractFloat,2}(undef,3,nf);
    ng=offf[2]-offf[1];
    for k in 1:nf
        i=incf[offf[k]:(offf[k+1]-1)];
        coordf[1,k]=sum(coord[1,i])/ng;
        coordf[2,k]=sum(coord[2,i])/ng;
        coordf[3,k]=sum(coord[3,i])/ng;
    end
    coordc=Array{AbstractFloat,2}(undef,3,nc);
    ngc=offc[2]-offc[1];
    for k in 1:nc
        i=incc[offc[k]:(offc[k+1]-1)];
        coordc[1,k]=sum(coord[1,i])/ngc;
        coordc[2,k]=sum(coord[2,i])/ngc;
        coordc[3,k]=sum(coord[3,i])/ngc;
    end
    if showann
        a=[(coord[1,i], coord[2,i],coord[3,i], text("$i", colorann[1], positionann[1], sizeann[1])) for i in 1:nv];
        ae=[(coorde[1,i], coorde[2,i],coorde[3,i], text("$i", colorann[2], positionann[2], sizeann[2])) for i in 1:ne];
        af=[(coordf[1,i], coordf[2,i],coordf[3,i], text("$i", colorann[3], positionann[3], sizeann[3])) for i in 1:nf];
        ac=[(coordc[1,i], coordc[2,i],coordc[3,i], text("$i", colorann[3], positionann[3], sizeann[3])) for i in 1:nc];
        append!(append!(append!(a,ae),af),ac);
        annotate!(a);
    end
    =#
    if showvertices
        plot!(coord[1,:], coord[2,:],coord[3,:], seriestype=:scatter, c=:black, markersize=sizevertices);
    end
    return p
end
