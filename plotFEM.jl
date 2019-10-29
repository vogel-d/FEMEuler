#Funktion zum Plotten eines zweidimesnionalen Meshes
#Input: p ist das FEM_Problem, in dem das zu plottende Mesh gespeichert ist,
#       showann gibt an, ob die Nummerierung angezeigt werden soll.
#Die restlichen Argumente sind keyword-Arguments, d.h. sie müssen, wenn sie
#geändert werden sollen, angegeben werden mit karg=wert (z.B.: schowvertices=false)
#Keyword Arguments: showver gibt an, ob die Knoten markiert werden sollen,
#                   showdegf gibt an, ob die Freiheitsgrade mit der Nummerierung angezeigt werden sollen,
#                   sizevertices bestimmt die Größe der Knotenpunkte,
#                   linecolor bestimmt die Farbe der Linien,
#                   colorann/positionann/sizeann ist jeweils ein Vektor, in dem
#                   die Farbe/Position/Größe der Nummerierung festgelegt wird,
#                   wobei die Anordnung [vertices, edges, facets,degrees of freedom] gilt.

function plotFEM(p::femProblem, showann::Bool=true; showboundary::Bool=false,keys::Array{Symbol,1}=collect(keys(p.degF)), showvertices::Bool=showann, showdegf::Bool=showann,
    sizevertices::Float64=2.5, linecolor::Symbol=:black,
    colorann::Array{Symbol,1}=[:grey, :darkgrey, :black],
    positionann::Array{Symbol,1}=[:bottom, :bottom, :bottom,:top],
    sizeann::Array{Int64,1}=[12,12,12,12])

    m=p.mesh;

    if showboundary
        degF=p.degFBoundary
    else
        degF=p.degF
    end

    ince=m.topology.incidence["10"];

    incf=m.topology.incidence["20"];
    off=m.topology.offset["20"];

    coord=m.geometry.coordinates;
    nv, ne, nf=m.topology.size[1:3];
    x=Array{Float64,2}(undef,2,ne);
    y=Array{Float64,2}(undef,2,ne);

    coorde=Array{Float64,2}(undef,2,ne);
    for k in 1:ne
        i=ince[[(2*k-1),2*k]];
        x[:,k]=coord[1,i];
        y[:,k]=coord[2,i];
        coorde[1,k]=0.5*sum(x[:,k]);
        coorde[2,k]=0.5*sum(y[:,k]);
    end

    pl=plot(x,y, c=linecolor, legend=false, xlabel="x-Achse", ylabel="z-Achse");
    if showann
        coordf=Array{Float64,2}(undef,2,nf);
        ng=off[2]-off[1];
        for k in 1:nf
            i=incf[off[k]:(off[k+1]-1)];
            coordf[1,k]=sum(coord[1,i])/ng;
            coordf[2,k]=sum(coord[2,i])/ng;
        end

        a=[(coord[1,i], coord[2,i], text("$i", colorann[1], positionann[1], sizeann[1])) for i in 1:nv];
        ae=[(coorde[1,i], coorde[2,i], text("$i", colorann[2], positionann[2], sizeann[2])) for i in 1:ne];
        af=[(coordf[1,i], coordf[2,i], text("$i", colorann[3], positionann[3], sizeann[3])) for i in 1:nf];
        append!(append!(a,ae),af);
        annotate!(a);
    end

    if showvertices
        plot!(coord[1,:], coord[2,:], seriestype=:scatter, c=:black, markersize=sizevertices);
    end

    if showdegf
        nc=0.8/length(keys);
        nbc=0.8/length(keys);
        z=0;
        for k in keys
            z+=1;
            num=degF[k].num;
            c=degF[k].coordinates[:,1:num];
            plot!(c[1,:], c[2,:], seriestype=:scatter, c=:black, markersize=sizevertices);
            #ad=[(c[1,i], c[2,i], text("$i", ColorGradient(:magma)[z*nc], positionann[4], sizeann[4])) for i in 1:size(c,2)];
            #annotate!(ad);

            bc=degF[k].coordinates[:,(num+1):end];
            if !isempty(bc)
                plot!(bc[1,:], bc[2,:], seriestype=:scatter, c=:black, markersize=sizevertices);
                #ad=[(bc[1,i], bc[2,i], text("$i", ColorGradient(:viridis)[(length(keys)-z)*nbc], positionann[4], sizeann[4])) for i in 1:size(bc,2)];
                #annotate!(ad);
            end

        end
    end

    return pl
end
