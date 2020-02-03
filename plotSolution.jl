function plotSolution(p::femProblem, key::Symbol, t::Float64, c::Symbol=:plasma)
    nf=p.mesh.topology.size[3];
    sol=getfield(p.solution[t],key);
    smin=minimum(sol);
    smax=maximum(sol);
    off=p.mesh.topology.offset["20"];
    inc=p.mesh.topology.incidence["20"];

    pc=plot();
    for i in 0:0.001:1
        color=ColorGradient(c)[i];
        plot!([i,i],[0,1],linecolor=color, legend=false, yaxis=nothing, xaxis=nothing, foreground_color=:white)
    end

    mw=0.5*(smin+smax);
    mr=0.5*(mw+smax);
    ml=0.5*(smin+mw);
    a=[(0.0,1.0,text("$(round(smin,digits=3))",:bottom, 6)),(0.25,1.0,text("$(round(ml,digits=3))",:bottom,6)),(0.5,1.0,text("$(round(mw,digits=3))",:bottom,6)),(0.75,1.0,text("$(round(mr,digits=3))",:bottom,6)),(1.0,1.0,text("$(round(smax,digits=3))",:bottom,6))];
    annotate!(a)

    fComp=getElementProperties(p.femType[key][1],p.mesh.meshType,0.5,0.5);

    if typeof(p.degFBoundary[p.femType[key][1]])==degF{1}
        pl=plot();

        for k in 1:nf
            sk=sol[l2g(p.degFBoundary[p.femType[key][1]], k)];
            skm=dot(fComp,sk);
            coord=p.mesh.geometry.coordinates[:,inc[off[k]:off[k+1]-1]];
            fillPoly(coord, ColorGradient(c)[(skm-smin)/(smax-smin)]);
        end

        pf=plot(pl,pc,layout=grid(2,1,heights=[0.98,0.02]))
        display(pf);
    else
        J=Array{Float64,2}(undef,2,2);
        dJ=0.0;
        coord=Array{Float64,2}(undef,2,p.mesh.meshType);
        for j in 1:2
            pl=plot();

            for k in 1:nf
                sk=sol[l2g(p.degFBoundary[p.femType[key][1]], k)];
                dJ=jacobi!(J,p.mesh,k,0.5,0.5,coord);
                fLoc=(1/dJ)*J*fComp
                skm=fLoc*sk;
                coord=p.mesh.geometry.coordinates[:,inc[off[k]:off[k+1]-1]];
                fillPoly(coord, ColorGradient(c)[(skm[j]-smin)/(smax-smin)]);
            end

            pf=plot(pl,pc,layout=grid(2,1,heights=[0.98,0.02]))
            display(pf);
        end
    end
end

#Funktion zum Plotten des Ergebnisses der Akkustikgleichung, dass den Plot nicht in
#der GUI anzeigt, sondern nur unter dem Dateinamen filename (z.B. acoustic.png) im momentanen Directionary speichert
function plotSolution(p::femProblem, key::Symbol, t::Float64, filename::String, c::Symbol=:plasma)
    nf=p.mesh.topology.size[3];
    sol=getfield(p.solution[t],key);
    smin=minimum(sol);
    smax=maximum(sol);
    off=p.mesh.topology.offset["20"];
    inc=p.mesh.topology.incidence["20"];

    pc=plot();
    for i in 0:0.001:1
        color=ColorGradient(c)[i];
        plot!([i,i],[0,1],linecolor=color, legend=false, yaxis=nothing, xaxis=nothing, foreground_color=:white)
    end

    mw=0.5*(smin+smax);
    mr=0.5*(mw+smax);
    ml=0.5*(smin+mw);
    a=[(0.0,1.0,text("$(round(smin,digits=3))",:bottom, 6)),(0.25,1.0,text("$(round(ml,digits=3))",:bottom,6)),(0.5,1.0,text("$(round(mw,digits=3))",:bottom,6)),(0.75,1.0,text("$(round(mr,digits=3))",:bottom,6)),(1.0,1.0,text("$(round(smax,digits=3))",:bottom,6))];
    annotate!(a)

    fComp=getElementProperties(p.femType[key][1],p.mesh.meshType,0.5,0.5);

    if typeof(p.degFBoundary[p.femType[key][1]])==degF{1}
        pl=plot();

        for k in 1:nf
            sk=sol[l2g(p.degFBoundary[p.femType[key][1]], k)];
            skm=dot(fComp,sk);
            coord=p.mesh.geometry.coordinates[:,inc[off[k]:off[k+1]-1]];
            fillPoly(coord, ColorGradient(c)[(skm-smin)/(smax-smin)]);
        end

        pf=plot(pl,pc,layout=grid(2,1,heights=[0.98,0.02]))
        savefig(pf, filename);
    else
        J=Array{Float64,2}(undef,2,2);
        dJ=0.0;
        coord=Array{Float64,2}(undef,2,p.mesh.meshType);
        for j in 1:2
            pl=plot();

            for k in 1:nf
                sk=sol[l2g(p.degFBoundary[p.femType[key][1]], k)];
                dJ=jacobi!(J,p.mesh,k,0.5,0.5,coord);
                fLoc=(1/dJ)*J*fComp
                skm=fLoc*sk;
                coord=p.mesh.geometry.coordinates[:,inc[off[k]:off[k+1]-1]];
                fillPoly(coord, ColorGradient(c)[(skm[j]-smin)/(smax-smin)]);
            end

            pf=plot(pl,pc,layout=grid(2,1,heights=[0.98,0.02]))
            savefig(pf, filename);
        end
    end
end

function fillPoly(c::Array{Float64,2}, fc::RGBA{Float64}; fillrange::Int64=0, fillalpha::Float64=1.0)
    x=push!(c[1,:],c[1,1]);
    y=push!(c[2,:],c[2,1]);
    p=Plots.plot!(x,y,fill=(fillrange,fillalpha,fc),legend=false, linecolor=fc, linewidth=0.0);
    return(p)
end
