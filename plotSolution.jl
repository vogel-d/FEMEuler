function plotSolution(p::femProblem, key::Symbol, t::Float64, c::Symbol=:plasma)
    nf=p.mesh.topology.size[3];
    sol=getfield(p.solution[t],key);
    smin=minimum(sol);
    smax=maximum(sol);
    trans(x)=(x-smin)/(smax-smin);
    i=p.degFBoundary[p.femType[key][1]].incidence;
    o=p.degFBoundary[p.femType[key][1]].offset;
    off=p.mesh.topology.offset["20"];
    inc=p.mesh.topology.incidence["20"];
    comp=p.degFBoundary[p.femType[key][1]].components;

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

    if comp[1]==0
        pl=plot();

        for k in 1:nf
            rstart1=o[k];
            rend1=o[k+1]-1;
            h=collect(rstart1:rend1);
            sk=sol[i[h]];
            skm=sum(sk)/length(sk);
            rstart2=off[k];
            rend2=off[k+1]-1;
            coord=p.mesh.geometry.coordinates[:,inc[rstart2:rend2]];
            fillPoly(coord, ColorGradient(c)[trans(skm)]);
        end

        pf=plot(pl,pc,layout=grid(2,1,heights=[0.98,0.02]))
        display(pf);
    else
        for j in 1:2
            sc=findall(comp.==j);
            pl=plot();

            for k in 1:nf
                rstart1=o[k];
                rend1=o[k+1]-1;
                h=collect(rstart1:rend1)[sc];
                sk=sol[i[h]];
                skm=sum(sk)/length(sk);
                rstart2=off[k];
                rend2=off[k+1]-1;
                coord=p.mesh.geometry.coordinates[:,inc[rstart2:rend2]];
                fillPoly(coord, ColorGradient(c)[trans(skm)]);
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
    trans(x)=(x-smin)/(smax-smin);
    i=p.degFBoundary[p.femType[key][1]].incidence;
    o=p.degFBoundary[p.femType[key][1]].offset;
    off=p.mesh.topology.offset["20"];
    inc=p.mesh.topology.incidence["20"]
    comp=p.degFBoundary[p.femType[key][1]].components;

    pc=plot()
    for i in 0:0.001:1
        color=ColorGradient(c)[i];
        plot!([i,i],[0,1],linecolor=color, legend=false, yaxis=nothing, xaxis=nothing, foreground_color=:white)
    end

    mw=0.5*(smin+smax);
    mr=0.5*(mw+smax);
    ml=0.5*(smin+mw);
    a=[(0.0,1.0,text("$(round(smin,digits=3))",:bottom, 6)),(0.25,1.0,text("$(round(ml,digits=3))",:bottom,6)),(0.5,1.0,text("$(round(mw,digits=3))",:bottom,6)),(0.75,1.0,text("$(round(mr,digits=3))",:bottom,6)),(1.0,1.0,text("$(round(smax,digits=3))",:bottom,6))];
    annotate!(a)

    if comp[1]==0
        pl=plot();

        for k in 1:nf
            rstart1=o[k];
            rend1=o[k+1]-1;
            h=collect(rstart1:rend1);
            sk=sol[i[h]];
            skm=sum(sk)/length(sk);
            rstart2=off[k];
            rend2=off[k+1]-1;
            coord=p.mesh.geometry.coordinates[:,inc[rstart2:rend2]];
            fillPoly(coord, ColorGradient(c)[trans(skm)]);
        end

        pf=plot(pl,pc,layout=grid(2,1,heights=[0.98,0.02]))
        savefig(pf, filename);
    else
        for j in 1:2
            sc=findall(comp.==j);
            pl=plot();

            for k in 1:nf
                rstart1=o[k];
                rend1=o[k+1]-1;
                h=collect(rstart1:rend1)[sc];
                sk=sol[i[h]];
                skm=sum(sk)/length(sk);
                rstart2=off[k];
                rend2=off[k+1]-1;
                coord=p.mesh.geometry.coordinates[:,inc[rstart2:rend2]];
                fillPoly(coord, ColorGradient(c)[trans(skm)]);
            end
            
            pf=plot(pl,pc,layout=grid(2,1,heights=[0.98,0.02]))
            savefig(pf, filename);
        end
    end
end

function fillPoly(c::Array{Float64,2}, fc::RGBA{Float64}; fillrange::Int64=0, fillalpha::Float64=1.0)
    x=push!(c[1,:],c[1,1]);
    y=push!(c[2,:],c[2,1]);
    p=Plots.plot!(x,y,fill=(fillrange,fillalpha,fc),legend=false, linecolor=:white, linewidth=0);
    return(p)
end
