#Nur für skalare Variablen
function plotSolutionGif(p::femProblem, key::Symbol, t::Array{T,1} where T=sort(collect(keys(p.solution))); fps::Int=25, filename::String="testgif.gif")
    (smax,smin)=calculateScale(p,key);
    anim = @animate for i in t
                        plotGif(p,key,i,smax,smin);
                    end
        gif(anim,filename,fps=fps);
end

function calculateScale(p::femProblem, key::Symbol)
    max=maximum(getfield(p.solution[0],key));
    min=minimum(getfield(p.solution[0],key));
    t=collect(keys(p.solution));
    for i in 2:length(t)
        hmax=maximum(getfield(p.solution[t[i]],key));
        hmin=minimum(getfield(p.solution[t[i]],key));
        if hmax>max
            max=hmax;
        end
        if hmin<min
            min=hmin;
        end
    end
    return (max, min)
end

function plotGif(p::femProblem, key::Symbol, t::AbstractFloat, smax::AbstractFloat, smin::AbstractFloat, c::Symbol=:plasma)
    nf=p.mesh.topology.size[3];
    sol=getfield(p.solution[t],key);
    smin=minimum(sol);
    smax=maximum(sol);
    off=p.mesh.topology.offset["20"];
    inc=p.mesh.topology.incidence["20"]
    fComp=getElementProperties(p.femType[key][1],p.mesh.meshType,0.5,0.5);
    pl=plot();
    for k in 1:nf
        sk=sol[l2g(p.degFBoundary[p.femType[key][1]], k)];
        skm=dot(fComp,sk);
        coord=p.mesh.geometry.coordinates[:,inc[off[k]:off[k+1]-1]];
        fillPoly(coord, ColorGradient(c)[(skm-smin)/(smax-smin)]);
    end
    pc=plot();
    for i in 0:0.001:1
        color=ColorGradient(c)[i];
        plot!([i,i],[0,1],linecolor=color, legend=false, yaxis=nothing, xaxis=nothing, foreground_color=:white);
    end
    mw=0.5*(smin+smax);
    mr=0.5*(mw+smax);
    ml=0.5*(smin+mw);
    a=[(0.0,1.0,text("$(round(smin,digits=3))",:bottom, 6)),(0.25,1.0,text("$(round(ml,digits=3))",:bottom,6)),(0.5,1.0,text("$(round(mw,digits=3))",:bottom,6)),(0.75,1.0,text("$(round(mr,digits=3))",:bottom,6)),(1.0,1.0,text("$(round(smax,digits=3))",:bottom,6))];
    annotate!(a);
    pf=plot(pl,pc,layout=grid(2,1,heights=[0.98,0.02]));

    return pf;
end
