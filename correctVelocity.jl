function correctVelocity!(p::femProblem,comp::Float64)
    type=p.femType[:v][1];
    if type==:RT0 || type==:RT0B
        h=0;
    elseif type==:VecP1 || type==:VecDG1
        h=1;
    end
    he=[1,4,2,3];
    le=p.mesh.edgeLength;
    inc=p.mesh.topology.incidence["21"];
    off=p.mesh.topology.offset["21"];
    incd=p.degFBoundary[type].incidence;
    offd=p.degFBoundary[type].offset;
    fac=ones(Float64,length(p.solution[comp].v));
    for k in 1:p.mesh.topology.size[p.mesh.topology.D+1]
        ie=inc[off[k]:off[k+1]-1];
        id=incd[offd[k]:offd[k+1]-1];
        z=1;
        for i in 1:length(ie)
            fac[id[z:z+h]].=le[he[i]];
            z+=h+1;
        end
    end
    p.solution[comp].v=(p.solution[comp].v)./fac;
    return nothing;
end

function correctVelocity!(p::femProblem)
    type=p.femType[:v][1];
    if type==:RT0 || type==:RT0B
        h=0;
    elseif type==:VecP1 || type==:VecDG1
        h=1;
    end
    he=[1,4,2,3];
    le=p.mesh.edgeLength;
    inc=p.mesh.topology.incidence["21"];
    off=p.mesh.topology.offset["21"];
    incd=p.degFBoundary[type].incidence;
    offd=p.degFBoundary[type].offset;
    fac=ones(Float64,length(p.solution[0.0].v));
    for k in 1:p.mesh.topology.size[p.mesh.topology.D+1]
        ie=inc[off[k]:off[k+1]-1];
        id=incd[offd[k]:offd[k+1]-1];
        z=1;
        for i in 1:length(ie)
            fac[id[z:z+h]].=le[he[i]];
            z+=h+1;
        end
    end
    for k in keys(p.solution)
        p.solution[k].v=(p.solution[k].v)./fac;
    end
    return nothing;
end
