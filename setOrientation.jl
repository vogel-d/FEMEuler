function setOrientation!(m::mesh)
    inc=m.topology.incidence["20"]
    off=m.topology.offset["20"]
    for k in 1:m.topology.size[m.topology.dim+1]
        coord=m.geometry.coordinates[:,inc[off[k]:off[k+1]-1]]
        e0=coord[:,2]-coord[:,1];
        e1=coord[:,4]-coord[:,1];
        n1=cross(e0,e1)
        n=transformation(m,coord,0.5,0.5)
        push!(m.orientation,sign(dot(n1,n)))
    end
    return nothing
end
function getOrientation(m::mesh)
    inc=m.topology.incidence["20"]
    off=m.topology.offset["20"]
    r=[];
    for k in 1:m.topology.size[m.topology.dim+1]
        coord=m.geometry.coordinates[:,inc[off[k]:off[k+1]-1]]
        e0=coord[:,2]-coord[:,1];
        e1=coord[:,4]-coord[:,1];
        n1=cross(e0,e1)
        #n=transformation(m,coord,0.5,0.5)
        n=coord[:,1]
        push!(r,sign(dot(n1,n)))
    end
    return r
end
