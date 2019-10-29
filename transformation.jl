function transformation(m::mesh, coord::Array{Float64,2}, x::Float64, y::Float64)
    #es wird von einem Gitter mit einheitlichen Elementen ausgegangen
    #davon ausgehend lässt sich n allgemein durch die ersten beiden Offset-Einträge feststellen
    n=m.meshType;
    if n==3
        r1=coord[:,1];
        return r1+(coord[:,2]-r1)*x+(coord[:,3]-r1)*y
    elseif n==4
        r1=coord[:,1];
        return r1+(coord[:,2]-r1)*x+(coord[:,4]-r1)*y+(coord[:,3]-coord[:,4]-coord[:,2]+r1)*x*y
    else
        error("Diese Funktion benötigt ein Mesh mit drei- oder viereckigen Elementen");
    end
end

function transformation(m::mesh, coord::SubArray{Float64,2,Array{Float64,2},Tuple{Base.Slice{Base.OneTo{Int64}},SubArray{Int64,1,Array{Int64,1},Tuple{UnitRange{Int64}},true}},false}, x::Float64, y::Float64)
    #es wird von einem Gitter mit einheitlichen Elementen ausgegangen
    #davon ausgehend lässt sich n allgemein durch die ersten beiden Offset-Einträge feststellen
    n=m.meshType;
    if n==3
        r1=coord[:,1];
        return r1+(coord[:,2]-r1)*x+(coord[:,3]-r1)*y
    elseif n==4
        r1=coord[:,1];
        return r1+(coord[:,2]-r1)*x+(coord[:,4]-r1)*y+(coord[:,3]-coord[:,4]-coord[:,2]+r1)*x*y
    else
        error("Diese Funktion benötigt ein Mesh mit drei- oder viereckigen Elementen");
    end
end
#=
function transformation(m::mesh, coord::SubArray{Float64,2,Array{Float64,2},Tuple{Base.Slice{Base.OneTo{Int64}},Array{Int64,1}},false}, x::Float64, y::Float64)
    #es wird von einem Gitter mit einheitlichen Elementen ausgegangen
    #davon ausgehend lässt sich n allgemein durch die ersten beiden Offset-Einträge feststellen
    n=m.meshType;
    if n==3
        r1=coord[:,1];
        return r1+(coord[:,2]-r1)*x+(coord[:,3]-r1)*y
    elseif n==4
        r1=coord[:,1];
        return r1+(coord[:,2]-r1)*x+(coord[:,4]-r1)*y+(coord[:,3]-coord[:,4]-coord[:,2]+r1)*x*y
    else
        error("Diese Funktion benötigt ein Mesh mit drei- oder viereckigen Elementen");
    end
end
=#
