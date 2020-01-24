function transformation(m::mesh, coord::Array{AbstractFloat,2}, x::AbstractFloat, y::AbstractFloat)
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

function transformation(m::mesh, coord::SubArray{AbstractFloat,2,Array{AbstractFloat,2},Tuple{Base.Slice{Base.OneTo{Int}},SubArray{Int,1,Array{Int,1},Tuple{UnitRange{Int}},true}},false}, x::AbstractFloat, y::AbstractFloat)
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
function transformation(m::mesh, coord::SubArray{AbstractFloat,2,Array{AbstractFloat,2},Tuple{Base.Slice{Base.OneTo{Int}},Array{Int,1}},false}, x::AbstractFloat, y::AbstractFloat)
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
