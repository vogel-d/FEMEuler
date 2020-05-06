function transformRecoveryCoord(n,t1,t2,kubPoint)
    length(kubPoint)==2 && return kubPoint.-n[1:2]
    kubPoint=lu([kubPoint (-1.0).*t1 (-1.0).*t2])\n
    return kubPoint[2:3]
    #=
    if length(kubPoint)==2
        push!(kubPoint,1.0)
        kubPoint=lu([kubPoint (-1.0).*t1 (-1.0).*t2])\n
        return kubPoint[2:3]
    else
        kubPoint=lu([kubPoint (-1.0).*t1 (-1.0).*t2])\n
        return kubPoint[2:3]
    end
    =#
end
