function transformRecoveryCoord(n,t1,t2,kubPoint)
    #=
    length(kubPoint)==2 && return kubPoint.-n
    kubPoint=lu([kubPoint (-1.0).*t1 (-1.0).*t2])\n
    return kubPoint[2:3]
    =#
    if length(kubPoint)==2
        push!(kubPoint,1.0)
        kubPoint=qr([kubPoint (-1.0).*t1 (-1.0).*t2])\n
        return kubPoint[1:2]
    else
        #kubPoint=qr([kubPoint (-1.0).*t1 (-1.0).*t2])\n
        kubPoint=lu([kubPoint (-1.0).*t1 (-1.0).*t2])\n
        return kubPoint[2:3]
    end
end
