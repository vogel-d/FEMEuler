function transformRecoveryCoord!(recPoint,n,t1,t2,kubPoint,A)
    if length(kubPoint)==2
        recPoint[:]=kubPoint.-n[1:2]
        return nothing
    end
    A[:,1]=kubPoint
    A[:,2]=(-1.0).*t1
    A[:,3]=(-1.0).*t2
    push!(recPoint,0.0)
    copyto!(recPoint,n)
    LinearAlgebra.LAPACK.gesv!(A,recPoint)
    popfirst!(recPoint)
    return nothing
end
