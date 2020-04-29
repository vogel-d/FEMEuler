function transformRecoveryCoord(n,kubPoint)

    ind=partialsortperm(abs.(n),1:2,rev=true)
    t1=zeros(length(n));
    t1[ind[1]]=n[ind[2]]
    t1[ind[2]]=-n[ind[1]]

    t2=cross(n,t1);

    kubPoint=lu([kubPoint (-1.0).*t1 (-1.0).*t2])\n
    return kubPoint[2:3]
end
