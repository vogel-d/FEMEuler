function getLambdaPhi(m::mesh)
    coord=m.geometry.coordinates;
    lc=size(coord,2);
    scoord=Array{Float64,2}(undef,2,lc);
    for i in 1:lc
        c=coord[:,i]
        theta=acos(c[3]/sqrt(c[1]^2+c[2]^2+c[3]^2))
        phi=atan(c[2],c[1]);
        scoord[:,i]=[phi,theta];
    end
    mG=meshGeometry(scoord,m.geometry.l,m.geometry.r)
    return mesh(m.topology,mG);
end
