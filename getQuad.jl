function getQuad(g::Int64)
    n=Int64(ceil((g+1)/2));
    gdots=Array{Float64,1};
    W=Array{Float64,1};

    if g>9
        a1=0.538469310105683;
        a2=0.906179845938664;
        gdots=[-a2,-a1,0.0,a1,a2];
        b1=0.236926885056189;
        b2=0.478628670499366;
        b3=0.568888888888889 ;
        W=[b1, b2, b3, b2, b1];
        @warn "Der eingegebene Genauigkeitsgrad ist nicht erreichbar. Dies ist eine Approximation vom Genauigkeitsgrad 9."
    elseif g>7
        a1=0.538469310105683;
        a2=0.906179845938664;
        gdots=[-a2,-a1,0.0,a1,a2];
        b1=0.236926885056189;
        b2=0.478628670499366;
        b3=0.568888888888889 ;
        W=[b1, b2, b3, b2, b1];
    elseif g>5
        a1=0.8611363115940526;
        a2=0.3399810435848563;
        gdots=[-a1,-a2,a2,a1];
        b1=0.347854845137454;
        b2=0.652145154862546;
        W=[b1, b2, b2, b1];
    elseif g>3
        a=0.7745966692414834;
        gdots=[-a,0.0,a];
        b1=5/9;
        b2=8/9;
        W=[b1, b2, b1];
    elseif g>1
        a=0.5773502691896258;
        gdots=[-a,a];
        W=[1.0, 1.0];
    else
        error("Geben Sie ein g>0 an.")
    end
    x=0.0;
    y=1.0
    hx=(y-x)/2;
    sx=(x+y)/2;

    for k in 1:n
        gdots[k]=hx*gdots[k]+sx;
    end
    
    return gdots, hx*W;
end
