import Base.intersect
function intersect(pcoord,p,J,x)
    x[1]=0.5;
    x[2]=0.5;
    x[3]=1.0;
    for i=1:5
      f=fun(x[1],x[2],x[3],pcoord[:,1],pcoord[:,2],pcoord[:,3],pcoord[:,4],p);
      jacobi!(J,x[1],x[2],x[3],pcoord[:,1],pcoord[:,2],pcoord[:,3],pcoord[:,4],p);
      x-=J\f
    end
    #Newton method
    return x[1],x[2]
end

function fun(ksi,eta,t,p1,p2,p3,p4,p)
    return (t*p-((1-ksi)*(1-eta)*p1+ksi*(1-eta)*p2+ksi*eta*p3+(1-ksi)*eta*p4))
end

function jacobi!(J,ksi,eta,t,p1,p2,p3,p4,p)
    J[:,1]=-((-1)*(1-eta)*p1+(1-eta)*p2+eta*p3+(-1)*eta*p4);
    J[:,2]=-((1-ksi)*(-1)*p1+ksi*(-1)*p2+ksi*p3+(1-ksi)*p4);
    J[:,3]=p;
end
