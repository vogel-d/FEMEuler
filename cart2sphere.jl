function cart2sphere(x,y,z)
  r=sqrt(x^2+y^2+z^2);
  lam=atan(y,x)
  if y<0.0
    lam+=8.0*atan(1.0)
  end
  phi=asin(z/r)
  return lam,phi,r
end

function distCircle(lam,phi,lam0,phi0,R)
  return R*acos(sin(phi)*sin(phi0)+cos(phi)*cos(phi0)*cos(lam-lam0))
end
