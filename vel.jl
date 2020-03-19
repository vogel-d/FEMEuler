function velCa(velSp::Array{Float64,1},lam::Float64,phi::Float64)
  rot=zeros(3,3)
  rot[1,1]=-sin(lam)
  rot[2,1]=-sin(phi)*cos(lam)
  rot[3,1]= cos(phi)*cos(lam)
  rot[1,2]=cos(lam)
  rot[2,2]=-sin(phi)*sin(lam)
  rot[3,2]=cos(phi)*sin(lam)
  rot[1,3]=0.0
  rot[2,3]=cos(phi)
  rot[3,3]=sin(phi)

  velCa=similar(velSp)
  velCa[1]=rot[1,1]*velSp[1]+rot[2,1]*velSp[2]+rot[3,1]*velSp[3]
  velCa[2]=rot[1,2]*velSp[1]+rot[2,2]*velSp[2]+rot[3,2]*velSp[3]
  velCa[3]=rot[1,3]*velSp[1]+rot[2,3]*velSp[2]+rot[3,3]*velSp[3]

  return velCa

end

function velSp(velCa::Array{Float64,1},lam::Float64,phi::Float64)
  rot=zeros(3,3)
  rot[1,1]=-sin(lam)
  rot[2,1]=-sin(phi)*cos(lam)
  rot[3,1]= cos(phi)*cos(lam)
  rot[1,2]=cos(lam)
  rot[2,2]=-sin(phi)*sin(lam)
  rot[3,2]=cos(phi)*sin(lam)
  rot[1,3]=0.0
  rot[2,3]=cos(phi)
  rot[3,3]=sin(phi)

  velSp=similar(velCa)
  velSp[1]=rot[1,1]*velCa[1]+rot[1,2]*velCa[2]+rot[1,3]*velCa[3]
  velSp[2]=rot[2,1]*velCa[1]+rot[2,2]*velCa[2]+rot[2,3]*velCa[3]
  velSp[3]=rot[3,1]*velCa[1]+rot[3,2]*velCa[2]+rot[3,3]*velCa[3]

  return velCa

end
