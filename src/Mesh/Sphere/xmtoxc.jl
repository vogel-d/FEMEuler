function xmtoxc(xm::Array{Float64,1},ipan::Int) #cuco.f90
    # From panel index, ipan, and map coordinates, xm, return earth-centered
    # cartesian vector, xc.

    #lnotinicuco && inicuco()
    mcuco=50
    rotp0,co1,co2=inicuco(mcuco);
    z=0.5*complex(xm[1]+1.0,xm[2]+1.0)
    xc=cuztoc(co1,co2,z)
    xc=rotp0[:,:,ipan]*xc;
    return xc;
end

function inicuco(mcuco::Int) #cuco.f90
    # Initialize conformal cubic expansion coefficients and the rotations, rotp0
    # for each of the six panel coordinate frames.

    co1,co2=set_cuco(mcuco)

    rotp0=zeros(3,3,6)
    rotp0[:,:,1]=[0 0 1; 1 0 0; 0 1 0];
    rotp0[:,:,2]=[-1 0 0; 0 0 1; 0 1 0];
    rotp0[:,:,3]=[0 0 -1; -1 0 0; 0 1 0];
    rotp0[:,:,4]=[1 0 0; 0 0 -1; 0 1 0];
    rotp0[:,:,5]=[0 -1 0; 1 0 0; 0 0 1];
    rotp0[:,:,6]=[0 1 0; 1 0 0; 0 0 -1];
    #lnotinicuco=false;
    return rotp0,co1,co2
end

function cuztoc(co1::Array{Float64,1},co2::Array{Float64,1},z::Complex{Float64}) #pcoco.f90
   # Conformal cubic mapping to the unit sphere.
   # Map from the unit square map face to cartesians oriented
   # with respect to the face center.
   #m=length(co1)

   r2=1.4142135623730950488016887242096980785696718753769480731766; #sqrt(2)
   r3=1.7320508075688772935274463415058723669428052538103806280558;  #sqrt(3)
   or2=1.0/r2;
   c1=complex(1.0,0.0);
   ci=complex(0.0,1.0);
   p045=complex(0.5,0.5*r3);
   p315=complex(0.5,-r3*0.5);
   z045=complex(or2,or2);
   z225=complex(-or2,-or2);
   z135=complex(-or2,or2);
   z315=complex(or2,-or2);

   #Int-Parameter
   f=4;v=3; # <- Face and vertex Symmetries
   #Float64-Parameter
   v0=(2*f)/(f-2);
   of=1.0/f;ov=1.0/v;ov0=1.0/v0;eps=1.e-10;war=r2/(1.0+r3);
   #Komplexe Parameter
   ca=complex(0.5,0.5);
   wa=complex(war,0.0);

   za=z-ca;
   k=0;
   if (real(za)>0.0)
      k=1;
   end
   if (imag(za)>0.0)
      k+=2;
   end
   if k==1
      za=-ci*za;
   elseif k==2
      za=ci*za;
   elseif k==3
      za=-za;
   end
   z0=ca+za
   za=z135*r2*za
   ka=(abs(za)<abs(z0));
   if (ka)         # Expand about the face center:
      kn=(imag(za)<0.0); # negative octant
      zp=za^(f);
      ww=tay(co2,zp);
      if (abs(ww)>0.0)
         w=ww^(of);
         kn1=(imag(w)<0.0);
         if (kn != kn1)
            w =conj(w);
         end
      else
         w=complex(0.0,0.0);
      end
   else                    # Expand about the vertex
      kn=(imag(z0)>real(z0))# negative sector of rotated z0 frame
      if (kn)
         z0 =-ci*z0     # <-perform the rotation
      end
      zp=z0^(v0);
      ww=tay(co1,zp);
      if (abs(ww)>0.0)
         w=ww^(ov);
         kn1=(imag(w)<0.0);
         if (kn != kn1)
            w =conj(w);
         end
      else
         w=complex(0.0,0.0);
      end
   # Convert from vertex-centered frame to face-centered frame:
      if (kn)
         w=p045*w;
      else
         w=p315*w;
      end
      w=(wa-w)/(wa*w+c1);
   end
   if k==0
      w=z225*w;
   elseif k==1
      w=z315*w;
   elseif k==2
      w=z135*w;
   elseif k==3
      w=z045*w;
   end
   # Then convert to cartesians
   xc=ztoc(w,false)
   return xc;
end

function ztoc(z::Complex{Float64},infz::Bool) #pmat4.f90
   infz && return [0.0, 0.0, -1.0];
   r=real(z); q=imag(z); rs=r*r+q*q
   rsc=1.0-rs
   rsbi=1.0/(1.0+rs)
   v=Array{Float64,1}(undef,3);
   v[1]=2*rsbi*r
   v[2]=2*rsbi*q
   v[3]=rsc*rsbi
   return v;
end
