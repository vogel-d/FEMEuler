function set_cuco(m::Int) #pcoco.f90
   # Cube, {4,3}
   r2=1.4142135623730950488016887242096980785696718753769480731766; #sqrt(2)
   r3=1.7320508075688772935274463415058723669428052538103806280558;  #sqrt(3)
   or2=1.0/r2;
   c1=complex(1.0,0.0);
   ci=complex(0.0,1.0);
   p270=complex(-0.5,-0.5*r3);
   p045=complex(0.5,0.5*r3);
   p315=complex(0.5,-r3*0.5);
   z045=complex(or2,or2);
   z225=complex(-or2,-or2);

   #Int-Parameter
   n=120; # <- number of test points in each test arc
   n2=n*2; # <- number of points that complete a minimal cycle
   nm=n-1;n2m=n2-1;nit=40;
   f=4;v=3;v0=Int((2*f)/(f-2)) # <- Face and vertex Symmetries
   #Float64-Parameter
   of=1.0/f;ov=1.0/v;ov0=1.0/v0;
   dalpha=pi*ov0/n;dbeta=pi*of/n;
   r0=.96;r0v0=r0^v0;r0f=r0^(f);
   r3p=r3+1.0;r2or3p=r2/r3p;mo2=-0.5
   #Komplexe Parameter
   ca=complex(0.5,0.5);cb=complex(0.5,mo2); wa=complex(r2or3p,0.0); w1=complex(or2,0.0);

   tco1=Float64[0.5222445411749E+00,-0.1349990960563E+00,-0.1970373552294E-01,-0.3167426867673E-02,-0.2797723788949E-02,-0.1720480734504E-02,
         -0.1164080730235E-02,-0.8325527861928E-03,-0.6217962129431E-03]; #(cu co1 seed)
   tco2=Float64[0.1094170834081E+00,-0.3352187479631E-01,-0.9955635663016E-03, -0.9915868353376E-03,-0.5372057577482E-03,-0.3367620231221E-03,
         -0.2279393461075E-03,-0.1632484731336E-03,-0.1219559186293E-03]; #(cu co2 seed)

   i=min(9,m);
   co1=zeros(m);
   co2=zeros(m);
   co1[1:i]=tco1[1:i];
   co2[1:i]=tco2[1:i];
   zarg1=OffsetArray{Complex{Float64},1}(undef,0:n);
   zarg2=OffsetArray{Complex{Float64},1}(undef,0:n);
   sect1=OffsetArray{Int,1}(undef,0:n);
   sect2=OffsetArray{Int,1}(undef,0:n);
   for i in 0:n
      z0=r0*exp(ci*i*dalpha); z1=c1-z0; za=z045*r2*(z0-ca)
      if (abs(z1)<=abs(za))
         zarg1[i]=z1^v0
         if (imag(z1)>=-real(z1))
            sect1[i]=1
         else
            sect1[i]=2
         end
      else
         zarg1[i]=za^(f)
         if (imag(z0)<0.5)
            sect1[i]=3
         else
            sect1[i]=4
         end
      end
      if (i==0 || i==n)
         zarg1[i]=complex(real(zarg1[i]),0.0)
      end
   end
   for i in 0:n
      z0=ca+z225*or2*r0*exp(ci*i*dbeta); zb=z225*r2*(z0-cb)
      if (abs(z0)<=abs(zb))
         zarg2[i]=z0^v0;
         sect2[i]=1
      else
         zarg2[i]=zb^(f);
         sect2[i]=2
      end
      if(i==0 || i==n)
         zarg2[i]=complex(real(zarg2[i]),0.0)
      end
   end

   ww=OffsetArray{Complex{Float64},1}(undef,0:n2m);
   for it in 1:nit
   # Recompute the w of each test point of each 45-degree arc:
      for i in 0:n
         if sect1[i]==1
            w=tay(co1,zarg1[i])
            w=w^ov
            w=(w1-w)/(w1*w+c1)
         elseif sect1[i]==2
            w=tay(co1,zarg1[i])
            w=w^ov
            w=p270*w; w=(w1-w)/(w1*w+c1)
         elseif sect1[i]==3
            w=tay(co2,zarg1[i])
            w=w^of
            w=ci*w; w=(wa-w)/(wa*w+c1); w=p045*w
         elseif sect1[i]==4
            w=tay(co2,zarg1[i])
            w=w^of
            w=-w;   w=(wa-w)/(wa*w+c1); w=p045*w
         end
         ww[i]=w^v
         if (i==0 || i==n)
            ww[i]=complex(real(ww[i]),0.0)
         else
            ww[n2-i]=conj(ww[i])
         end
      end
      co1=tay_reset(m,ww,r0v0)
      for i in 0:n
         if sect2[i]==1
            w=tay(co1,zarg2[i])
            w=w^ov
            w=p315*w; w=(wa-w)/(wa*w+c1)
         elseif sect2[i]==2
            w=tay(co2,zarg2[i]);
            if (i==n)
               w=complex(real(w),0.0)
            end
            w=w^of;
            if (i==n)
               w=conj(w)
            end
            w=w*z045; w=(c1-w)/(w+c1); w=w*z045
         end
         ww[i]=w^(f)
         if (i==0 || i==n)
            ww[i]=complex(real(ww[i]),0.0)
         else
            ww[n2-i]=conj(ww[i])
         end
      end
      co2=tay_reset(m,ww,r0f)
   end
   return co1, co2;
end
