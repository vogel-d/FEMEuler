function xmtoxc(xm::Array{Float64,1},ipan::Int)
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

function cuztoc(co1::Array{Float64,1},co2::Array{Float64,1},z::Complex{Float64})
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
   of=1.0/f;ov=1.0/v;ov0=1.0/v0;eps=1.e-10;war=r2/(1.0+r3)
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
      za= ci*za;
   elseif k==3
      za=-za;
   end
   z0=ca+za
   za =z135*r2*za
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

function ztoc(z::Complex{Float64},infz::Bool)
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

function inicuco(mcuco::Int)
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

function set_cuco(m::Int)
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
   nm=n-1;nit=40;
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
   zarg1=Array{Complex{Float64},1}(undef,n+1);
   zarg2=Array{Complex{Float64},1}(undef,n+1);
   sect1=Array{Int,1}(undef,n+1);
   sect2=Array{Int,1}(undef,n+1);
   for i in 1:n+1
      z0=r0*exp(ci*(i-1)*dalpha); z1=c1-z0; za=z045*r2*(z0-ca)
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
      if (i==1 || i==n+1)
         zarg1[i]=complex(real(zarg1[i]),0.0)
      end
   end
   for i in 1:n+1
      z0=ca+z225*or2*r0*exp(ci*(i-1)*dbeta); zb=z225*r2*(z0-cb)
      if (abs(z0)<=abs(zb))
         zarg2[i]=z0^v0;
         sect2[i]=1
      else
         zarg2[i]=zb^(f);
         sect2[i]=2
      end
      if(i==1 || i==n+1)
         zarg2[i]=complex(real(zarg2[i]),0.0)
      end
   end

   ww=Array{Complex{Float64},1}(undef,n2);
   for it in 1:nit
   # Recompute the w of each test point of each 45-degree arc:
      for i in 1:n+1
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
         if (i==1 || i==n+1)
            ww[i]=complex(real(ww[i]),0.0)
         else
            ww[n2-i]=conj(ww[i])
         end
      end
      co1=tay_reset(m,ww,r0v0)
      for i in 1:n+1
         if sect2[i]==1
            w=tay(co1,zarg2[i])
            w=w^ov
            w=p315*w; w=(wa-w)/(wa*w+c1)
         elseif sect2[i]==2
            w=tay(co2,zarg2[i]);
            if (i==n+1)
               w=complex(real(w),0.0)
            end
            w=w^of;
            if (i==n+1)
               w=conj(w)
            end
            w=w*z045; w=(c1-w)/(w+c1); w=w*z045
         end
         ww[i]=w^(f)
         if (i==1 || i==n+1)
            ww[i]=complex(real(ww[i]),0.0)
         else
            ww[n2-i]=conj(ww[i])
         end
      end
      co2=tay_reset(m,ww,r0f)
   end
   return co1, co2;
end

function tay(aco::Array{Float64,1},z::Complex{Float64})
   w=0.0;
   zp=complex(1.0,0.0);
   for k in 1:length(aco)
      zp=zp*z
      w+=aco[k]*zp
   end
   return w;
end

function tay_reset(nco::Int,zzf::Array{Complex{Float64},1},r::Float64)
   # Reset the Taylor series coefficients, assumed real, from the results
   # of a complex fourier transformation of cycle nf. The nominal radius of
   # the circuit is r.
   zzf=dfft(zzf)
   ri=1.0/r;
   rp=1;
   co=Array{Float64,1}(undef,nco);
   for i in 1:nco
      rp=rp*ri;
      co[i]=real(zzf[i])*rp;
   end
   return co;
end


function dfft(z::Array{Complex{Float64},1})
   # Like dfft, but with r and q combined into complex z.

   # complex(dpc),dimension(0:n-1),intent(INOUT):: z
   r=real(z);
   q=imag(z);

   n=length(z);
   rfac=1/n;
   #  for fourier synthesis, scale, and reverse the order of wavenumbers:
   r[1]=r[1]*rfac; r[2:n]=r[n:-1:2]*rfac
   q[1]=q[1]*rfac; q[2:n]=q[n:-1:2]*rfac

   #real(sp),dimension(0:n-1),intent(INOUT):: r,q

   #integer                   :: ln2,ln3,ln4,ln5,ln7,lnb,lnd
   #integer, dimension(0:n-1) :: jumble
   #real(sp),dimension(0:n-1) :: w

   jumble,w,ln2,ln3,ln4,ln5,ln7,lnb,lnd=fftco(n);
   r,q=gfft(n,jumble,w,ln2,ln3,ln4,ln5,ln7,lnb,lnd,r,q);

   for i in 1:length(r)
      z[i]=complex(r[i],q[i])
   end
   return z;
end


function fftco(n::Int)
   #   Provide the FFT integer coefficients, j, and real coefficients, w,
   # corresponding to a period, n, which factors into powers of 2, 3, 5, 7, 11, 13
   # only. For short period transforms, retain (in js, ws) up to two previously
   # used sets of these coefficients to avoid the cost of recalculating them
   # in repeated applications. Then the first element of the relevant
   # column of js is always the period, n, making it easy to recognize
   # whether or not the new requirements for j and w are met in previously
   # recorded pairs, js and ws.
   #   For long period transforms, it is worth supplying an integer array, js,
   # and a real array, ws, both of size n, to store the transform parameters
   # between applications so that they do not have to be recomputed each time.
   # This is found to save significant time in large scale repeated applications.
   # In that case, the relevant version of this routine is SFFTCOP.
   # --> n:      period of data
   # --> j:      array of n permutation indices to unscramble the fft output.
   # --> w:      array of n real trigonometric coefficients for the fft.

   bsize=2048;
   ult=1; # Column index of js and ws for latest fft coefficients used.
   #data js/bsize2*0/,ws/bsize2*zero/ bsize2=bsize*2;
   js=zeros(Int,bsize,2);
   ws=zeros(bsize,2);
   ln2,ln3,ln4,ln5,ln7,lnb,lnd,ln,nt,nm,nh=sget2357bd(n);
   j=Array{Int,1}(undef,n);
   w=Array{Float64,1}(undef,n);
   if (n <= bsize)
      if (n == js[1,ult])
         j[1:n]=js[1:n,ult]; w[1:n]=ws[1:n,ult] # copy existing values
         return
      end
      ult=3-ult;                            # Try the alternative location
      if (n == js[1,ult])
         j[1:n]=js[1:n,ult]; w[1:n]=ws[1:n,ult] # copy existing values
         return
      end
      j,w=rumble(ln2,ln3,ln4,ln5,ln7,lnb,lnd,ln,nt,nm,nh)
      js[1:n,ult]=j[1:n]; ws[1:n,ult]=w[1:n]
   else
      j,w=rumble(ln2,ln3,ln4,ln5,ln7,lnb,lnd,ln,nt,nm,nh) # Too big to record in js and ws
   end
   return j,w,ln2,ln3,ln4,ln5,ln7,lnb,lnd;
end

function sget2357bd(nin::Int)
   #   Factorize NIN in terms of 2's, 3's, 4's, 5's, 7's, 11's and 13's  and
   #   verify that these include the only prime factors. Set LN2, LN3, LN4, LN5,
   #   LN7, LNB, LND respectively to the number of powers of: 2, 3, 4, 5, 7,
   #   11, 13; set LN to be the sum,
   #   LN2+LN3+LN4*2+LN5+LN7+LNB+LND, and initialize n, nm, nh,
   #   all in fft23457bd.
   # --> NIN      Number of data along the line of Fourier transformation
   k=nin;
   k,lnd=snfftln(k,13);
   k,lnb=snfftln(k,11);
   k,ln7=snfftln(k,7);
   k,ln5=snfftln(k,5);
   k,ln4=snfftln(k,4);
   k,ln3=snfftln(k,3);
   k,ln2=snfftln(k,2);

   ln=ln2+ln4*2+ln3+ln5+ln7+lnb+lnd;
   n=(2^ln2)*(3^ln3)*(4^ln4)*(5^ln5)*(7^ln7)*(11^lnb)*(13^lnd);
   nm=n-1; nh=Int(n/2);
   (n != nin) && error("prime factors of fft period are not only 2, 3, 5, 7, 11, 13");
   (k!=1) && error("prime factors are not only 2, 3, 5, 7, 11, 13");
   return ln2,ln3,ln4,ln5,ln7,lnb,lnd,ln,n,nm,nh;
end

function snfftln(n1::Int,m::Int)
   # Divide out as many factors, m, of given n1 as is possible, returning
   # this number of factors as nfftln, and returning n1 as the quotient.
   nfftln=0;
   for i in 1:30
      n2=Int(floor(n1/m));
      (n2*m != n1) && return n1,nfftln;
      n1=n2
      nfftln=i;
   end
   return n1,nfftln;
end

function rumble(ln2::Int,ln3::Int,ln4::Int,ln5::Int,ln7::Int,lnb::Int,lnd::Int,ln::Int,n::Int,nm::Int,nh::Int)
   #   Initialize coefficient arrays JUMBLE and TUMBLE for use with the
   #   fast-Fourier-transform routines when the number N of data has the prime-
   #   factorization 2**LN2*3**LN3*5**LN5*7**LN7*11**LNB*13**LND
   # <-- JUMBLE:  Permutation of N data encoded as a sequence of transpositions
   # <-- TUMBLE:  Trigonometric coefficients for use in FFT. The first half are
   #              the cosines, the second half the sines, of uniformly increasing
   #              relevant angles.

   ml=30;
   pi2=6.2831853071795864769252867665590057683943387987502116419498;
   pi2on=pi2/n;
   jumble=Array{Int,1}(undef,n);
   tumble=Array{Float64,1}(undef,n);
   for i in 1:nh
      ang=pi2on*i;
      tumble[i]=cos(ang); tumble[i+nh]=sin(ang);
   end
   id=1;  is=0;
   lnvec=[lnd,lnb,ln7,ln5,ln3,ln2+ln4*2];
   lnnumb=[13,11,7,5,3,2];
   nd=Array{Int,1}(undef,ml);
   md=Array{Int,1}(undef,ml)
   for k in 1:length(lnvec)
      for i in 1:lnvec[k]
         is+=1; md[is]=id; id=id*lnnumb[k];
      end
   end
   id=1
   for k in length(lnvec):-1:1
      for i in 1:lnvec[k]
         nd[is]=id; id=id*lnnumb[k]; is-=1;
      end
   end
   jumble[1]=n;
   for i in 1:nm
      ir=i; j=0
      for l in 1:ln
         kd=Int(floor(ir/nd[l])); ir-=kd*nd[l]; j+=kd*md[l];
      end
      jumble[i+1]=j
   end
   for i in 1:nm
      j=jumble[i+1]
      if (j < i)
         j=jumble[j+1];
         while (j < i)
            j=jumble[j+1];
         end
         jumble[i+1]=j
      end
   end
   return jumble,tumble;
end

function gfft(n::Int,jumble::Array{Int,1},w::Array{Float64,1},ln2::Int,ln3::Int,ln4::Int,ln5::Int,ln7::Int,lnb::Int,lnd::Int,r::Array{Float64,1},q::Array{Float64,1})
s30=0.5;
ms30=-0.5;
s60=1.7320508075688772935274463415058723669428052538103806280558*0.5;
nh=Int(n/2);
# PERMUTE THE DATA:
for i=2:n
   j=jumble[i]
   if  (j > i)
      t1=r[i]; r[i]=r[j]; r[j]=t1
      t1=q[i]; q[i]=q[j]; q[j]=t1
   end
end

#  TRANSFORM THE DATA:
ma=1; mb=n
#  RADIX 4
for l=1:ln4
   mb=Int(floor(mb/4)); ma4=ma*4; mb2=mb*2
   for j=1:ma
      jmb=j*mb; jmb2=j*mb2
      rf1=w[jmb];          qf1=w[nh+jmb]        # f**1
      rf2=w[jmb2];         qf2=w[nh+jmb2]       # f**2
      rf3=rf1*rf2-qf1*qf2; qf3=rf1*qf2+qf1*rf2  # f**3
      for i=1:ma4:n
         k0=(i-1)+j;    k1=k0+ma;   k2=k1+ma;   k3=k2+ma
         r0=r[k0]; r1=r[k1];  r2=r[k2];  r3=r[k3]
         q0=q[k0]; q1=q[k1];  q2=q[k2];  q3=q[k3]
         t1    =r3*rf3-q3*qf3 # q13
         q3    =r3*qf3+q3*rf3 # r13
         r3    =r2*rf1-q2*qf1 # r12
         r2    =r2*qf1+q2*rf1 # q12
         q2    =q3    -r2     # r23
         r2    =q3    +r2     # q22
         q3    =r3    +t1     # r22
         t1    =r3    -t1     # q23
         r3    =r1*rf2-q1*qf2 # r11
         q1    =r1*qf2+q1*rf2 # q11
         r1    =r0    -r3     # r21
         r0    =r0    +r3     # r20
         r[k3] =r1    -q2     # r43
         r[k1] =r1    +q2     # r41
         q2    =q0    +q1     # q20
         q1    =q0    -q1     # q21
         q[k0] =q2    +r2     # q40
         q[k2] =q2    -r2     # q42
         r[k2] =r0    -q3     # r42
         r[k0] =r0    +q3     # r40
         q[k3] =q1    -t1     # q43
         q[k1] =q1    +t1     # q41
      end
   end
   ma=ma4
end
if (ln2==1)
#  RADIX 2
   mb=Int(floor(mb/2)); ma2=ma*2
   for j=1:ma
      jmb=j*mb
      for i=1:ma2:n
         k0=j+(i-1);        k1=k0+ma
         rf1=w[jmb];    qf1=w[nh+jmb] # f**1
         r0=r[k0];     q0=q[k0]
         r1=r[k1];     q1=q[k1]
         t1   =r1*qf1+q1*rf1 # q11
         q1   =r1*rf1-q1*qf1 # r11
         r[k1]=r0    -q1     # r41
         r[k0]=r0    +q1     # r40
         q[k1]=q0    -t1     # q41
         q[k0]=q0    +t1     # q40
      end
   end
   ma=ma2
end
#  RADIX 3
rep=ms30;  rec=3*s30;  qep=s60 # <- epsilon
for l=1:ln3
   mb=Int(floor(mb/3)); ma3=ma*3
   for j=1:ma
      jmb=j*mb
      rf1=w[jmb];          qf1=w[nh+jmb] # f**1
      rf2=rf1*rf1-qf1*qf1; qf2=2*rf1*qf1 # f**2
      for i=1:ma3:n
         k0=(i-1)+j;     k1=k0+ma;    k2=k1+ma
         r1=r[k1];  q1=q[k1]
         r2=r[k2];  q2=q[k2]
         t1    = r2*qf2+q2 *rf2 # r12
         q2    = r2*rf2-q2 *qf2 # q12
         r2    = r1*qf1+q1 *rf1 # q11
         r1    = r1*rf1-q1 *qf1 # r11
         q1    = r2    +t1      # q21
         r2    =(r2    -t1)*qep # r22
         t1    = r1    +q2      # r21
         r1    =(r1    -q2)*qep # q22
         r[k0] = r[k0]+t1       # r40
         q[k0] = q[k0]+q1       # q40
         t1    = r[k0]-t1 *rec  # r21
         q1    = q[k0]-q1 *rec  # q21
         q[k2] = q1    -r1      # q42
         q[k1] = q1    +r1      # q41
         r[k1] = t1    -r2      # r41
         r[k2] = t1    +r2      # r42
      end
   end
   ma=ma3
end

if (ln5 > 0)
#  RADIX 5
   nze=Int(floor(n/5));  rze=w[nze]; qze=w[nh+nze]; rzc=1.0-rze # <- zeta
   ret=rze*rze-qze*qze;  qet=2*rze*qze; rec=1.0-ret # <- eta
   for l=1:ln5
      mb=Int(floor(mb/5)); ma5=ma*5
      for j=1:ma
         jmb=j*mb;             jmb2=jmb*2
         rf1=w[jmb];           qf1=w[nh+jmb]       # f**1
         rf2=w[jmb2];          qf2=w[nh+jmb2]      # f**2
         rf3=rf1*rf2-qf1*qf2;  qf3=rf1*qf2+qf1*rf2 # f**3
         rf4=rf2*rf2-qf2*qf2;  qf4=2*rf2*qf2       # f**4
         for i=1:ma5:n
            k0=(i-1)+j;   k1=k0+ma;  k2=k1+ma;  k3=k2+ma;  k4=k3+ma
            r1=r[k1]; r2=r[k2]; r3=r[k3]; r4=r[k4]
            q1=q[k1]; q2=q[k2]; q3=q[k3]; q4=q[k4]
            t1    =r1*qf1+q1*rf1       # q11
            r1    =r1*rf1-q1*qf1       # r11
            q1    =r4*rf4-q4*qf4       # r14
            r4    =r4*qf4+q4*rf4       # q14
            q4    =r1    -q1           # q24
            r1    =r1    +q1           # r21
            q1    =t1    +r4           # q21
            r4    =t1    -r4           # q24
            t1    =r3*rf3-q3*qf3       # r13
            r3    =r3*qf3+q3*rf3       # q13
            q3    =r2*qf2+q2*rf2       # q12
            r2    =r2*rf2-q2*qf2       # r12
            q2    =q3    +r3           # q22
            r3    =q3    -r3           # q23
            q3    =r2    -t1           # r23
            r2    =r2    +t1           # r22
            r[k0] =r[k0]+r1    +r2     # r40
            q[k0] =q[k0]+q1    +q2     # q40
            t1    =r4*qze+r3*qet       # q34
            r3    =r4*qet-r3*qze       # q33
            r4    =r[k0]-r1*rec-r2*rzc # r32
            r1    =r[k0]-r1*rzc-r2*rec # r31
            r[k2] =r4    -r3           # r42
            r[k3] =r4    +r3           # r43
            r[k4] =r1    +t1           # r44
            r[k1] =r1    -t1           # r41
            t1    =q[k0]-q1*rzc-q2*rec # q31
            q2    =q[k0]-q1*rec-q2*rzc # q32
            q1    =q4*qet-q3*qze       # r33
            q4    =q4*qze+q3*qet       # r34
            q[k3] =q2    -q1           # q43
            q[k2] =q2    +q1           # q42
            q[k1] =t1    +q4           # q41
            q[k4] =t1    -q4           # q44
         end
      end
      ma=ma5
   end
end

if (ln7 > 0)
#  RADIX 7
   nze=Int(floor(n/7))
   rze=w[nze];          qze=w[nh+nze];       rzc=1.0-rze # <- zeta
   ret=rze*rze-qze*qze; qet=2*rze*qze;       rec=1.0-ret # <- eta
   rth=rze*ret-qze*qet; qth=rze*qet+qze*ret; rtc=1.0-rth # <- theta
   for L=1:Ln7
      mb=Int(floor(mb/7)); ma7=ma*7
      for j=1:ma
         jmb=j*mb; jmb2=jmb*2; jmb3=jmb*3
         rf1=w[jmb];          qf1=w[nh+jmb]       # f**1
         rf2=w[jmb2];         qf2=w[nh+jmb2]      # f**2
         rf3=w[jmb3];         qf3=w[nh+jmb3]      # f**3
         rf4=rf2*rf2-qf2*qf2; qf4=2*rf2*qf2       # f**4
         rf5=rf2*rf3-qf2*qf3; qf5=rf2*qf3+qf2*rf3 # f**5
         rf6=rf3*rf3-qf3*qf3; qf6=2*rf3*qf3       # f**6
         for i=1:ma7:n
            k0=(i-1)+j
            k1=k0+ma; k2=k1+ma; k3=k2+ma; k4=k3+ma; k5=k4+ma; k6=k5+ma
            r0=r[k0]
            r1=r[k1];r2=r[k2];r3=r[k3];r4=r[k4];r5=r[k5];r6=r[k6]
            q0=q[k0]
            q1=q[k1];q2=q[k2];q3=q[k3];q4=q[k4];q5=q[k5];q6=q[k6]
            t1=r1*qf1+q1*rf1           # q11
            r1=r1*rf1-q1*qf1           # r11
            q1=r6*rf6-q6*qf6           # r16
            r6=r6*qf6+q6*rf6           # q16
            q6=r1-q1                   # r26
            r1=r1+q1                   # r21
            q1=t1+r6                   # q21
            r6=t1-r6                   # q26
            t2=r2*qf2+q2*rf2           # q12
            r2=r2*rf2-q2*qf2           # r12
            q2=r5*rf5-q5*qf5           # r15
            r5=r5*qf5+q5*rf5           # q15
            q5=r2-q2                   # r25
            r2=r2+q2                   # r22
            q2=t2+r5                   # q22
            r5=t2-r5                   # q25
            t3=r3*qf3+q3*rf3           # q13
            r3=r3*rf3-q3*qf3           # r13
            q3=r4*rf4-q4*qf4           # r14
            r4=r4*qf4+q4*rf4           # q14
            q4=r3-q3                   # r24
            r3=r3+q3                   # r23
            q3=t3+r4                   # q23
            r4=t3-r4                   # q24
            r0=r0+r1+r2+r3             # r40
            q0=q0+q1+q2+q3             # q40
            t1=   r6*qze+r5*qet+r4*qth # q36
            t2=   r6*qet-r5*qth-r4*qze # q35
            t3=   r6*qth-r5*qze+r4*qet # q34
            r6=r0-r1*rzc-r2*rec-r3*rtc # r31
            r5=r0-r1*rec-r2*rtc-r3*rzc # r32
            r4=r0-r1*rtc-r2*rzc-r3*rec # r33
            r[k0]=r0                   # r40
            r[k1]=r6-t1                # r41
            r[k6]=r6+t1                # r46
            r[k2]=r5-t2                # r42
            r[k5]=r5+t2                # r45
            r[k3]=r4-t3                # r43
            r[k4]=r4+t3                # r44
            t1=   q6*qze+q5*qet+q4*qth # r36
            t2=   q6*qet-q5*qth-q4*qze # r35
            t3=   q6*qth-q5*qze+q4*qet # r34
            q6=q0-q1*rzc-q2*rec-q3*rtc # q31
            q5=q0-q1*rec-q2*rtc-q3*rzc # q32
            q4=q0-q1*rtc-q2*rzc-q3*rec # q33
            q[k0]=q0                   # q40
            q[k1]=q6+t1                # q41
            q[k6]=q6-t1                # q46
            q[k2]=q5+t2                # q42
            q[k5]=q5-t2                # q45
            q[k3]=q4+t3                # q43
            q[k4]=q4-t3                # q44
         end
      end
      ma=ma7
   end
end
if (lnb > 0)
#  RADIX 11
   nze=Int(floor(n/11));
   rze=w[nze];          qze=w[nh+nze];       rzc=1.0-rze # <- zeta
   ret=rze*rze-qze*qze; qet=2*rze*qze;       rec=1.0-ret # <- eta
   rth=rze*ret-qze*qet; qth=rze*qet+qze*ret; rtc=1.0-rth # <- theta
   rio=ret*ret-qet*qet; qio=2*ret*qet;       ric=1.0-rio # <- iota
   rka=ret*rth-qet*qth; qka=ret*qth+qet*rth; rkc=1.0-rka # <- kappa
   for L=1:Lnb
      mb=Int(floor(mb/11)); mab=ma*11
      for j=1:ma
         jmb=j*mb; jmb2=jmb*2; jmb3=jmb*3; jmb4=jmb*4; jmb5=jmb*5
         rf1=w[jmb];          qf1=w[nh+jmb]       # f**1
         rf2=w[jmb2];         qf2=w[nh+jmb2]      # f**2
         rf3=w[jmb3];         qf3=w[nh+jmb3]      # f**3
         rf4=w[jmb4];         qf4=w[nh+jmb4]      # f**4
         rf5=w[jmb5];         qf5=w[nh+jmb5]      # f**5
         rf6=rf3*rf3-qf3*qf3; qf6=2*rf3*qf3       # f**6
         rf7=rf3*rf4-qf3*qf4; qf7=rf3*qf4+qf3*rf4 # f**7
         rf8=rf4*rf4-qf4*qf4; qf8=2*rf4*qf4       # f**8
         rf9=rf4*rf5-qf4*qf5; qf9=rf4*qf5+qf4*rf5 # f**9
         rfa=rf5*rf5-qf5*qf5; qfa=2*rf5*qf5       # f**10
         for i=1:mab:n
            k0=(i-1)+j
            k1=k0+ma; k2=k1+ma; k3=k2+ma; k4=k3+ma; k5=k4+ma; k6=k5+ma
            k7=k6+ma; k8=k7+ma; k9=k8+ma; ka=k9+ma
            r0=r[k0]
            r1=r[k1];r2=r[k2];r3=r[k3];r4=r[k4];r5=r[k5];r6=r[k6]
            r7=r[k7];r8=r[k8];r9=r[k9];ra=r[ka]
            q0=q[k0]
            q1=q[k1];q2=q[k2];q3=q[k3];q4=q[k4];q5=q[k5];q6=q[k6]
            q7=q[k7];q8=q[k8];q9=q[k9];qa=q[ka]
            t1=r1*qf1+q1*rf1           # q11
            r1=r1*rf1-q1*qf1           # r11
            q1=ra*rfa-qa*qfa           # r1a
            ra=ra*qfa+qa*rfa           # q1a
            qa=r1-q1                   # r2a
            r1=r1+q1                   # r21
            q1=t1+ra                   # q21
            ra=t1-ra                   # q2a
            t2=r2*qf2+q2*rf2           # q12
            r2=r2*rf2-q2*qf2           # r12
            q2=r9*rf9-q9*qf9           # r19
            r9=r9*qf9+q9*rf9           # q19
            q9=r2-q2                   # r29
            r2=r2+q2                   # r22
            q2=t2+r9                   # q22
            r9=t2-r9                   # q29
            t3=r3*qf3+q3*rf3           # q13
            r3=r3*rf3-q3*qf3           # r13
            q3=r8*rf8-q8*qf8           # r18
            r8=r8*qf8+q8*rf8           # q18
            q8=r3-q3                   # r28
            r3=r3+q3                   # r23
            q3=t3+r8                   # q23
            r8=t3-r8                   # q28
            t4=r4*qf4+q4*rf4           # q14
            r4=r4*rf4-q4*qf4           # r14
            q4=r7*rf7-q7*qf7           # r17
            r7=r7*qf7+q7*rf7           # q17
            q7=r4-q4                   # r27
            r4=r4+q4                   # r24
            q4=t4+r7                   # q24
            r7=t4-r7                   # q27
            t5=r5*qf5+q5*rf5           # q15
            r5=r5*rf5-q5*qf5           # r15
            q5=r6*rf6-q6*qf6           # r16
            r6=r6*qf6+q6*rf6           # q16
            q6=r5-q5                   # r26
            r5=r5+q5                   # r25
            q5=t5+r6                   # q25
            r6=t5-r6                   # q26
            r0=r0+r1+r2+r3+r4+r5       # r40
            q0=q0+q1+q2+q3+q4+q5       # q40
            t1=   ra*qze+r9*qet+r8*qth+r7*qio+r6*qka # q3a
            t2=   ra*qet+r9*qio-r8*qka-r7*qth-r6*qze # q39
            t3=   ra*qth-r9*qka-r8*qet+r7*qze+r6*qio # q38
            t4=   ra*qio-r9*qth+r8*qze+r7*qka-r6*qet # q37
            t5=   ra*qka-r9*qze+r8*qio-r7*qet+r6*qth # q36
            ra=r0-r1*rzc-r2*rec-r3*rtc-r4*ric-r5*rkc # r31
            r9=r0-r1*rec-r2*ric-r3*rkc-r4*rtc-r5*rzc # r32
            r8=r0-r1*rtc-r2*rkc-r3*rec-r4*rzc-r5*ric # r33
            r7=r0-r1*ric-r2*rtc-r3*rzc-r4*rkc-r5*rec # r34
            r6=r0-r1*rkc-r2*rzc-r3*ric-r4*rec-r5*rtc # r35
            r[k0]=r0                   # r40
            r[k1]=ra-t1                # r41
            r[ka]=ra+t1                # r4a
            r[k2]=r9-t2                # r42
            r[k9]=r9+t2                # r49
            r[k3]=r8-t3                # r43
            r[k8]=r8+t3                # r48
            r[k4]=r7-t4                # r44
            r[k7]=r7+t4                # r47
            r[k5]=r6-t5                # r45
            r[k6]=r6+t5                # r46
            t1=   qa*qze+q9*qet+q8*qth+q7*qio+q6*qka # r3a
            t2=   qa*qet+q9*qio-q8*qka-q7*qth-q6*qze # r39
            t3=   qa*qth-q9*qka-q8*qet+q7*qze+q6*qio # r38
            t4=   qa*qio-q9*qth+q8*qze+q7*qka-q6*qet # r37
            t5=   qa*qka-q9*qze+q8*qio-q7*qet+q6*qth # r36
            qa=q0-q1*rzc-q2*rec-q3*rtc-q4*ric-q5*rkc # q31
            q9=q0-q1*rec-q2*ric-q3*rkc-q4*rtc-q5*rzc # q32
            q8=q0-q1*rtc-q2*rkc-q3*rec-q4*rzc-q5*ric # q33
            q7=q0-q1*ric-q2*rtc-q3*rzc-q4*rkc-q5*rec # q34
            q6=q0-q1*rkc-q2*rzc-q3*ric-q4*rec-q5*rtc # q35
            q[k0]=q0                   # q40
            q[k1]=qa+t1                # q41
            q[ka]=qa-t1                # q4a
            q[k2]=q9+t2                # q42
            q[k9]=q9-t2                # q49
            q[k3]=q8+t3                # q43
            q[k8]=q8-t3                # q48
            q[k4]=q7+t4                # q44
            q[k7]=q7-t4                # q47
            q[k5]=q6+t5                # q45
            q[k6]=q6-t5                # q46
         end
      end
      ma=mab
   end
end
if (lnd > 0)
#  RADIX 13
   nze=Int(floor(n/13))
   rze=w[nze];          qze=w[nh+nze];       rzc=1.0-rze # <- zeta
   ret=rze*rze-qze*qze; qet=2*rze*qze;       rec=1.0-ret # <- eta
   rth=rze*ret-qze*qet; qth=rze*qet+qze*ret; rtc=1.0-rth # <- theta
   rio=ret*ret-qet*qet; qio=2*ret*qet;       ric=1.0-rio # <- iota
   rka=ret*rth-qet*qth; qka=ret*qth+qet*rth; rkc=1.0-rka # <- kappa
   rla=rth*rth-qth*qth; qla=2*rth*qth;       rlc=1.0-rla # <- lambda
   for L=1:Lnd
      mb=Int(floor(mb/13)); mad=ma*13
      for j=1:ma
         jmb=j*mb; jmb2=jmb*2; jmb3=jmb*3; jmb4=jmb*4; jmb5=jmb*5; jmb6=jmb*6
         rf1=w[jmb];          qf1=w[nh+jmb]       # f**1
         rf2=w[jmb2];         qf2=w[nh+jmb2]      # f**2
         rf3=w[jmb3];         qf3=w[nh+jmb3]      # f**3
         rf4=w[jmb4];         qf4=w[nh+jmb4]      # f**4
         rf5=w[jmb5];         qf5=w[nh+jmb5]      # f**5
         rf6=w[jmb6];         qf6=w[nh+jmb6]      # f**6
         rf7=rf3*rf4-qf3*qf4; qf7=rf3*qf4+qf3*rf4 # f**7
         rf8=rf4*rf4-qf4*qf4; qf8=2*rf4*qf4       # f**8
         rf9=rf4*rf5-qf4*qf5; qf9=rf4*qf5+qf4*rf5 # f**9
         rfa=rf5*rf5-qf5*qf5; qfa=2*rf5*qf5       # f**10
         rfb=rf5*rf6-qf5*qf6; qfb=rf5*qf6+qf5*rf6 # f**11
         rfc=rf6*rf6-qf6*qf6; qfc=2*rf6*qf6       # f**12
         for i=1:mad:n
            k0=(i-1)+j
            k1=k0+ma; k2=k1+ma; k3=k2+ma; k4=k3+ma; k5=k4+ma; k6=k5+ma
            k7=k6+ma; k8=k7+ma; k9=k8+ma; ka=k9+ma; kb=ka+ma; kc=kb+ma
            r0=r[k0]
            r1=r[k1];r2=r[k2];r3=r[k3];r4=r[k4];r5=r[k5];r6=r[k6]
            r7=r[k7];r8=r[k8];r9=r[k9];ra=r[ka];rb=r[kb];rc=r[kc]
            q0=q[k0]
            q1=q[k1];q2=q[k2];q3=q[k3];q4=q[k4];q5=q[k5];q6=q[k6]
            q7=q[k7];q8=q[k8];q9=q[k9];qa=q[ka];qb=q[kb];qc=q[kc]
            t1=r1*qf1+q1*rf1           # q11
            r1=r1*rf1-q1*qf1           # r11
            q1=rc*rfc-qc*qfc           # r1c
            rc=rc*qfc+qc*rfc           # q1c
            qc=r1-q1                   # r2c
            r1=r1+q1                   # r21
            q1=t1+rc                   # q21
            rc=t1-rc                   # q2c
            t2=r2*qf2+q2*rf2           # q12
            r2=r2*rf2-q2*qf2           # r12
            q2=rb*rfb-qb*qfb           # r1b
            rb=rb*qfb+qb*rfb           # q1b
            qb=r2-q2                   # r2b
            r2=r2+q2                   # r22
            q2=t2+rb                   # q22
            rb=t2-rb                   # q2b
            t3=r3*qf3+q3*rf3           # q13
            r3=r3*rf3-q3*qf3           # r13
            q3=ra*rfa-qa*qfa           # r1a
            ra=ra*qfa+qa*rfa           # q1a
            qa=r3-q3                   # r2a
            r3=r3+q3                   # r23
            q3=t3+ra                   # q23
            ra=t3-ra                   # q2a
            t4=r4*qf4+q4*rf4           # q14
            r4=r4*rf4-q4*qf4           # r14
            q4=r9*rf9-q9*qf9           # r19
            r9=r9*qf9+q9*rf9           # q19
            q9=r4-q4                   # r29
            r4=r4+q4                   # r24
            q4=t4+r9                   # q24
            r9=t4-r9                   # q29
            t5=r5*qf5+q5*rf5           # q15
            r5=r5*rf5-q5*qf5           # r15
            q5=r8*rf8-q8*qf8           # r18
            r8=r8*qf8+q8*rf8           # q18
            q8=r5-q5                   # r28
            r5=r5+q5                   # r25
            q5=t5+r8                   # q25
            r8=t5-r8                   # q28
            t6=r6*qf6+q6*rf6           # q16
            r6=r6*rf6-q6*qf6           # r16
            q6=r7*rf7-q7*qf7           # r17
            r7=r7*qf7+q7*rf7           # q17
            q7=r6-q6                   # r27
            r6=r6+q6                   # r26
            q6=t6+r7                   # q26
            r7=t6-r7                   # q27
            r0=r0+r1+r2+r3+r4+r5+r6    # r40
            q0=q0+q1+q2+q3+q4+q5+q6    # q40
            t1=   rc*qze+rb*qet+ra*qth+r9*qio+r8*qka+r7*qla # q3c
            t2=   rc*qet+rb*qio+ra*qla-r9*qka-r8*qth-r7*qze # q3b
            t3=   rc*qth+rb*qla-ra*qio-r9*qze+r8*qet+r7*qka # q3a
            t4=   rc*qio-rb*qka-ra*qze+r9*qth-r8*qla-r7*qet # q39
            t5=   rc*qka-rb*qth+ra*qet-r9*qla-r8*qze+r7*qio # q38
            t6=   rc*qla-rb*qze+ra*qka-r9*qet+r8*qio-r7*qth # q37
            rc=r0-r1*rzc-r2*rec-r3*rtc-r4*ric-r5*rkc-r6*rlc # r31
            rb=r0-r1*rec-r2*ric-r3*rlc-r4*rkc-r5*rtc-r6*rzc # r32
            ra=r0-r1*rtc-r2*rlc-r3*ric-r4*rzc-r5*rec-r6*rkc # r33
            r9=r0-r1*ric-r2*rkc-r3*rzc-r4*rtc-r5*rlc-r6*rec # r34
            r8=r0-r1*rkc-r2*rtc-r3*rec-r4*rlc-r5*rzc-r6*ric # r35
            r7=r0-r1*rlc-r2*rzc-r3*rkc-r4*rec-r5*ric-r6*rtc # r36
            r[k0]=r0                   # r40
            r[k1]=rc-t1                # r41
            r[kc]=rc+t1                # r4c
            r[k2]=rb-t2                # r42
            r[kb]=rb+t2                # r4b
            r[k3]=ra-t3                # r43
            r[ka]=ra+t3                # r4a
            r[k4]=r9-t4                # r44
            r[k9]=r9+t4                # r49
            r[k5]=r8-t5                # r45
            r[k8]=r8+t5                # r48
            r[k6]=r7-t6                # r46
            r[k7]=r7+t6                # r47
            t1=   qc*qze+qb*qet+qa*qth+q9*qio+q8*qka+q7*qla # r3c
            t2=   qc*qet+qb*qio+qa*qla-q9*qka-q8*qth-q7*qze # r3b
            t3=   qc*qth+qb*qla-qa*qio-q9*qze+q8*qet+q7*qka # r3a
            t4=   qc*qio-qb*qka-qa*qze+q9*qth-q8*qla-q7*qet # r39
            t5=   qc*qka-qb*qth+qa*qet-q9*qla-q8*qze+q7*qio # r38
            t6=   qc*qla-qb*qze+qa*qka-q9*qet+q8*qio-q7*qth # r37
            qc=q0-q1*rzc-q2*rec-q3*rtc-q4*ric-q5*rkc-q6*rlc # q31
            qb=q0-q1*rec-q2*ric-q3*rlc-q4*rkc-q5*rtc-q6*rzc # q32
            qa=q0-q1*rtc-q2*rlc-q3*ric-q4*rzc-q5*rec-q6*rkc # q33
            q9=q0-q1*ric-q2*rkc-q3*rzc-q4*rtc-q5*rlc-q6*rec # q34
            q8=q0-q1*rkc-q2*rtc-q3*rec-q4*rlc-q5*rzc-q6*ric # q35
            q7=q0-q1*rlc-q2*rzc-q3*rkc-q4*rec-q5*ric-q6*rtc # q36
            q[k0]=q0                   # q40
            q[k1]=qc+t1                # q41
            q[kc]=qc-t1                # q4c
            q[k2]=qb+t2                # q42
            q[kb]=qb-t2                # q4b
            q[k3]=qa+t3                # q43
            q[ka]=qa-t3                # q4a
            q[k4]=q9+t4                # q44
            q[k9]=q9-t4                # q49
            q[k5]=q8+t5                # q45
            q[k8]=q8-t5                # q48
            q[k6]=q7+t6                # q46
            q[k7]=q7-t6                # q47
         end
      end
      ma=mad
   end
end
return r,q;
end
