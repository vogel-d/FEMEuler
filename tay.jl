function tay(aco::Array{Float64,1},z::Complex{Float64}) #pfft1.f90 + subroutine lau aus pfft1.f90
   w=complex(0.0,0.0);
   zp=complex(1.0,0.0);
   for k in 1:length(aco)
      zp=zp*z
      w+=aco[k]*zp
   end
   return w;
end

function tay_reset(nco::Int,zzf::OffsetArray{Complex{Float64},1},r::Float64) #pfft1.f90
   # Reset the Taylor series coefficients, assumed real, from the results
   # of a complex fourier transformation of cycle nf. The nominal radius of
   # the circuit is r.

   #complex(dpc),dimension(0:nf-1),intent(INOUT):: zzf => nf=length(zzf)
   zzf=dfft(zzf)
   ri=1.0/r;
   rp=1.0;
   co=Array{Float64,1}(undef,nco);
   for i in 1:nco
      rp=rp*ri;
      co[i]=real(zzf[i])*rp;
   end
   return co;
end
