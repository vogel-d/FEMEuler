struct MIS #Multirate Inﬁnitesimal Step method
    nStage::Int64;

    #Die Koeffizientenmarizen von alpha, beta, gamma haben jeweils die Dimension nStage+1 x nStage
    alpha::Array{Float64,2};
    beta::Array{Float64,2};
    gamma::Array{Float64,2};

    #c und d haben die Länge nStage+1
    c::Array{Float64,1};
    d::Array{Float64,1};

    betaS::Array{Float64,2};
    gammaS::Array{Float64,2};
    nPhiStage::Array{Int64,1};
    nPhi::Int64
    nPhiStageBeta::Array{Int64,1};
    nPhiStageGamma::Array{Int64,1};
end

function MIS(s::Symbol)
    if s==:MIS5_4
      nStage=5;
      nPhiStage=zeros(Int64,nStage+1);
      nPhi=1;

      nPhiStageBeta=ones(Int64,nStage+1);
      nPhiStageBeta[1]= 0;

      nPhiStageGamma=ones(Int64,nStage+1);
      nPhiStageGamma[1]= 0;

      beta=zeros(nStage+1, nStage);
      beta[2,1]=  0.529563257126527;
      beta[3,1]=  0.209002691927241;
      beta[3,2]= -0.172349890927723;
      beta[4,1]= -0.314103776563258;
      beta[4,2]=  0.832388085305425;
      beta[4,3]= -0.518284308724052;
      beta[5,1]= -1.180526334435395;
      beta[5,2]=  1.077092096618435;
      beta[5,3]=  0.519427132306836;
      beta[5,4]=  0.395577371153868;
      beta[6,1]= -0.280934261471696;
      beta[6,2]=  0.419312573639456;
      beta[6,3]= -0.257246871763194;
      beta[6,4]= -0.092439012970596;
      beta[6,5]=  0.399033266416755;

      c=zeros(nStage+1);
      c[2]=  0.529563257126527;
      c[3]=  0.036610765667902;
      c[4]=  0.027460619651280;
      c[5]=  0.811935172198554;
      c[6]=  0.999999999999404;

      alpha=zeros(nStage+1,nStage);
      alpha[3,2]=   0.001067565501309;
      alpha[4,2]=  -0.442781118636665;
      alpha[4,3]=   0.991850184207169;
      alpha[5,2]=   0.000145652123851;
      alpha[5,3]=   0.230748323989429;
      alpha[5,4]=  -0.272045335417296;
      alpha[6,2]=  -0.003334820890679;
      alpha[6,3]=   0.507413848225913;
      alpha[6,4]=  -0.673178314864792;
      alpha[6,5]=   1.002677604378304;

      gamma=zeros(nStage+1, nStage);
      gamma[3,2]=  -0.001146942858123;
      gamma[4,2]=   0.441818586053289;
      gamma[4,3]=  -0.227857926388024;
      gamma[5,2]=  -0.002816964038452;
      gamma[5,3]=   0.189471581599903;
      gamma[5,4]=  -0.223392767297830;
      gamma[6,2]=   0.001973069409795;
      gamma[6,3]=   0.251301260882666;
      gamma[6,4]=  -0.328491753699942;
      gamma[6,5]=  -0.001705104076025;

      d=zeros(nStage+1);
      for i=2:nStage+1
         for j=1:nStage
            d[i]=d[i]+beta[i,j];
         end
      end

   elseif s==:MIS4_4
      nStage=4;
      nPhiStage=zeros(Int64,nStage+1);
      nPhi=1;

      nPhiStageBeta=ones(Int64,nStage+1);
      nPhiStageBeta[1]= 0;

      nPhiStageGamma=ones(Int64,nStage+1);
      nPhiStageGamma[1]= 0;

      beta=zeros(nStage+1, nStage);
      beta[2,1]=  0.38758444641450318;
      beta[3,1]=  -2.5318448354142823E-002;
      beta[3,2]=  0.38668943087310403;
      beta[4,1]=  0.20899983523553325;
      beta[4,2]= -0.45856648476371231;
      beta[4,3]=  0.43423187573425748;
      beta[5,1]= -0.10048822195663100;
      beta[5,2]= -0.46186171956333327;
      beta[5,3]=  0.83045062122462809;
      beta[5,4]=  0.27014914900250392;

      c=zeros(nStage+1);
      c[2]=  0.38758444641450318;
      c[3]=  0.61521685655017821;
      c[4]=  0.23254717315441453;
      c[5]=   1.0000000000000002;

      alpha=zeros(nStage+1,nStage);
      alpha[3,2]=  0.52349249922385610;
      alpha[4,2]=   1.1683374366893629;
      alpha[4,3]= -0.75762080241712637;
      alpha[5,2]=  -3.6477233846797109E-002;
      alpha[5,3]=  0.56936148730740477;
      alpha[5,4]=  0.47746263002599681;

      gamma=zeros(nStage+1, nStage);
      gamma[3,2]=  0.13145089796226542;
      gamma[4,2]= -0.36855857648747881;
      gamma[4,3]=  0.33159232636600550;
      gamma[5,2]=  -6.5767130537473045E-002;
      gamma[5,3]=   4.0591093109036858E-002;
      gamma[5,4]=   6.4902111640806712E-002;

      d=zeros(nStage+1);
      for i=2:nStage+1
         for j=1:nStage
            d[i]=d[i]+beta[i,j];
         end
      end

   elseif s==:MIS4
      nStage=4;
      nPhiStage=zeros(Int64,nStage+1);
      nPhi=1;

      nPhiStageBeta=ones(Int64,nStage+1);
      nPhiStageBeta[1]= 0;

      nPhiStageGamma=ones(Int64,nStage+1);
      nPhiStageGamma[1]= 0;

      beta=zeros(nStage+1, nStage);
      beta[2,1]=  0.13629647842266179;
      beta[3,1]=  0.28046239897933556;
      beta[3,2]= -1.60351333596248577E-002;
      beta[4,1]=  0.90471335520843155;
      beta[4,2]=  -1.0401118315403626;
      beta[4,3]=  0.65233756348866223;
      beta[5,1]=  6.71969845545695721E-002;
      beta[5,2]= -0.36562186260961194;
      beta[5,3]= -0.15486147083521187;
      beta[5,4]=  0.97036244446880304;


      c=zeros(nStage+1);
      c[2]=  0.38758444641450318;
      c[3]=  0.61521685655017821;
      c[4]=  0.23254717315441453;
      c[5]=   1.0000000000000002;

      alpha=zeros(nStage+1,nStage);
      alpha[3,2]=  0.52349249922385610;
      alpha[4,2]=   1.1683374366893629;
      alpha[4,3]= -0.75762080241712637;
      alpha[5,2]=  -3.6477233846797109E-002;
      alpha[5,3]=  0.56936148730740477;
      alpha[5,4]=  0.47746263002599681;

      gamma=zeros(nStage+1, nStage);
      gamma[3,2]=  0.13145089796226542;
      gamma[4,2]= -0.36855857648747881;
      gamma[4,3]=  0.33159232636600550;
      gamma[5,2]=  -6.5767130537473045E-002;
      gamma[5,3]=   4.0591093109036858E-002;
      gamma[5,4]=   6.4902111640806712E-002;

      d=zeros(nStage+1);
      for i=2:nStage+1
         for j=1:nStage
            d[i]=d[i]+beta[i,j];
         end
      end

   elseif s==:MIS2
       nStage=3;
       nPhiStage=zeros(nStage+1);
       nPhi=1;

       nPhiStageBeta=ones(Int64,nStage+1);
       nPhiStageBeta[1]= 0;

       nPhiStageGamma=ones(Int64,nStage+1);
       nPhiStageGamma[1]= 0;

       beta=zeros(nStage+1, nStage);
       beta[2,1]=  0.12684849455255601;
       beta[3,1]= -0.78483827882640156;
       beta[3,2]=   1.3744267526826737;
       beta[4,1]= -4.56727081748555391e-002;
       beta[4,2]= -8.75082271190387971e-003;
       beta[4,3]=  0.52477578862897312;

       c=zeros(nStage+1);
       c[2]=1/3;
       c[3]=1/6;
       c[4]=1.0;

       alpha=zeros(nStage+1,nStage);
       alpha[3,2]=  0.53694656671020691;
       alpha[4,2]=  0.48089296855085184;
       alpha[4,3]=  0.50056116356635882;

       gamma=zeros(nStage+1, nStage);
       gamma[3,2]=  0.65246512600423157;
       gamma[4,2]= -7.32769849456572780E-002;
       gamma[4,3]=  0.14490243042028150;

       d=zeros(nStage+1);
       for i in 2:nStage+1
         for j=1:nStage
           d[i]=d[i]+beta[i,j];
         end
       end
   elseif s==:MIS_Euler

      nStage=1;
      nPhiStage=zeros(nStage+1);
      nPhi=1;

      nPhiStageBeta=ones(Int64,nStage+1);
      nPhiStageBeta[1]= 0;

      nPhiStageGamma=ones(Int64,nStage+1);
      nPhiStageGamma[1]= 0;

      beta=zeros(nStage+1, nStage);
      beta[2,1] = 1.0;

      c=zeros(nStage+1);
      c[2]=1.0;

      alpha=zeros(nStage+1,nStage);

      gamma=zeros(nStage+1, nStage);

      d=zeros(nStage+1);
      for i in 2:nStage+1
        for j=1:nStage
          d[i]=d[i]+beta[i,j];
        end
      end

   elseif s==:MISRK2

      nStage=2;
      nPhiStage=zeros(nStage+1);
      nPhi=1;

      nPhiStageBeta=zeros(Int64,nStage+1);
      nPhiStageBeta[2]= 0;
      nPhiStageBeta[3]= 0;

      nPhiStageGamma=zeros(Int64,nStage+1);
      nPhiStageGamma[2]= 0;
      nPhiStageGamma[3]= 0;

      beta=zeros(nStage+1, nStage);
      beta[2,1] = 0.5;
      beta[3,2] = 1.0;

      c=zeros(nStage+1);
      c[2]=0.5;
      c[3]=1.0;

      alpha=zeros(nStage+1,nStage);

      gamma=zeros(nStage+1, nStage);

      d=c;

   else
      error("Keine implementierte MIS-Methode.")
   end

   betaS=copy(beta);
   gammaS=copy(gamma);
   for i in 2:(nStage+1)
      for j in 1:(i-1)
         kFac=1;
         for k in 1:nPhi
            kFac=kFac*max(k-1,1);
            betaS[i,j]=beta[i,j]/(kFac*d[i]);
            beta[i,j]=beta[i,j]/d[i];
            gammaS[i,j]=gamma[i,j]/(kFac*d[i]);
            gamma[i,j]=gamma[i,j]/d[i];
         end
      end
   end
   MIS(nStage, alpha,beta, gamma, c, d, betaS, gammaS, nPhiStage, nPhi, nPhiStageBeta, nPhiStageGamma)
end
