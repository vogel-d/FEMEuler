function getKub(g::Int, mt::Int)
    if mt==4
        #Viereckseiten werden als parallel zu den Koordinatenachsen vorrausgesetzt
        x=[0.0,1.0];
        y=[0.0, 1.0];
        n=Int(ceil((g+1)/2));
        gdots=Array{AbstractFloat,1};
        W=Array{AbstractFloat,2};

        if g>9
            #a1=1/3*sqrt(5-2*sqrt(10/7));
            #a2=1/3*sqrt(5+2*sqrt(10/7));
            a1=0.538469310105683;
            a2=0.906179845938664;
            gdots=[-a2,-a1,0.0,a1,a2];
            b00=0.32363456790123457;
            b11=0.05613434886242864;
            b22=0.22908540422399112;
            b01=0.13478507238752090;
            b02=0.27228653255075070;
            b12=0.1134;
            W=[b11 b12 b01 b12 b11;
               b12 b22 b02 b22 b12;
               b01 b02 b00 b02 b01;
               b12 b22 b02 b22 b12;
               b11 b12 b01 b12 b11];
            @warn "Der eingegebene Genauigkeitsgrad ist nicht erreichbar. Dies ist eine Approximation vom Genauigkeitsgrad 9."
        elseif g>7
            #a1=1/3*sqrt(5-2*sqrt(10/7));
            #a2=1/3*sqrt(5+2*sqrt(10/7));
            a1=0.538469310105683;
            a2=0.906179845938664;
            gdots=[-a2,-a1,0.0,a1,a2];
            b00=0.32363456790123457;
            b11=0.05613434886242864;
            b22=0.22908540422399112;
            b01=0.13478507238752090;
            b02=0.27228653255075070;
            b12=0.1134;
            W=[b11 b12 b01 b12 b11;
               b12 b22 b02 b22 b12;
               b01 b02 b00 b02 b01;
               b12 b22 b02 b22 b12;
               b11 b12 b01 b12 b11];
        elseif g>5
            #a1=sqrt(3/7+(2/7)*sqrt(6/5));
            #a2=sqrt(3/7-(2/7)*sqrt(6/5));
            a1=0.8611363115940526;
            a2=0.3399810435848563;
            gdots=[-a1,-a2,a2,a1];
            #b1=354-36*sqrt(30);
            #b2=354+36*sqrt(30);
            b1=156.8198792981401992;
            b2=551.1801207018598008;
            W=(1/1296)*[b1 294 294 b1;
                        294 b2 b2 294;
                        294 b2 b2 294;
                        b1 294 294 b1;];
        elseif g>3
            #a=sqrt(3/5);
            a=0.7745966692414834;
            gdots=[-a,0.0,a];
            W=(1/81)*[25 40 25;
                       40 64 40;
                       25 40 25];
        elseif g>1
            #a=1/sqrt(3);
            a=0.5773502691896258;
            gdots=[-a,a];
            W=[1.0 1.0; 1.0 1.0];
        else
            error("Geben Sie ein g>0 an.")
        end
        hx=(x[2]-x[1])/2;
        hy=(y[2]-y[1])/2;
        sx=(x[2]+x[1])/2;
        sy=(y[2]+y[1])/2;

        kubPoints=Array{AbstractFloat,2}(undef,2,n);
        for k in 1:n
            kubPoints[1,k]=hx*gdots[k]+sx;
            kubPoints[2,k]=hy*gdots[k]+sy;
        end

        return kubPoints, hx*hy*W
    elseif mt==3
        if g<3
            n=3;
        elseif g>=3
            n=7;
        end
        kubPoints=Array{AbstractFloat,2}(undef,2,n);
        W=Array{AbstractFloat,1};

        if g<=2
            kubPoints=[0.5 0.5 0.0;
                       0.0 0.5 0.5];

            W=[0.16666666666666666 0.16666666666666666 0.16666666666666666];
        elseif g>2
            kubPoints=[0.0 1.0 0.0 0.5 0.5 0.0 0.3333333333333333;
                       0.0 0.0 1.0 0.0 0.5 0.5 0.3333333333333333];

            a=0.06666666666666667;
            W=[0.025 0.025 0.025 a a a 0.225]

        elseif g>3
            @warn "Der eingegebene Genauigkeitsgrad ist nicht erreichbar. Dies ist eine Approximation vom Genauigkeitsgrad 3."
            kubPoints=[0.0 1.0 0.0 0.5 0.5 0.0 0.3333333333333333;
                       0.0 0.0 1.0 0.0 0.5 0.5 0.3333333333333333];

           a=0.06666666666666667;
           W=[0.025, 0.025, 0.025, a, a, a, 0.225]
        end

        return kubPoints, W;
    end
end
