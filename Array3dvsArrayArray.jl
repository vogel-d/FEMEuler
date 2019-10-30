A1=Array{Array{Float64,2},1}(undef,10)
B1=Array{Array{Float64,2},1}(undef,10)

A2=Array{Float64,3}(undef,10,1000,1000)
B2=Array{Float64,3}(undef,10,1000,1000)

for i in 1:10
    A1[i]=rand(1000,1000)
    B1[i]=rand(1000,1000)

    A2[i,:,:]=A1[i]
    B2[i,:,:]=B1[i]
end



function arrayarray(A::Array{Array{Float64,2},1}, B::Array{Array{Float64,2},1})
    sum=0.0;
    for i in 1:10
        for j in 1:10
            for k=1:1000
                for l=1:1000
                    sum+=A[i][l,k]*B[j][l,k]
                end
            end
        end
    end
 return sum;
end

function array3d(A::Array{Float64,3}, B::Array{Float64,3})
    sum=0.0;
    for k in 1:1000
        for l in 1:1000
            for i=1:10
                for j=1:10
                    sum+=A[i,l,k]*B[j,l,k]
                end
            end
        end
    end
 return sum;
end

arrayarray(A1,B1);
array3d(A2,B2);

sum1=0.0;
testarrayarray = @benchmark sum1=arrayarray(A1,B1);

sum2=0.0;
test3d= @benchmark sum2=array3d(A2,B2);

#sum3=A2*(B2')
println(isequal(sum1,sum2));
#println("test, ", isequal(sum1,sum3));
display(test3d);
display(testarrayarray);
