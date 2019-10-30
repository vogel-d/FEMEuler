using BenchmarkTools

#meantime: 1.903ms
#allocs estimate: 2
#memory estimate: 7.63MiB
function testview()
    vector=rand(1000000+5)
    sum1=0.0;

    for i in 1:1000000
        a= @view vector[i:i+5];
        #sum1+=sum(a);
    end
end

#meantime: 9.429ms
#allocs estimate: 3
#memory estimate: 7.63MiB
function testpreallocatearray()
    vector=rand(1000000+5)
    sum1=0.0
    a=ones(6);

    for i in 1:1000000
        for j in 0:5
            a[j+1]=vector[i+j]
        end
        #a=vector[i:i+5];
        #sum1+=sum(a);
    end
end

#meantime: 1.871ms
#allocs estimate: 2
#memory estimate: 7.63MiB
function testpreallocateview()
    vector=rand(1000000+5)
    sum1=0.0
    a = @views vector[1:6]
    
    for i in 1:1000000
        a= @views vector[i:i+5];
        #sum1+=sum(a);
    end
end
