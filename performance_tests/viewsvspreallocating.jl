using BenchmarkTools

function testest(a,b)
    return a*b;
end

#meantime: 547.947 μs
#allocs estimate: 0
#memory estimate: 0
function testview(vector)
    sum1=0.0;

    for i in 1:1000000
        a= @view vector[i:i+5];
        #sum1+=sum(a);
    end
end

#meantime: 15.206ns ########### WINNER
#allocs estimate: 0
#memory estimate: 0
function testviewinbounds(vector)
    sum1=0.0;

    for i in 1:1000000
        a= @inbounds @views vector[i:i+5];
        #sum1+=sum(a);
    end
end

#meantime: 6.004ms
#allocs estimate: 1
#memory estimate: 128 bytes
function testpreallocatearray(vector)
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

#meantime: 415.998 μs
#allocs estimate: 0
#memory estimate: 0 bytes
function testpreallocateview(vector)
    sum1=0.0
    a = @views vector[1:6]

    for i in 1:1000000
        a= @views vector[i:i+5];
        #sum1+=sum(a);
    end
end
