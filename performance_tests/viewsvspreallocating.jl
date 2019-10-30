using BenchmarkTools

vector=rand(1000000+5)

testview = @benchmark begin
    sum1=0.0;
    for i in 1:1000000
        a= @views vector[i:i+5];
        #sum1+=sum(a);
    end
end


testpreallocatearray = @benchmark begin
    sum1=0.0
    a=ones(6);
    for i in 1:1000000
        a=vector[i:i+5];
        #sum1+=sum(a);
    end
end

testpreallocateview = @benchmark begin
    sum1=0.0
    a = @views vector[1:6]
    for i in 1:1000000
        a= @views vector[i:i+5];
        #sum1+=sum(a);
    end
end

println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
println("Old Version: creating a new view in every loop-iteration\n")
display(testview);
println("\n\n");

println("New try #1: preallocating an array, not using @views\n")
display(testpreallocatearray);
println("\n\n");

println("New try #2: preallocating a view\n")
display(testpreallocateview);
println("\n\n");
