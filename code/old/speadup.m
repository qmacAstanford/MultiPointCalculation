
v1=rand(1,5);
v2=rand(5,1);

ntimes=200000;

tic()

for j=1:ntimes
    
    out=v2*v1;
    
end
toc()
tic()
for j=1:ntimes
    
    out=kron(v2,v1);
    
end

toc()