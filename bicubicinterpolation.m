%function a00 = bicubicinterpolation(x,y,f)
A=zeros(16); fh = reshape(f',16,1); [xx yy] =meshgrid(x,y);
vx = reshape(xx',16,1);
vy = reshape(yy',16,1);

for j = 1:4
    f(j)= 0;
    for i = 1:4
        tmp = 0;
        for k = 1:i
            tmp = tmp + f(k,j); 
        end
        f(j) = f(j) + f(i,j)
        
    end
end
s1 = 49/144;s2 = 7/144;s3=1/144;

s = s1*(f(2,2)+f(3,2)+f(2,3)+f(3,3));
s = s - s2* (f(1,2)+f(1,3)+f(4,2)+f(4,3)+f(2,1)+f(3,1)+f(2,4)+f(3,4));
s = s + s3* (f(1,4)+f(4,1)+f(1,1)+f(4,4));
