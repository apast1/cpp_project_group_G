
% data = rand(2,4);
% data= [95,0.2,0.06,0.5;105,0.4,0.1,1];

function test_matrix = trans4_16(data)


i=1;
for a=1:2
    for b=1:2
        for c=1:2
            for d=1:2
                    temp = [data(a,1),data(b,2),data(c,3),data(d,4)];
                    mat(i,:)=temp;
                    i=i+1;
            end
        end
    end
end
    
test_matrix=mat;