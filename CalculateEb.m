function [Eb] = CalculateEb(x,size,E_S,parameter)

S=zeros(6*size*size,1);
field=parameter(14);

for i=1:1:size
    for j=1:1:size
        
        index_x=((i-1)*size+j-1)*4;
        index_S=((i-1)*size+j-1)*6;

        theta1=x(index_x+1);
        phi1=x(index_x+2);
        theta2=x(index_x+3);
        phi2=x(index_x+4);

        S(index_S+1:index_S+6)=[sin(theta1)*cos(phi1);sin(theta1)*sin(phi1);cos(theta1);sin(theta2)*cos(phi2);sin(theta2)*sin(phi2);cos(theta2)];
    end
end

field_vec=zeros(6*size*size,1);

for i=1:1:size
    for j=1:1:size

        index=((i-1)*size+j-1)*6;
        field_vec(index+2)=field;
        field_vec(index+5)=field;

    end
end


Eb=dot(S,E_S*S)+dot(S,field_vec);


end

