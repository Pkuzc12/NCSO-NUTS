function [gradEb] = GradEb(x,size,ES,parameter)

field=parameter(14);

S=zeros(6*size*size,1);

grad_field=zeros(6*size*size,1);

trans=zeros(4*size*size,6*size*size);

for i=1:1:size
    for j=1:1:size
        
        index_x=((i-1)*size+j-1)*4;
        index_S=((i-1)*size+j-1)*6;

        theta1=x(index_x+1);
        phi1=x(index_x+2);
        theta2=x(index_x+3);
        phi2=x(index_x+4);

        S(index_S+1:index_S+6)=[sin(theta1)*cos(phi1);sin(theta1)*sin(phi1);cos(theta1);sin(theta2)*cos(phi2);sin(theta2)*sin(phi2);cos(theta2)];
        trans(index_x+1:index_x+2,index_S+1:index_S+3)=[cos(theta1)*cos(phi1),cos(theta1)*sin(phi1),-sin(theta1);-sin(theta1)*sin(phi1),sin(theta1)*cos(phi1),0];
        trans(index_x+3:index_x+4,index_S+4:index_S+6)=[cos(theta2)*cos(phi2),cos(theta2)*sin(phi2),-sin(theta2);-sin(theta2)*sin(phi2),sin(theta2)*cos(phi2),0];
    end
end

for i=1:1:size
    for j=1:1:size

        index=((i-1)*size+j-1)*6;
        grad_field(index+2)=field;
        grad_field(index+5)=field;

    end
end

grad_S=2*ES*S+grad_field;

gradEb=trans*grad_S;




end

