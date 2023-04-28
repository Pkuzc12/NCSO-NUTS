function [ES] = GradEbHelper(chart,size,parameter)

ES=zeros(size*size*6);
[J1,J1p,J2,J2p,J3,K,G,Gs,Kp,Gp,Gps,B,D,field]=deal(parameter(1),parameter(2),parameter(3),parameter(4),parameter(5),parameter(6),parameter(7),parameter(8),parameter(9),parameter(10),parameter(11),parameter(12),parameter(13),parameter(14));

U=[[1/sqrt(6),-1/sqrt(2),1/sqrt(3)];[1/sqrt(6),1/sqrt(2),1/sqrt(3)];[-sqrt(2/3),0,1/sqrt(3)]];
J=[[0,G,Gs];[G,0,Gs];[Gs,Gs,K]];
Jp=[[0,Gp,Gps];[Gp,0,Gps];[Gps,Gps,Kp]];
order=[[[0,1,0];[0,0,1];[1,0,0]];[[0,0,1];[1,0,0];[0,1,0]];[[1,0,0];[0,1,0];[0,0,1]]];

for i=1:1:size
    for j=1:1:size
        index1=((i-1)*size+j-1)*6;
        ES(index1+2,index1+2)=B;
        ES(index1+3,index1+3)=D;
        ES(index1+5,index1+5)=B;
        ES(index1+6,index1+6)=D;


        for k=1:1:2
            nei=chart(i,j).left.neighbor(k,:);
            index2=((nei(1)-1)*size+nei(2)-1)*6;
            for t=1:1:3
                ES(index1+t,index2+t+3)=J1p/2;
            end
            nei=chart(i,j).right.neighbor(k,:);
            index2=((nei(1)-1)*size+nei(2)-1)*6;
            for t=1:1:3
                ES(index1+t+3,index2+t)=J1p/2;
            end
        end
        nei=chart(i,j).left.neighbor(3,:);
        index2=((nei(1)-1)*size+nei(2)-1)*6;
        for t=1:1:3
            ES(index1+t,index2+t+3)=J1/2;
        end
        nei=chart(i,j).right.neighbor(3,:);
        index2=((nei(1)-1)*size+nei(2)-1)*6;
        for t=1:1:3
            ES(index1+t+3,index2+t)=J1/2;
        end


        nei=chart(i,j).left.neighbor2(1,:);
        index2=((nei(1)-1)*size+nei(2)-1)*6;
        for t=1:1:3
            ES(index1+t,index2+t)=J2p/2;
        end
        
        nei=chart(i,j).right.neighbor2(1,:);
        index2=((nei(1)-1)*size+nei(2)-1)*6;
        for t=1:1:3
            ES(index1+t+3,index2+t+3)=J2p/2;
        end
        
        nei=chart(i,j).left.neighbor2(4,:);
        index2=((nei(1)-1)*size+nei(2)-1)*6;
        for t=1:1:3
            ES(index1+t,index2+t)=J2p/2;
        end
        
        nei=chart(i,j).right.neighbor2(4,:);
        index2=((nei(1)-1)*size+nei(2)-1)*6;
        for t=1:1:3
            ES(index1+t+3,index2+t+3)=J2p/2;
        end


        for k=2:1:3
            nei=chart(i,j).left.neighbor2(k,:);
            index2=((nei(1)-1)*size+nei(2)-1)*6;
            for t=1:1:3
                ES(index1+t,index2+t)=J2/2;
            end
            nei=chart(i,j).right.neighbor2(k,:);
            index2=((nei(1)-1)*size+nei(2)-1)*6;
            for t=1:1:3
                ES(index1+t+3,index2+t+3)=J2/2;
            end
        end
        
        for k=5:1:6
            nei=chart(i,j).left.neighbor2(k,:);
            index2=((nei(1)-1)*size+nei(2)-1)*6;
            for t=1:1:3
                ES(index1+t,index2+t)=J2/2;
            end
            nei=chart(i,j).right.neighbor2(k,:);
            index2=((nei(1)-1)*size+nei(2)-1)*6;
            for t=1:1:3
                ES(index1+t+3,index2+t+3)=J2/2;
            end
        end



        for k=1:1:3
            nei=chart(i,j).left.neighbor3(k,:);
            index2=((nei(1)-1)*size+nei(2)-1)*6;
            for t=1:1:3
                ES(index1+t,index2+t+3)=J3/2;
            end
            nei=chart(i,j).right.neighbor3(k,:);
            index2=((nei(1)-1)*size+nei(2)-1)*6;
            for t=1:1:3
                ES(index1+t+3,index2+t)=J3/2;
            end
        end

       

        nei=chart(i,j).left.neighbor(1,:);
        index2=((nei(1)-1)*size+nei(2)-1)*6;
        ES(index1+1:index1+3,index2+4:index2+6)=ES(index1+1:index1+3,index2+4:index2+6)+inv(U)*inv(order(1:3,:))*Jp*order(1:3,:)*U/2;
        nei=chart(i,j).right.neighbor(1,:);
        index2=((nei(1)-1)*size+nei(2)-1)*6;
        ES(index1+4:index1+6,index2+1:index2+3)=ES(index1+4:index1+6,index2+1:index2+3)+inv(U)*inv(order(1:3,:))*Jp*order(1:3,:)*U/2;
        nei=chart(i,j).left.neighbor(2,:);
        index2=((nei(1)-1)*size+nei(2)-1)*6;
        ES(index1+1:index1+3,index2+4:index2+6)=ES(index1+1:index1+3,index2+4:index2+6)+inv(U)*inv(order(4:6,:))*Jp*order(4:6,:)*U/2;
        nei=chart(i,j).right.neighbor(2,:);
        index2=((nei(1)-1)*size+nei(2)-1)*6;
        ES(index1+4:index1+6,index2+1:index2+3)=ES(index1+4:index1+6,index2+1:index2+3)+inv(U)*inv(order(4:6,:))*Jp*order(4:6,:)*U/2;
        nei=chart(i,j).left.neighbor(3,:);
        index2=((nei(1)-1)*size+nei(2)-1)*6;
        ES(index1+1:index1+3,index2+4:index2+6)=ES(index1+1:index1+3,index2+4:index2+6)+inv(U)*inv(order(7:9,:))*J*order(7:9,:)*U/2;
        nei=chart(i,j).right.neighbor(3,:);
        index2=((nei(1)-1)*size+nei(2)-1)*6;
        ES(index1+4:index1+6,index2+1:index2+3)=ES(index1+4:index1+6,index2+1:index2+3)+inv(U)*inv(order(7:9,:))*J*order(7:9,:)*U/2;


    end
end
end

