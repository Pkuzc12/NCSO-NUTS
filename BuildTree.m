function [var_1,var_2,var_3,var_4,var_5,var_6,var_7,var_8,var_9] = BuildTree(x,p,u,v,j,epsilon,x0,p0,kT,size,ES,parameter)

if j==1

    [x_temp,p_temp]=Leapfrog(x,p,v*epsilon,size,ES,parameter);

    n_temp=(u<=Boltzmann(x_temp,kT,size,ES,parameter));
    s_temp=(u<BoltzmannLimit(x_temp,kT,size,ES,parameter,-1.7639));

    [var_1, var_2, var_3, var_4, var_5, var_6, var_7, var_8, var_9] = deal(x_temp,p_temp,x_temp,p_temp,x_temp,n_temp,s_temp,min([1,Boltzmann(x_temp,kT,size,ES,parameter)/Boltzmann(x0,kT,size,ES,parameter)]),1);

else

    [x_neg,p_neg,x_pos,p_pos,x_temp,n_temp,s_temp,alpha_temp,n_alpha_temp]=BuildTree(x,p,u,v,j-1,epsilon,x0,p0,kT,size,ES,parameter);

    if s_temp==1

        if v==-1

            [x_neg,p_neg,~,~,x_temp_temp,n_temp_temp,s_temp_temp,alpha_temp_temp,n_alpha_temp_temp]=BuildTree(x_neg,p_neg,u,v,j-1,epsilon,x0,p0,kT,size,ES,parameter);
        
        else

            [~,~,x_pos,p_pos,x_temp_temp,n_temp_temp,s_temp_temp,alpha_temp_temp,n_alpha_temp_temp]=BuildTree(x_pos,p_pos,u,v,j-1,epsilon,x0,p0,kT,size,ES,parameter);

        end

        if rand(1)<(n_temp_temp)/(n_temp+n_temp_temp)

            x_temp=x_temp_temp;

        end

        alpha_temp=alpha_temp+alpha_temp_temp;
        n_alpha_temp=n_alpha_temp+n_alpha_temp_temp;

        s_temp=s_temp_temp*(dot((x_pos-x_neg),p_neg)>=0)*(dot((x_pos-x_neg),p_pos)>=0);
        n_temp=n_temp+n_temp_temp;

    end

    [var_1,var_2,var_3,var_4,var_5,var_6,var_7,var_8,var_9] = deal(x_neg,p_neg,x_pos,p_pos,x_temp,n_temp,s_temp,alpha_temp,n_alpha_temp);


end

end
