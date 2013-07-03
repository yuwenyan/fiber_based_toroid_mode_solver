
function L_mat=calc_L_mat(fiber_geom,fiber_modes,fem,lambda_guess,m_guess,scale_fac)

    x_mat=[fem.mesh.p(1,:);fem.mesh.p(1,:)]+fiber_geom.toroid_R;
    comp_indx=[4 5];
    num_eig_modes=fiber_geom.mode_count;
    L_mat=zeros(num_eig_modes,num_eig_modes);
    for count1=1:num_eig_modes
        for count2=1:num_eig_modes
            int_term_1=sum(area_integral(fem,conj(fiber_modes.Field{count1}(comp_indx,:)).*(-((m_guess./x_mat).^2)).*fiber_modes.Field{count2}(comp_indx,:)));
            int_term_2=sum(area_integral(fem,conj(fiber_modes.Field{count1}(comp_indx,:)).*(1./x_mat).*fiber_modes.Fieldx{count2}(comp_indx,:)));
            int_term_3=sum(area_integral(fem,conj(fiber_modes.Field{count1}(comp_indx,:)).*((fiber_modes.Neff(count2)*2*pi/lambda_guess)^2).*fiber_modes.Field{count2}(comp_indx,:)));
            L_mat(count1,count2)=L_mat(count1,count2)+int_term_1+int_term_2+int_term_3;
        end %for count2
    end %for count1
    L_mat=L_mat./scale_fac;
	