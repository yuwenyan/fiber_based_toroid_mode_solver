%$Id: calc_T_matrix.m,v 1.9 2006/10/04 07:42:42 taolu Exp $
%$Revision: 1.9 $
%$Author: taolu $
%$Date: 2006/10/04 07:42:42 $
%function T_mat=calc_T_matrix(fiber_geom,fiber_modes,x_mat)
function T_mat=calc_T_matrix_2(fiber_geom,fiber_modes,fem)
if (1==0)
    % Create node points
    mesh_number=300;
    mesh_number_x=mesh_number; % number of x elements in the mesh
    mesh_number_y=mesh_number; % number of y elements in the mesh
    wnd_sz_X=3e-5;
    wnd_sz_Y=3e-5;
    x_vec=linspace(-wnd_sz_X/2.1,wnd_sz_X/2.1,mesh_number_x);
    y_vec=linspace(-wnd_sz_Y/2.1,wnd_sz_Y/2.1,mesh_number_y);
    mesh_spacing_x=x_vec(2)-x_vec(1);
    mesh_spacing_y=y_vec(2)-y_vec(1);
    [x,y]=meshgrid(x_vec,y_vec);
    node_points = [x(:)'; y(:)'];
    num_spmodes=fiber_modes.num_modes;
    fiber_comp_str={'Ex','Ey','Ez','Hx','Hy','Hz'};
    num_components=length(fiber_comp_str);
    fiber_fields=cell(num_spmodes,num_components);
    mode_id=fiber_modes.Neff_indx;
    M_Guess=2*pi*fiber_geom.toroid_R*fiber_geom.n_guide/fiber_geom.lambda;
    for count=1:num_spmodes
        for field_id=1:num_components
            fiber_fields{count,field_id}=postinterp(fem,fiber_comp_str{field_id},node_points,'Solnum',mode_id(count));
            fiber_fields{count,field_id}=reshape(fiber_fields{count,field_id},mesh_number_y,mesh_number_x);
        end
    end
    for count=1:num_spmodes
        for count_i=1:mesh_number
            for count_j=1:mesh_number
                if isnan(fiber_fields{count,6}(count_i,count_j))
                    fiber_fields{count,6}(count_i,count_j)=0;
                end
            end
        end
    end
    partial_Hz=zeros(mesh_number_y,mesh_number_x,num_spmodes);
    for count_i=1:num_spmodes
        for count=2:mesh_number_x-1
            partial_Hz(:,count,count_i)=(fiber_fields{count_i,6}(:,count+1)-fiber_fields{count_i,6}(:,count-1))./(2*mesh_spacing_x);
        end
    end
    matrix_A=zeros(num_spmodes,num_spmodes);
    %     matrix_B=zeros(num_spmodes,num_spmodes);
    for count_j=1:num_spmodes
        for count_i=1:num_spmodes
            matrix_A(count_j,count_i)=sum(sum(conj(fiber_fields{count_j,6}(:,:))./(x+fiber_geom.toroid_R).*(partial_Hz(:,:,count_i))))*mesh_spacing_x*mesh_spacing_y...
                -sum(sum(conj(fiber_fields{count_j,6}(:,:)).*fiber_fields{count_i,6}(:,:).*(M_Guess^2)./((x+fiber_geom.toroid_R).^2)))*mesh_spacing_x*mesh_spacing_y;
            %             if count_j==count_i
            %                 matrix_B(count_j,count_i)=(2*pi*n_eff(mode_id(count_i)))^2*sum(sum(conj(fiber_fields{count_j,6}(:,:)).*fiber_fields{count_i,6}(:,:)))*mesh_spacing_x*mesh_spacing_y;
            %             end
        end
    end
    %     T_mat=matrix_A;
    T_mat=-matrix_A./fiber_geom.toroid_R;
else
    x_mat=[fem.mesh.p(1,:);fem.mesh.p(1,:)]/fiber_geom.toroid_R;
    comp_indx=[4 5];
    num_eig_modes=fiber_geom.mode_count;
    T_mat=zeros(num_eig_modes,num_eig_modes);
    Neff=zeros(num_eig_modes,1);
    for count1=1:num_eig_modes
        for count2=1:num_eig_modes
            %                 T_mat(count1,count2)=T_mat(count1,count2)+sum(sum(conj(fiber_modes.Field{count1}(:,:,field_count)).*x_mat.*fiber_modes.Field{count2}(:,:,field_count)));
            int_factor=sum(area_integral(fem,conj(fiber_modes.Field{count1}(comp_indx,:)).*x_mat.*fiber_modes.Field{count2}(comp_indx,:)));
            T_mat(count1,count2)=T_mat(count1,count2)+int_factor;
        end %for count2
        Neff(count1)=fiber_modes.Neff(count1);
    end %for count1
    T_mat=T_mat*diag(Neff);
    T_mat=diag(Neff)+T_mat;
end
