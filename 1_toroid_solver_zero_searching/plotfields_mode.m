%$Id: plotfields_mode.m,v 1.5 2011/06/13 03:47:47 taolu Exp $
%$Revision: 1.5 $
%$Author: taolu $
%$Date: 2011/06/13 03:47:47 $
%function figid=plotfields_mode(mesh_X,mesh_Y,Field,Phi,figid,levels,file_label)
function figid=plotfields_mode(fem,Field,figid,file_label)
        
% field_str={'E_r';'E_z';'E_\theta';'Z_0H_r';'Z_0H_z';'Z_0H_\theta'};   
field_str={'E_r';'E_z';'E_\theta';'H_r';'H_z';'H_\theta'};
indx=[1 2 4 5 3 6];

% for count=1:length(indx)
%     figure(figid);
%     figid=figid+1;
%     clf
%     %  subplot(3,2,count);
%     geomplot(fem,'pointmode','off');
%     hold on;
%     Field(count,isnan(Field(count,:)))=0;
%     clim_min=min(real(Field(count,:)));
%     clim_max=max(real(Field(count,:)));
%     pdeplot(fem.mesh.p,[],fem.mesh.t,'xydata',real(Field(count,:)),'contour','off','colormap','jet','colorbar','on');
%     %     set(gca,'clim',[clim_min clim_max]);
%     title(field_str{count});
%     xlabel('X [\mum]');ylabel('Y [\mum]');
%     
% end %for count

% figure(figid);
% figid=figid+1;
% clf
% geomplot(fem,'pointmode','off');
% hold on;
% pdemesh(fem.mesh.p,fem.mesh.e,fem.mesh.t);
% title('mesh');


% if exist('file_label')
%     if ~isempty(file_label)
%         saveas(gcf,[file_label '_mesh.fig'], 'fig');
%     end
% end

I=abs((Field(1,:).*conj(Field(5,:)))-(Field(2,:).*conj(Field(4,:))));
pic=figure(figid);
clf;
figid=figid+1;
geomplot(fem,'pointmode','off');
hold on;
pdeplot(fem.mesh.p,[],fem.mesh.t,'xydata',I,'contour','off','colormap','jet','colorbar','on');
set(gca,'clim',[min(I(:)) max(I(:))]);
title('Intensity');
% if exist('file_label')
%     if ~isempty(file_label)
%         saveas(gcf,[file_label int2str(figid-1) '_I_ave.fig'], 'fig');
%     end
% end
% if exist('file_label')
%     if ~isempty(file_label)
%         saveas(gcf,['Intensity_' file_label int2str(figid-1) '.png'], 'png');
%     end
% end

% saveppt('test.ppt',['fig(' num2str(figid-1) ')'],pic);

