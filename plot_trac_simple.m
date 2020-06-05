function [] = plot_trac_simple(tri,trac,points)
% Simple plotting function to generate a coneplot of traction on a deformed
% surface configuration for diagnostics
%
%
%--- INPUTS ---
% tri    : Surface triangulization (e.g. from 'calculateNormals')
% trac   : traction vectors from funCalculateTractions
%           - trac{i} is the vector component in the i'th direction (x,y,z)
% points : Points where the tractions are sampled
%
%--- OUTPUTS ---
% *None*
%
% NOTES
% ----------------------------------------------------------------------
% AKL, MP, 2020-04-17

%need faces and vertices from the DT
% fv.vertices = DT{1}{t}.Points;
% fv.faces = DT{1}{t}.ConnectivityList;
ti(:,1) = trac{1};
ti(:,2) = trac{2};
ti(:,3) = trac{3};
x(:,1) = points{1};
x(:,2) = points{2};
x(:,3) = points{3};

tmag = sqrt(sum(ti.^2,2));
% fv.vertices = fv.vertices(:,[2,1,3]); % convert to xyz format

%%
figure;
hold all
hc = coneplot(x(:,1),x(:,2),x(:,3),ti(:,1),ti(:,2),ti(:,3),0.03,'nointerp');
fvc = repmat(tmag(:)',[42 1]);
set(hc, 'facecolor', 'flat', 'facevertexcdata', fvc(:))
hc.EdgeColor = 'none';
hc.AmbientStrength = 0.6;
hc.DiffuseStrength = 0.75;
hc.SpecularStrength = 0.4;
hc.SpecularExponent = 3;
colormap(turbo)
% caxis([0 max(cmax)])
hCBc = colorbar('location','west','FontSize',14,'Color','white');
% freezeColors
% cbfreeze(hCBc)

% axis image;
hl = light;
lightangle(hl,160, 20)
view([132.61 16.5]);
camva('manual')
lighting gouraud

h = gcf;
set(gca,'colormap',gray(255),'Visible','off')
caxis([min(x(:,3)),max(x(:,3))])
hCBt = colorbar('location','east','FontSize',14,'Color','white');
set(h,'color','k');

hCBc.Label.String = 'Displacement Mag';
hCBt.Label.String = 'Traction Mag';
hCBc.Label.FontSize = 18;
hCBt.Label.FontSize = 18;
hCBc.Label.Color = 'white';
hCBt.Label.Color = 'white';


% colorbar(hCBt,'hide')
% colorbar(hc)

% % hs = patch(fv);
% % hs.FaceColor = [60 60 60]/255;
% % hs.EdgeColor = 'none';
% % hs.AmbientStrength = 0.40;
% % hs.DiffuseStrength = 0.50;
% % hs.SpecularStrength = 0.4;
% % hs.SpecularExponent = 3;
% % axis image;
% % hl = light;
% % lightangle(hl,160, 20)
% % view([132.61 6.5]);
% % camva('manual')
% % lighting gouraud
% %
% % h = gcf;
% % set(gca,'color',0.9*[1,1,1],'Visible','off')
% % set(h,'color','k');
% %

end
