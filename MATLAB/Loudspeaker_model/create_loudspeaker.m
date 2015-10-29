function [ loudspeaker ] = create_loudpseaker( Loudspeaker, N )

R_diaph    = Loudspeaker.Diaphragm_radius;
r_dustcap  = Loudspeaker.DustCap_radius;
r_vc       = Loudspeaker.VoiceCoil_radius;
h_diaph    = Loudspeaker.Diaphragm_height;
R_dustcap  = Loudspeaker.DustCap_sphere_radius;
w_susp     = Loudspeaker.SuspensionWidth;
%%
% Create diaphragm
% Create a circle
diaphragm_w_susp = create_circle(R_diaph + w_susp ,N);

r = sqrt(diaphragm_w_susp.Nodes(:,2).^2 + diaphragm_w_susp.Nodes(:,3).^2);
slope = h_diaph / ( R_diaph - r_vc ) ;
diaphragm_w_susp.Nodes(:,4) = (r - R_diaph)*slope;

[~, elements] = mesh_select(diaphragm_w_susp, sprintf('r < (%g+0.000)', r_vc), 'ID', 'any');
diaphragm_w_susp.Elements(elements,:) = [];
diaphragm_w_susp = drop_unused_nodes(diaphragm_w_susp);
diaphragm_w_susp.Elements = drop_IDs(diaphragm_w_susp);
diaphragm_w_susp.Nodes(1:end,1) = (1:length(diaphragm_w_susp.Nodes));

r = sqrt(diaphragm_w_susp.Nodes(:,2).^2 + diaphragm_w_susp.Nodes(:,3).^2);
[node_suspension, element_suspension] = mesh_select(diaphragm_w_susp, sprintf('r > (%g+0.000)', R_diaph), 'ID', 'any');

% Suspension: flat or circular
%diaphragm.Nodes(node_suspension,4) = 0;
diaphragm_w_susp.Nodes(node_suspension,4) = sqrt( (w_susp/2 + 1E-6)^2 - (r(node_suspension)-(R_diaph + w_susp/2) ).^2 );

diaphragm_w_susp.Elements(element_suspension,3:4) = 1*ones(size(element_suspension,1),2);

plot_mesh(diaphragm_w_susp);
%%
% Dust cap
eps = 5E-3;
dustcap = create_circle(R_diaph+ w_susp,N);
r = sqrt(dustcap.Nodes(:,2).^2 + dustcap.Nodes(:,3).^2);
dustcap.Nodes(:,4) = (r - R_diaph)*slope;
[~, elements] = mesh_select(dustcap, sprintf('r > (%g+0.000)', r_dustcap+eps), 'ID', 'any');
z0 = -slope*(R_diaph - r_dustcap) - sqrt(R_dustcap^2 - r_dustcap^2);
ind1 = find( r <= r_dustcap );
z_coord = real(sqrt(R_dustcap^2-r.^2))+z0;
dustcap.Nodes(ind1,4) = z_coord(ind1);
dustcap.Elements(elements,:) = []; dustcap = drop_unused_nodes(dustcap);
dustcap.Elements = drop_IDs(dustcap); 
dustcap.Nodes(1:end,1) = (1:length(dustcap.Nodes));


plot_mesh(dustcap);
%%
loudspeaker  = join_meshes(diaphragm_w_susp,dustcap);
loudspeaker = merge_coincident_nodes(loudspeaker,0.0001);
loudspeaker = drop_unused_nodes(loudspeaker);
loudspeaker.Elements = drop_IDs(loudspeaker);
loudspeaker.Nodes(1:end,1) = (1:length(loudspeaker.Nodes));

[node_suspension, element_suspension] = mesh_select(loudspeaker, sprintf('r > (%g+0.000)', R_diaph), 'ID', 'any');
end

