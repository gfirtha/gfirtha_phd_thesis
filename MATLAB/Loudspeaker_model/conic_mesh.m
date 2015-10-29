function [ diaphragm ] = conic_mesh( r_membrane, r_dustcap, r_vc, h_membrane, R_dustcap, N, w_surr)

w_surr = w_surr/1000;       % convert the surround width into meters
% h0 = h_membrane*( r_dustcap/R - 1 );
% x0 = sqrt(R_dustcap^2-r_dustcap^2);
% orig = -(x0 + abs(h0));
% 
% NiHu_mesh1 = create_circle(R,N);
% 
% r = sqrt(NiHu_mesh1.Nodes(:,2).^2+NiHu_mesh1.Nodes(:,3).^2);
% NiHu_mesh1.Nodes(:,4) = h_membrane*( r/R - 1 );
% 
% temp = find(r<=r_dustcap);
% y = real(sqrt(R_dustcap^2-r.^2))+orig;
% NiHu_mesh1.Nodes(temp,4) = y(temp);
% 
% NiHu_mesh1.Nodes(NiHu_mesh1.Nodes(:,4)> -h_membrane/50,4) = -w_susp;
% NiHu_mesh1.Nodes(:,4) = NiHu_mesh1.Nodes(:,4)+w_susp;
% 
% plot_mesh(NiHu_mesh1);

%%
% Create a circle
diaphragm = create_circle(r_membrane,N);

% Create diapraghm shape ()
r = sqrt(diaphragm.Nodes(:,2).^2+diaphragm.Nodes(:,3).^2);
diaphragm.Nodes(:,4) = h_membrane*(r/r_membrane-1);

% Delete nodes and elements which would be inside the voice coil 
%(choose elements, delete elements, delete nodes and rearrange nodes)
[node_voice_coil elements] = mesh_select(diaphragm, sprintf('r < (%g+0.000)', r_vc), 'ID', 'any');
diaphragm.Elements(elements,:) = [];
diaphragm = drop_unused_nodes(diaphragm);
diaphragm.Elements = drop_IDs(diaphragm);
diaphragm.Nodes(1:end,1) = (1:length(diaphragm.Nodes));

% Create suspension (surround)
diaphragm.Nodes(diaphragm.Nodes(:,4)> -w_surr,4) = -w_surr;
diaphragm.Nodes(:,4) = diaphragm.Nodes(:,4)+w_surr;

NiHu_mesh3 = join_meshes(diaphragm,diaphragm);
NiHu_mesh3.Nodes(NiHu_mesh3.Nodes(:,4)> -w_surr,4) = -w_surr;
NiHu_mesh3.Nodes(:,4) = NiHu_mesh3.Nodes(:,4)+w_surr;

% 
diaphragm = merge_coincident_nodes(NiHu_mesh3,0.00035);

% Again, drop unused/duplicated nodes and elements, and rearrange
diaphragm = drop_unused_nodes(diaphragm);
diaphragm.Elements = drop_IDs(diaphragm);
diaphragm.Nodes(1:end,1) = (1:length(diaphragm.Nodes));
end

