function [ NiHu_mesh, v_z ] = post_process_mesh( NiHu_mesh, v, DOF, zind, Loudspeaker )

R = Loudspeaker.Diaphragm_radius;
h = Loudspeaker.Diaphragm_height;
r_vc = Loudspeaker.VoiceCoil_radius;
r_dc = Loudspeaker.DustCap_radius;
w = Loudspeaker.MembraneThickness;

orig = NiHu_mesh.Nodes(:,4);
node_ind = floor(DOF(zind));
v_z = v(zind);

slope = h / ( R - r_vc ) ;
z0 = slope*(R-r_dc);
[node_exc,elements] = mesh_select(NiHu_mesh, sprintf('z < (%g+0.0001)', -z0-1E-3), 'ID', 'any');
NiHu_mesh.Nodes(node_exc,:) = [];
NiHu_mesh.Elements(elements,:) = [];
NiHu_mesh = drop_unused_nodes(NiHu_mesh);
NiHu_mesh.Elements = drop_IDs(NiHu_mesh);
NiHu_mesh.Nodes(1:end,1) = (1:length(NiHu_mesh.Nodes));
v_z(node_exc)= [];

end