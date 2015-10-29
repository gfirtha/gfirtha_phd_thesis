function [ node_sequence ] = arrange_it( node_per,NiHu_mesh )

x = NiHu_mesh.Nodes(node_per,2);
y = NiHu_mesh.Nodes(node_per,3);
unordered = node_per;

coord = [x y];

orig = coord(1,:);
unordered(1,:) = [];
coord(1,:) = [];

node_sequence = zeros(length(node_per),1);
node_sequence(1,1) = node_per(1);

for n = 1:length(node_per)-1
    temp = coord - repmat(orig,length(coord(:,1)),1);
    [val, ind] = min(sqrt(temp(:,1).^2+temp(:,2).^2));
    node_sequence(n+1,1) = unordered(ind);
    unordered(ind,:) = [];
    coord(ind,:) = [];

end
end

