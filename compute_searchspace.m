function searchspace = compute_searchspace(Gstar, node, depth)
% Compute constrained search space for a given node from
% the edges E^star.
% Juri Kuronen (2018)
%
% Input:
% - Gstar: logical dxd matrix, graph G^star with edges E^star computed using plr.m
% - node: index of the node
% - depth: the search space will contain all nodes within graph-
%   theoretic-distance 'depth' from given node with respect to E^star
%
% Output:
% - searchspace: constrained search space for a given node.

    searchspace = Gstar(:, node);
    for i = 2:depth
        searchspace = (searchspace + sum(Gstar(:, find(searchspace)), 2)) > 0;
    end
    searchspace(node) = 0;
    searchspace = find(searchspace);

end
