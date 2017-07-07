function ConvertToDot(L, fileName)
% function ConvertToDot(L, filename)
%
% L = (nx3) output from the linkage function
% fileName = output filename where the dot file output is written
%
[n,d] = size(L);

if (d ~= 3)
  fprintf(1, 'L (%d,%d) is not the output of a linkage function\n', n, d);
  return
else
end

if (0 == ischar(fileName))
  fprintf(1, 'fileName should be a string\n');
end

fprintf(1, 'Converting linkage of %d original points, output in %s\n', ... 
           n+1, fileName);

% 1. Create an array where all original nodes have weight zero
V = zeros(n+1, 1);
E = zeros(2*n, 2);

% 2. Go over each linkage and create node W 
newNodeIndex = (n+2);
edgeIndex = 1;
for i=1:n

  % 2.1 Each link row is two edges with an implicit parent
  leftChild = L(i, 1);
  rightChild = L(i, 2);
  linkWeight = L(i,3);

  % 2.2 Parent is always a new node
  parent = newNodeIndex;
  newNodeIndex = newNodeIndex + 1;

  % 2.3 Save unatlered link weight of parent
  V(parent) = linkWeight;

  % 2.4 Subtract the (max) height of the children from the current height
  %     and spread evenly amongst the children
  linkWeight = linkWeight - max(V(leftChild), V(rightChild));
  V(leftChild) = linkWeight/2;
  V(rightChild) = linkWeight/2;

  % 2.5 Save the edges for posterity
  E(edgeIndex,:) = [parent leftChild];
  E(edgeIndex+1,:) = [parent rightChild];
  edgeIndex = edgeIndex + 2;

end

% 5. The root is the last node, which has weight 0, technically
V(size(V,1)) = 0;

% 4. Print out everything in the requested file
of = fopen(fileName, 'w');
fprintf(of, 'digraph G {\n');
for i=1:size(V,1) fprintf(of, '  %d [weight=%f] ;\n', i-1, V(i)); end
for i=1:size(E,1) fprintf(of, '  %d -> %d ;\n', E(i,1)-1, E(i,2)-1); end
fprintf(of, '}');
fclose(of);

