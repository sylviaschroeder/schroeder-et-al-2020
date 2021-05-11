function matrix = makeLambdaMatrix(krnlDims, smoothPerDim)

allDims = [krnlDims, prod(krnlDims)];
% create kernel of given size
matrix = zeros(allDims);
% create matrix of kernel size with 1 at each 
tmp = reshape(eye(prod(krnlDims)), [krnlDims, prod(krnlDims)]);
for dim = 1:length(krnlDims)
    if smoothPerDim(dim) == 0
        continue
    end
    prm = [dim, setdiff(1:length(allDims),dim)];
    m = zeros(krnlDims(dim), prod(allDims(prm(2:end))));
    t = reshape(permute(tmp, prm), krnlDims(dim), []);
    m(1:end-1,:) = m(1:end-1,:) - t(2:end,:) .* smoothPerDim(dim);
    m(2:end,:) = m(2:end,:) - t(1:end-1,:) .* smoothPerDim(dim);
    m = reshape(m, [krnlDims(prm(1:end-1)), prod(krnlDims)]);
    [~,backprm] = sort(prm);
    m = permute(m, backprm);
    matrix = matrix + m;
end
% reshape matrix to have one column per center pixel
matrix = reshape(matrix, [], prod(krnlDims));
% normalise so that all smoothing values per parameter sum to -1 (except on edges)
matrix = matrix ./ abs(median(sum(matrix,1)));
% set center value for each parameter to one
matrix = matrix + reshape(tmp, [], allDims(end));