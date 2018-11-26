function mat2utrivec{T}(m::Matrix{T})
    assert(size(m,1) == size(m,2))
    n = size(m,1)
    [m[i,j] for i=1:n, j=1:n if i < j]
end

function mat2ltrivec{T}(m::Matrix{T})
    assert(size(m,1) == size(m,2))
    n = size(m,1)
    [m[j,i] for i=1:n, j=1:n if i < j]
end
