using SparseArrays
using BandedMatrices
using LinearAlgebra

size = 5000
band = 500
b∇² = BandedMatrix(randn(size,size), (band,band))
a = randn(size)
∇² = Array(b∇²)

s∇² = sparse(b∇²)
dropzeros!(s∇²)

println("banded")
@btime c = b∇² \ a;
println("full")
@btime c = ∇² \ a;
println("sparse")
@btime c = s∇² \ a;

#compute factorization
qrb = qr(b∇²)
lus = lu(s∇²)
luf = lu(∇²)

println("factored banded")
@btime c = qrb \ a;
println("factored sparse")
@btime c = lus \ a;
println("factored full")
@btime c = luf \ a;
