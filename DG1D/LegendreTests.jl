include("dg_utils.jl")

α = 0
β = 0

for n in [1,2,4]
    println("!!!    for n = $n")

    r = jacobiGL(α, β, n)
    println("r = ")
    show(stdout, "text/plain", r)
    println()

    V = ones(length(r),length(r))
    vandermonde!(V, r, α, β)
    println("V = ")
    show(stdout, "text/plain", V)
    println()

    dV = similar(V)
    dvandermonde!(dV, r, α, β)
    println("dV = ")
    show(stdout, "text/plain", dV)
    println()

    D = dmatrix(r, α, β)
    println("D = ")
    show(stdout, "text/plain", D)
    println()

end
