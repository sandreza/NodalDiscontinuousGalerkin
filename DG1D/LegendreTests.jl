include("../src/utils.jl")

α = 0
β = 0

for n in [1, 2, 4]
    println("!!!    for n = $n")

    r = jacobiGL(α, β, n)
    println("r = ")
    show(stdout, "text/plain", r)
    println()

    V = vandermonde(r, α, β, n)
    println("V = ")
    show(stdout, "text/plain", V)
    println()

    dV = dvandermonde(r, α, β, n)
    println("dV = ")
    show(stdout, "text/plain", dV)
    println()

    D = dmatrix(r, α, β, n)
    println("D = ")
    show(stdout, "text/plain", D)
    println()

end
