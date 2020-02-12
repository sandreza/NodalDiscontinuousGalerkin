include("conjugate_gradient.jl")


# simple tests
A = [4.0 1.0; 1.0 3.0]
b = [1.0; 2.0]
x⁰ = [2.0; 1.0]
solution = A \b

tmp(x) = A*x
B = inv(A)
pre_tmp(x) = B*x


x⁰ = [2.0; 1.0]
conjugate_gradient!(tmp, x⁰, b, maximum_iterations = 1)
println("the relative error after one iteration is ")
println(norm(x⁰ - solution) / norm(solution))
conjugate_gradient!(tmp, x⁰, b, maximum_iterations = 2)
println("the relative error after two iterations  is ")
println(norm(x⁰ - solution) / norm(solution))
x⁰ = [2.0; 1.0]
conjugate_gradient!(tmp, x⁰, b, maximum_iterations = 1, P = pre_tmp)
println("the relative error after one iteration with a perfect preconditioner is ")
println(norm(x⁰ - solution) / norm(solution))

###
# More complex text using 1D DG stuff
