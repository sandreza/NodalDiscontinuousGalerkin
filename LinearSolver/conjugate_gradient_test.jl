include("conjugate_gradient.jl")


# simple tests (from wikipedia)
A = [4.0 1.0; 1.0 3.0]
b = [1.0; 2.0]
x⁰ = [2.0; 1.0]
solution = A \b

A_tmp(x) = A*x
B = inv(A)
pre_tmp(x) = B*x

x⁰ = [2.0; 1.0]
conjugate_gradient!(A_tmp, x⁰, b, maximum_iterations = 1)
println("the relative error after one iteration is ")
println(norm(x⁰ - solution) / norm(solution))
conjugate_gradient!(A_tmp, x⁰, b, maximum_iterations = 2)
println("the relative error after two iterations  is ")
println(norm(x⁰ - solution) / norm(solution))
x⁰ = [2.0; 1.0]
conjugate_gradient!(A_tmp, x⁰, b, maximum_iterations = 1, P = pre_tmp)
println("the relative error after one iteration with a perfect preconditioner is ")
println(norm(x⁰ - solution) / norm(solution))

###
# More complex text using 1D DG stuff
# include("../DG1D/runPoisson.jl")

solution = copy(sol[:])
b = copy(f[:])
∇²_tmp(x) = s∇² * x
G = inv(∇²) # Greens function
G_tmp(x) = G * x
###
# Laplacian
norm(∇²_tmp(solution) - b)/norm(b)
x⁰ = randn(length(sol))
r = conjugate_gradient!(∇²_tmp, x⁰, b, track_residual = true)
println("The relative error is")
println(norm(x⁰-solution)/norm(solution))
scatter(log.(r)/log(10), ylabel = "log10 residual norm", xlabel = "iterations", title = "Convergence of Conjugate Gradient with the Laplacian")
###
# Greens function
x⁰ = randn(length(x))
Gb = copy(solution)
r = conjugate_gradient!(G_tmp, x⁰, Gb,  track_residual = true)
println("The relative error is")
println(norm(x⁰-b)/norm(b))
scatter(log.(r)/log(10), ylabel = "log10 residual norm", xlabel = "iterations", title = "Convergence of Conjugate Gradient with the Green's function")

###
# Preconditioned Laplacian, Perfect
x⁰ = randn(length(sol))
P_tmp(x) = G * x
r = conjugate_gradient!(∇²_tmp, x⁰, b, track_residual = true, P = P_tmp)
println("The relative error is")
println(norm(x⁰-solution)/norm(solution))
scatter(log.(r)/log(10), ylabel = "log10 residual norm", xlabel = "iterations", title = "Laplacian with Perfect Preconditioner ")

###
# Preconditioned Laplacian, Inverse Diagonal
x⁰ = randn(length(sol))
prec = 1.0 ./ diag(s∇²)
P_tmp(x) = prec .* x
r = conjugate_gradient!(∇²_tmp, x⁰, b, track_residual = true, P = P_tmp)
println("The relative error is")
println(norm(x⁰-solution)/norm(solution))
scatter(log.(r)/log(10), ylabel = "log10 residual norm", xlabel = "iterations", title = "Laplacian with Bad Preconditioner ")

###
# Preconditioned Laplacian, Discrete Tridiagonal inverse Laplacian (-1 2 -1)
x⁰ = randn(length(sol))
prec = zeros(length(x⁰), length(x⁰))
for i in eachindex(x⁰)
    prec[i,i] = -2
    if i < length(x⁰)
        prec[i+1,i] = 1
        prec[i,i+1] = 1
    end
end
# this is a tridiagonal matrix thus doing this is a bad idea
prec = inv(prec)
P_tmp(x) = prec * x
r = conjugate_gradient!(∇²_tmp, x⁰, b, track_residual = true, P = P_tmp)
println("The relative error is")
println(norm(x⁰-solution)/norm(solution))
scatter(log.(r)/log(10), ylabel = "log10 residual norm", xlabel = "iterations", title = "Laplacian with Bad Preconditioner ")
