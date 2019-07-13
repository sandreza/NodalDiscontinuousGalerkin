 include("../DG1D/run_advection.jl")
# make sure to run the not periodic case
i = 160

p1 = plot(x, sol.u[i], xlims=(0,L), ylims = (-1.1,1.1), marker = 3,    leg = false)
display(p1)
