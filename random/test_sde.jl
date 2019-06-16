using DifferentialEquations
const κ = 1
function lorenz(du,u,p,t)
  du[1] = -u[1]
  du[2] = -κ * u[2] + u[1] * u[3]
  du[3] = -κ * u[3] - u[1] * u[2] + u[1]
end

function σ_lorenz(du,u,p,t)
  du[1] = sqrt(2)
  du[2] = 0
  du[3] = 0
end
#http://docs.juliadiffeq.org/latest/solvers/sde_solve.html for different methods
prob_sde_lorenz = SDEProblem(lorenz,σ_lorenz,[0.0,1.0,0.0],(0.0,100.0))
#sol = solve(prob_sde_lorenz)
monte_prob = MonteCarloProblem(prob_sde_lorenz)
dt = 0.1
#SRA(), EM(), SOSRA
ensemble_members = 1000
sol = solve(monte_prob, SOSRA(), num_monte=ensemble_members, parallel_type=:threads, dt=dt, adaptive = false)

summ = MonteCarloSummary(sol,0:0.01:1)
plot(summ,labels="Middle 95%")
summ = MonteCarloSummary(sol,0:0.01:1;quantiles=[0.25,0.75])
plot!(summ,labels="Middle 50%",legend=true)

t = sol.u[1].t #since all the same with same adaptive timestepping
tlength = length(t)
mem_number = length(sol.u[1].u[1])
u = ones(ensemble_members, tlength, mem_number)
#sol.u[ensemble_member][time][member]
for k in 1:mem_number
  for j in 1:tlength
    for i in 1:ensemble_members
      u[i,j,k] = sol.u[i].u[j][k]
    end
  end
end

function flatten(u)
  reshape(u,length(u))
end
