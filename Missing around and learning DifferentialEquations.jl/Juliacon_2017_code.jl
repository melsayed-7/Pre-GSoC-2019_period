using DifferentialEquations

f(u,t) = 0.98*u
u0 = 1.0
tspan = (0.0,5.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob)

using Plots; gr()
plot(sol)



plot(sol,linewidth=5,title="Solution to the linear ODE with a thick line",
     xaxis="Time (t)",yaxis="u(t) (in μm)",label="My Thick Line!")

plot!(sol.t,t->0.5*exp(1.01t),lw=3,ls=:dash,label="True Solution!")


###############
# Solving system of Equaitons

function lorenz(du,u,p,t)
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3] - u[2])
    du[3] = u[1] * u[2] - (8.0/3.0) * u[3]
end

u0 = [1.0;0.0;0.0]
tspan = (0.0,1.0)

prob = ODEProblem(lorenz,u0,tspan)
sol = solve(prob)
plot(sol , vars=[1,2,3])
# animate(sol,lw=3,every=4)
function lorenz(du,u,p,t)
 du[1] = 10.0*(u[2]-u[1])
 du[2] = u[1]*(28.0-u[3]) - u[2]
 du[3] = u[1]*u[2] - (8/3)*u[3]
end

u0 = [1.0;0.0;0.0]
tspan = (0.0,100.0)
prob = ODEProblem(lorenz,u0,tspan)
gr()
sol = solve(prob,vars=[1,2,3])

using StaticArrays

#######################
using DiffEqOperators

order = 2
deriv = 2

Δx =  0.1
N = 9
A = DerivativeOperator{Float64}(order,deriv,Δx,N,:Dirichlet0,:Dirichlet0)

full(DerivativeOperator{Float64}(order,deriv,Δx,N,:Dirichlet0,:Dirichlet0))
