using DifferentialEquations

f(u,t) = 0.98*u
u0 = 1.0
tspan = (0.0,5.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob)

using Plots
plot(sol)

gui()

plot(sol,linewidth=5,title="Solution to the linear ODE with a thick line",
     xaxis="Time (t)",yaxis="u(t) (in μm)",label="My Thick Line!")

plot!(sol.t,t->0.5*exp(1.01t),lw=3,ls=:dash,label="True Solution!")
