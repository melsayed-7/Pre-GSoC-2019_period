using DifferentialEquations
using DiffEqOperators
using Plots; gr()


function SemiLinearHeat(;non_homo::Function,tspan,u0,B1::Symbol,B2::Symbol, params = [1.0,0.5])

    dx = (u0[end]-u0[1])/(length(u0)-1) # Deciding the step in the x variable

    N = length(u0)       # Deciding the length of the u0 to be able to build Our Derivative


    # Dirichlet Boundary Conditions

    if B1 in  Set([:Dirichlet0,:Dirichlet]) && B2 in Set([:Dirichlet0,:Dirichlet])
        A= DerivativeOperator{Float64}(2,2,dx,N,B1,B2;BC=(u0[1],u0[end]))
    end


    # Neumann Boundary Conditions
    if B1 in  Set([:Neumann,:Neumann0]) && B2 in  Set([:Neumann,:Neumann0])
        B = DerivativeOperator{Float64}(1,2,dx,N,:None,:None)
        deriv_start, deriv_end = (B*u0)[1], (B*u0)[end]
        A =  DerivativeOperator{Float64}(2,2,dx,N,B1,B2;BC=(deriv_start , deriv_end))
    end


    # Robin Boundary Condition
    if B1 in  Set([:Robin]) && B2 in  Set([:Robin])
        B = DerivativeOperator{Float64}(1,2,dx,N,:None,:None)
        deriv_start, deriv_end = (B*u0)[1], (B*u0)[end]
        Param = params
        left_RBC = params[1]*u0[1] - params[2]*deriv_start
        right_RBC = params[1]*u0[end] + params[2]*deriv_end
        A= DerivativeOperator{Float64}(2, 2, dx, N, B1, B2;
        BC = ((params[1],params[2],left_RBC),(params[1],params[2],right_RBC)))
    end

    f(u,p,t) = A*u + non_homo(u)  # This is the heat equation with a non-linear Part
    prob = ODEProblem(f , u0 ,tspan)
end



# A problem example

dx = 2*pi/15
x = collect(-pi : dx : pi)
u0 = @. -(x - 0.5).^2 + 1/12





non_homo(u) =3*u
prob = SemiLinearHeat(non_homo = non_homo ,tspan = (0.0,5.0),u0=u0,B1=:Neumann,B2=:Neumann)
plot(solve(prob))

    
