using DifferentialEquations
using DiffEqOperators
using Plots; gr()


function SemiLinearHeat(;non_homo::Function,tspan,u0,B1::Symbol,B2::Symbol, params = [1.0,0.5])

    N = length(u0)       # Deciding the length of the u0 to be able to build Our Derivative
    dx = 1/N # Deciding the step in the x variable

    # Dirichlet Boundary Conditions

    if B1 in  Set([:Dirichlet0,:Dirichlet]) && B2 in Set([:Dirichlet0,:Dirichlet])
        A= DerivativeOperator{Float64}(2,2,dx,N,:None,:None)
        function heat(du,u,L,t)
            mul!(du,L,u)
        end
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

    f(u,p,t) = A*u .+ non_homo(u)  # This is the heat equation with a non-linear Part
    prob = ODEProblem(heat , u0 ,tspan,A)


end







#######################################################
# More generalized function





u0= T*cos.(cos.(x.-0.1))

D2 = Derivative(S,2)
L = D2[1:n,1:n]
tspan = (0.0,10.0)


function heat(du,u,L,t)
    mul!(du,L,u)
end


function SemiLinearHeat(domain, discretization, u0, non_linear)
    if discretization == "ApproxFun"

        S = Fourier()
        n=100
        x = points(S,n)
        T = ApproxFun.plan_transform(S,n)
        Ti = ApproxFun.plan_itransform(S,n)
        u0= T*u0
        D2 = Derivative(S,2)
        L = D2[1:n,1:n]

        function heat(du,u,L,t)
            mul!(du,L,u)
        end

        prob = ODEProblem(heat, u0,domain,L)

    end

end



u0 = cos.(cos.(x.-0.1))
non_homo(u) =0;

prob = SemiLinearHeat((0.0,10.0), "ApproxFun", u0, non_homo)
sol = solve(prob)
plot(x,Ti*sol(10.0))
plot!(x,Ti*sol(0.5))
plot!(x,Ti*sol(2.0))
plot!(x,Ti*sol(3.0))

    
