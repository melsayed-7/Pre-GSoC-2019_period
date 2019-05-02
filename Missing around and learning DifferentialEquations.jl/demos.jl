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




######################
# This is SemiLinearHeat but using FFTW

using FFTW
x=range(0,stop=2π,length=100)
u(x) = sin(x)
freqs = fft(u.(x))[ 1 : length(x) ÷ 2 + 1]  ## ÷ this means integer division
c= 2*abs.(freqs/length(x))


######################
# instead of the above section
using ApproxFun
s = Fourier()
n = 100
T = ApproxFun.plan_transform(s,n)
Ti = ApproxFun.plan_itransform(s,n)

x = points(s,n)    # 0  to 2*pi
r = (T*(cos.(cos.(x.-0.1))))


der = Derivative(s,2)







####################
# PDE-derived ODE


using Sundials
using Plots; gr()
using ApproxFun
using Sundials
using LinearAlgebra
using LinearMaps

import Base.At_mul_B!

S = Fourier()
n=100
x = points(S,n)
T = ApproxFun.plan_transform(S,n)
Ti = ApproxFun.plan_itransform(S,n)


u0= T*cos.(cos.(x.-0.1))

D2 = Derivative(S,2)
L = D2[1:n,1:n]
tspan = (0.0,10.0)


function heat(du,u,L,t)
    mul!(du,L,u)
end

prob = ODEProblem(heat, u0,tspan,L)

plot(x,Ti*sol(0.0))
plot!(x,Ti*sol(0.5))
plot!(x,Ti*sol(2.0))
plot!(x,Ti*sol(3.0))
plot!(x,Ti*sol(10.0))
