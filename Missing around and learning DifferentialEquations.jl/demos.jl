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
