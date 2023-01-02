using Pkg

Pkg.activate(".")
Pkg.add("Plots")
Pkg.add("DifferentialEquations")
Pkg.update()

using Plots

# Setting up parameters (3 params)
k = 1.7 # rate of maximum population growth
M = 10000 # Carrying Capacity
r = 49 # (M/(y_0 - 1)) in the formula

# Packaging the parameters into one object, p::Vector..
p = [k, M, r] # set k = p[1], M = p[2], and r = p[3]

# Defining function according to the problem
f(P,p,t) = p[2]/(1 + (p[3] * exp( -p[1] * t )))

# Setting t = 0.0 since f(P,p,t) is independent of time
PRange = range(1, 10; length = 10)
dP_dt = [ f(P,p,0.0) for P in PRange ]

plt = plot( PRange, dP_dt # plot with basic feature
    ,linecolor=:black, marker=(:circle,:black)
    ,label="k=1.7"
)

# Enhancing resulting plot
plot!(plt
    ,title="Spread of Disease"
    ,xlabel="Time, t"
    ,ylabel="Population, P"
)

struct MyProblem
    f::Function
    u0
    tspan
    p
end

P0 = 200
trange = (0.0, 10.0)
prob0 = MyProblem(f, P0, trange, p)

typeof(prob0)

@show prob0.p
@show prob0.u0;
@show prob0.tspan
@show prob0.f( 80.0, prob0.p, 0.0 );

struct MySolution
    t
    u
end

t = range(0.0,10.0; length=101)
soln0 = MySolution( t, sin.(t) )

soln0

plt = plot( soln0.t, soln0.u
    ,lc=:black
    ,marker=(:circle,:black)
    ,label="k"
    ,xlabel="Time, t"
    ,ylabel="Population, P"
    ,title="Just a sample(Euler Method)"
)

function solvePls( prob::MyProblem; npoints=100 )
    # Ensuring one-to-one correspondence..
    u = [ prob.u0 ] # first entry of solution
    t = [ prob.tspan[1] ] # first entry of time param
    dt = (prob.tspan[2] - prob.tspan[1])/(npoints+1) #ensures clean end
    
    # Initializing..
    uold = prob.u0
    told = prob.tspan[1]
    
    # Looping steps 2 and 3..
    for _ in 1:npoints
        # Computing Euler step..
        unew = uold + dt*prob.f( uold, prob.p, told )
        tnew = told + dt
        
        # Collecting values..
        u = vcat(u,unew)
        t = vcat(t,tnew)
        
        # Exchanging values..
        uold = unew
        told = tnew
    end
    
    return MySolution(t,u)
end

soln0=solvePls(prob0; npoints=40)

plt = plot( soln0.t, soln0.u
    ,lc=:black
    ,marker=(:circle,:black)
    ,label="k=1.7"
    ,xlabel="Time, t"
    ,ylabel="Population, P"
    ,title="Spread of Disease"
)

using DifferentialEquations

P0 = 200
tspan = (0.0, 10.0)
prob = ODEProblem(f, P0, tspan, p)

soln = solve(prob)

plt = plot(soln.t[2:end], soln.u[2:end]
    ,lc=:black,marker=(:circle,:black)
    ,label="numerical"
)

plot!(plt
    ,title="ODE Solver results"
    ,xlabel="Time, t"
    ,ylabel="Population, P"
)

plt = plot(soln
    ,lc=:black
    ,label="ODE Solver (interpolated)"
)

plot!(plt, soln.t[2:end], soln.u[2:end]
    ,lc=:black,marker=(:circle,:black)
    ,label="(numerical data)"
)

plot!(plt
    ,title="ODE Solver results"
    ,xlabel="Time, t"
    ,ylabel="Temperature in °C"
)

plt = plot(soln
    ,lc=:black
    ,label="ODE Solver (interpolated)"
)

scatter!(plt, soln0.t, soln0.u
    ,markershape=:circle,mc=:black
    ,label="Euler method"
)

plot!(plt
    ,xlabel="Time (s)"
    ,ylabel="Temperature (°C)"
#    ,yaxis=:log
)

Pkg.add("Profile")
Pkg.add("ProfileView")
Pkg.update()

using Profile
using ProfileView

@profile soln0 = solvePls(prob0; npoints=40)

ProfileView.view()


