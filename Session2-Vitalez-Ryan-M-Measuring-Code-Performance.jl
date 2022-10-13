using Pkg;
Pkg.add("BenchmarkTools");
Pkg.add("Plots");
Pkg.add("Profile");
Pkg.add("ProfileView");
Pkg.status()

# I will be using a kinematic equation: t = √(2d/a) where vᵢ = 0.
function myroot(x::Float64; nloop=5)
    t = 0.5x
    for _ in 1:nloop
        t = sqrt(2d/a)
    end
    return t
end

d = 4.0

@show myroot(a) - sqrt(a)

@show eps()

@time myroot(500.0, nloop=10)

using BenchmarkTools

mark1 = @benchmark for _ in 1:10_000_000 myroot(500.0, nloop=10) end

mark0 = @benchmark for _ in 1:10_000_000 sqrt(500.0) end

speedx = median(mark0.times)/median(mark1.times)
bestx = minimum(mark0.times)/median(mark1.times)

@show speedx;
@show bestx;

using Profile

@profile myroot

Profile.print()
