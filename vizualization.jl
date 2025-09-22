using Plots
using LinearAlgebra
using ColorVectorSpace

default(background_color = :white)

# Params
R = 3.0
levels = 5                 # includes poles
n_long = 9
dθ = 180.0/(levels-1)
dφ = 360.0/n_long

az_deg = 131.0
el_deg = 60.0

# Keep grid/points on (or slightly inside) the surface so they don't poke out
R_grid   = R * 0.999
R_points = R * 0.999

to_cart(θdeg, φdeg, r) = (
    r * sin(deg2rad(θdeg)) * cos(deg2rad(φdeg)),
    r * sin(deg2rad(θdeg)) * sin(deg2rad(φdeg)),
    r * cos(deg2rad(θdeg)),
)

# Great-circle interpolation between two Cartesian points on the sphere
function slerp(a::NTuple{3,Float64}, b::NTuple{3,Float64}, t::Float64)
    ra = sqrt(a[1]^2 + a[2]^2 + a[3]^2)
    rb = sqrt(b[1]^2 + b[2]^2 + b[3]^2)
    r = (ra + rb)/2              # both radii ~ R_grid
    ua = (a[1]/ra, a[2]/ra, a[3]/ra)
    ub = (b[1]/rb, b[2]/rb, b[3]/rb)
    ω = acos(clamp(ua[1]*ub[1] + ua[2]*ub[2] + ua[3]*ub[3], -1.0, 1.0))
    if isapprox(ω, 0.0; atol=1e-10)
        return (r*ua[1], r*ua[2], r*ua[3])
    end
    s0 = sin((1-t)*ω)/sin(ω)
    s1 = sin(t*ω)/sin(ω)
    p = (s0*ua[1] + s1*ub[1], s0*ua[2] + s1*ub[2], s0*ua[3] + s1*ub[3])
    return (r*p[1], r*p[2], r*p[3])
end

p = plot(legend=false, grid=false, framestyle=:none, size=(1200, 900),
         camera=(az_deg, el_deg), xlims=(-R,R), ylims=(-R,R), zlims=(-R,R))

# --- Sphere surface: plot the full sphere, let z-buffer hide the back ---
θs = range(0, stop=π, length=160)
φs = range(0, stop=2π, length=320)
X = [R*sin(θ)*cos(φ) for φ in φs, θ in θs]
Y = [R*sin(θ)*sin(φ) for φ in φs, θ in θs]
Z = [R*cos(θ)        for φ in φs, θ in θs]

# simple view-dependent light for nicer look (optional)
v̂ = let az = deg2rad(az_deg), el = deg2rad(el_deg)
    # unit vector from camera toward origin (approx)
    (-cos(el)*cos(az), -cos(el)*sin(az), -sin(el))
end
shade = [max(0, (X[i,j]/R, Y[i,j]/R, Z[i,j]/R)⋅v̂) for i in eachindex(φs), j in eachindex(θs)]
surface!(p, X, Y, Z, zcolor=shade, c=:greys, alpha=0.12, linealpha=0.0, cbar=false)

# --- Points ---
scatter!(p, [0.0],[0.0],[ R_points], mc=:black, ms=3.5)
scatter!(p, [0.0],[0.0],[-R_points], mc=:black, ms=3.5)
for r in 1:(levels-2)
    θ = r*dθ
    xs=Float64[]; ys=Float64[]; zs=Float64[]
    for ℓ in 0:(n_long-1)
        φ = ℓ*dφ
        x,y,z = to_cart(θ, φ, R_points)
        push!(xs,x); push!(ys,y); push!(zs,z)
    end
    scatter!(p, xs, ys, zs, mc=:black, ms=3.5)
end

# --- Longitudes / Latitudes (on-surface; back side will be hidden by the sphere) ---
for ℓ in 0:(n_long-1)
    φ = ℓ*dφ
    thetas = range(0.0, stop=180.0, length=180)
    xs=Float64[]; ys=Float64[]; zs=Float64[]
    for θ in thetas
        x,y,z = to_cart(θ, φ, R_grid)
        push!(xs,x); push!(ys,y); push!(zs,z)
    end
    plot!(p, xs, ys, zs, lc=:black, lw=0.5, la=0.6)
end

for r in 1:(levels-2)
    θ = r*dθ
    phis = range(0.0, stop=360.0, length=361)
    xs=Float64[]; ys=Float64[]; zs=Float64[]
    for φ in phis
        x,y,z = to_cart(θ, φ, R_grid)
        push!(xs,x); push!(ys,y); push!(zs,z)
    end
    plot!(p, xs, ys, zs, lc=:black, lw=0.5, la=0.6)
end

# --- Great-circle diagonals (slerp), e.g. for the two belts r=1 and r=2 ---
for r in (1,2)
    θA = r*dθ
    θB = (r+1)*dθ
    for ℓ in 0:(n_long-1)
        φA = (ℓ+1)*dφ
        φB =  ℓ    *dφ
        A = to_cart(θA, φA, R_grid)
        B = to_cart(θB, φB, R_grid)
        ts = range(0.0, stop=1.0, length=100)
        xs=Float64[]; ys=Float64[]; zs=Float64[]
        for t in ts
            x,y,z = slerp(A,B,t)
            push!(xs,x); push!(ys,y); push!(zs,z)
        end
        plot!(p, xs, ys, zs, lc=:black, lw=0.5, la=0.6)
    end
end

savefig(p, "sphere_grid.png")

