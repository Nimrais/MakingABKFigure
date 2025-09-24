using Plots
using LinearAlgebra
using ColorVectorSpace
using ColorSchemes
default(background_color = :white)

# Params
R = 1.0
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







# --- Create 2D projection plot ---
# For each point on the sphere, compute (arcsin(z/R), atan(y,x))



# Define the domain
ϕ = range(-π/2, π/2, length=300)
θ = range(-π, π, length=300)

# Create a 2D array of cos(ϕ) values (constant along θ)
Z = [cos(phi) for phi in ϕ, _ in θ]

# Plot the heatmap with reversed grayscale
p2 = heatmap(θ, ϕ, Z,
    color = cgrad(:grays, rev=true),
    xlabel = "θ",
    ylabel = "φ",
    title = "Shading of the density",
    colorbar_title = "cos(φ)", 
    xticks = (-π:π/2:π, ["-π", "-π/2", "0", "π/2", "π"]), 
    yticks = (-π/2:π/4:π/2, ["-π/2", "-π/4", "0", "π/4", "π/2"]),
    framestyle=:box, size=(900, 600),
    legend = false,)

# # Add red horizontal line at ϕ = 0
# hline!([0], color = :red, linewidth = 2, alpha = 0.7, legend = false, linestyle = :dot)

# # Add blue vertical line at θ = 0
# vline!([0], color = :blue, linewidth = 2, alpha = 0.7, legend = false)




# p2 = plot(legend=false, grid=true, framestyle=:box, size=(900, 600),
#           xlabel="θ = atan(y,x)", ylabel="φ = arcsin(z)",
#           xlims=(-π+0.02, π+0.02), ylims=(-π/2-0.05, π/2+0.05),
#           aspect_ratio=:equal, 
#           colorbar_title = "cos(φ)",
#           xticks = (-π:π/2:π, ["-π", "-π/2", "0", "π/2", "π"]), 
#           yticks = (-π/2:π/4:π/2, ["-π/2", "-π/4", "0", "π/4", "π/2"]),)

# Project poles
scatter!( [0.0], [π/2], mc=:green, ms=4.0)  # North pole: z/R = 1
scatter!( [0.0], [-π/2], mc=:green, ms=4.0) # South pole: z/R = -1


rings = [] # To store the projection of each ring of points

for r in 1:(levels-2)
    ring = []
    θ = r*dθ
    proj_x = Float64[]
    proj_y = Float64[]
    for ℓ in 0:(n_long-1)
        φ = ℓ*dφ
        x,y,z = to_cart(θ, φ, R_points)
        # Compute projection: (arcsin(z/R), atan(y,x))
        θ_proj = asin(z)
        φ_proj = atan(y, x)
        # Store projected points for scatter plot
        push!(proj_x, φ_proj)
        push!(proj_y, θ_proj)
        # Reporting the different values on the abcissa and ordinate axes
        push!(ring, (φ_proj, θ_proj))
    end
    push!(rings, ring)
    scatter!( proj_x, proj_y, mc=:green, ms=4.0)
end

# North and south pole
north_pole = (0.0, π/2)
south_pole = (0.0, -π/2)

color_pole = false
# Fill triangles between north pole and first ring
# for i in 1:n_long
#     point1 = north_pole
#     point2 = rings[1][i]
#     point3 = rings[1][mod1(i+1, n_long)]
#     plot!(p2, [point2[1], point2[1], point3[2], point3[1]], [point2[2], π/2, π/2, point2[2]], seriestype=:shape, fillcolor=color_pole, linecolor="green", label=false)
# end

# Fill triangles between last ring and south pole
for i in 1:n_long-1
    point1 = rings[1][i]
    point2 = rings[1][mod1(i+1, n_long)]
    plot!(p2, [point1[1], point1[1], point2[1], point2[1]], [point1[2], π/2, π/2, point2[2]], seriestype=:shape, fillcolor=color_pole, linecolor="green", label=false)
end


# Fill triangles between last ring and south pole
for i in 1:n_long
    point1 = rings[end][i]
    point2 = rings[end][mod1(i+1, n_long)]
    plot!(p2, [point1[1], point1[1], point2[1], point2[1]], [point1[2], -π/2, -π/2, point2[2]], seriestype=:shape, fillcolor=color_pole, linecolor="green", label=false)
end

# color_middle = RGBA(0.1,0.1,0.1,0.7)
# color_external = RGBA(0.4,0.4,0.4,0.7)
color_middle = false
color_external = false
# Fill triangles between rings
for r in 1:(length(rings)-1)
    ringA = rings[r]
    ringB = rings[r+1]
    for i in 1:n_long
        if abs(ringA[i][1] - ringA[mod1(i+1, n_long)][1]) < π/4.0
            a1 = ringA[i]
            a2 = ringA[mod1(i+1, n_long)]
            b1 = ringB[i]
            b2 = ringB[mod1(i+1, n_long)]
            # Two triangles per quad
            plot!(p2, [a1[1], b1[1], b2[1]], [a1[2], b1[2], b2[2]], seriestype=:shape, fillcolor=color_middle, linecolor="green", label=false)
            plot!(p2, [a1[1], b2[1], a2[1]], [a1[2], b2[2], a2[2]], seriestype=:shape, fillcolor=color_external, linecolor="green", label=false)
        end
    end
end


#Fill a fake triangles to cover the discontinuity at φ = ±π
for r in 1:(length(rings)-1)
    ringA = rings[r]
    ringB = rings[r+1]
    max = maximum([rings[r][i][1] for i in 1:n_long])
    min = minimum([rings[r][i][1] for i in 1:n_long])
    fake_positive = max + abs(ringA[1][1] - ringA[2][1])
    fake_negative = min -abs(ringA[1][1] - ringA[2][1])
    ringA_y = ringA[1][2]
    ringB_y = ringB[1][2]
    # Two triangles per quad
    plot!(p2,[max, fake_positive, fake_positive], [ringA_y, ringA_y, ringB_y], seriestype=:shape, fillcolor=color_external, linecolor="green", label=false)
    plot!(p2, [max, fake_positive, max], [ringA_y, ringB_y, ringB_y], seriestype=:shape, fillcolor=color_middle, linecolor="green", label=false)
    plot!(p2, [min, fake_negative, min], [ringA_y, ringA_y, ringB_y], seriestype=:shape, fillcolor=color_external, linecolor="green", label=false)
    plot!(p2, [fake_negative,min, fake_negative], [ringA_y, ringB_y, ringB_y], seriestype=:shape, fillcolor=color_middle, linecolor="green", label=false)
end

#Additional square on the corners 
Max = maximum([rings[1][i][1] for i in 1:n_long])
Min = minimum([rings[1][i][1] for i in 1:n_long])
fake_positive = Max + abs(rings[1][1][1] - rings[1][2][1])
fake_negative = Min -abs(rings[1][1][1] - rings[1][2][1])
plot!(p2, [Min, Min, fake_negative,fake_negative], [rings[end][1][2], -π/2, -π/2, rings[end][1][2]], seriestype=:shape, fillcolor=color_middle, linecolor="green", label=false)
plot!(p2, [Min, Min, fake_positive,fake_positive], [rings[end][1][2], -π/2, -π/2, rings[end][1][2]], seriestype=:shape, fillcolor=color_middle, linecolor="green", label=false)
plot!(p2, [Max, Max, fake_negative,fake_negative], [rings[end][1][2], π/2, π/2, rings[end][1][2]], seriestype=:shape, fillcolor=color_middle, linecolor="green", label=false)
plot!(p2, [Max, Max, fake_positive,fake_positive], [rings[end][1][2], π/2, π/2, rings[end][1][2]], seriestype=:shape, fillcolor=color_middle, linecolor="green", label=false)







# Only project the points, no grid lines or geodesics
savefig("sphere_projection.png")
