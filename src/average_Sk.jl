import Pkg; Pkg.activate(".")


using DelimitedFiles
using CairoMakie
using Statistics

files = filter(x -> !(contains(x, "Equilibration")), readdir("Processed_Data/Sk/"))
for rho in [0.8, 0.94]
    files_rho = filter(x -> contains(x, "rho_$(rho)"), files)
    Sk_total = []
    k = []

    for file in files_rho
        seed = parse(Int, split(split(file, "seed_")[2], ".")[1])
        savefile = "Processed_Data/Sk/Sk_rho_$(rho)_seed_$(seed).txt"
        if isfile(savefile)
            data = readdlm(savefile)
            push!(Sk_total, data[:, 2])
            push!(k, data[:, 1])
            println("Loaded $(savefile)")
        else
            println("File not found: $savefile")
        end
    end

    Sk_total_mean = sum(Sk_total; init=zeros(length(Sk_total[1]))) ./ length(Sk_total)
    std_of_mean_Sk = std(reduce(hcat, Sk_total); dims=2) ./ sqrt(length(Sk_total))

    @assert allequal(k)
    k_mean = k[1]
    println("averaged Sk for ρ = $(rho) for $(length(Sk_total)) seeds")
    open("Processed_Data/Sk/mean/Sk_rho_$(rho).txt", "w") do io
        writedlm(io, [k_mean Sk_total_mean std_of_mean_Sk])
    end
end

f = Figure(size=(1000, 1000))
ax = Axis(f[1, 1], xlabel="k", ylabel="S(k)")
for rho in [0.8, 0.94]
    data = readdlm("Processed_Data/Sk/mean/Sk_rho_$(rho).txt")
    k = data[:, 1]
    Sk = data[:, 2]
    lines!(ax, k, Sk, label="ρ = $(rho)")
end
axislegend(ax)
#zoom in on the low k behaviour

ax2 = Axis(f[1, 2], xlabel="k", ylabel="S(k)", limits=(0,2,0.015,0.045), title="Zoomed in on low k")

for rho in [0.8, 0.94]
    data = readdlm("Processed_Data/Sk/mean/Sk_rho_$(rho).txt")
    k = data[:, 1]
    Sk = data[:, 2]
    scatterlines!(ax2, k, Sk, label="ρ = $(rho)")
end

#zoom of peak

ax3 = Axis(f[2, 1], xlabel="k", ylabel="S(k)", limits=(6,8,1.5,3.0), title="zoomed in on peak")

for rho in [0.8, 0.94]
    data = readdlm("Processed_Data/Sk/mean/Sk_rho_$(rho).txt")
    k = data[:, 1]
    Sk = data[:, 2]
    scatterlines!(ax3, k, Sk, label="ρ = $(rho)")
end

#zoom of second peak

ax4 = Axis(f[2, 2], xlabel="k", ylabel="S(k)", limits=(11,15,1.0,1.5), title="Zoomed in on second peak")

for rho in [0.8, 0.94]
    data = readdlm("Processed_Data/Sk/mean/Sk_rho_$(rho).txt")
    k = data[:, 1]
    Sk = data[:, 2]
    scatterlines!(ax4, k, Sk, label="ρ = $(rho)")
end

display(f)


function c_from_s(Sk, rho)
    ρCk = @. (Sk-1)/Sk
    return ρCk/rho
end

f = Figure(size=(1000, 1000))
ax = Axis(f[1, 1], xlabel="k", ylabel="ρ c(k)")
for rho in [0.8, 0.94]
    data = readdlm("Processed_Data/Sk/mean/Sk_rho_$(rho).txt")
    k = data[:, 1]
    Sk = data[:, 2]
    lines!(ax, k, rho*c_from_s(Sk, rho), label="ρ = $(rho)")
end
axislegend(ax)

#zoom in on the low k behaviour

ax2 = Axis(f[1, 2], xlabel="k", ylabel="ρ c(k)", limits=(0.0,1,-59, -47), title="Zoomed in on low k")

for rho in [0.8, 0.94]
    data = readdlm("Processed_Data/Sk/mean/Sk_rho_$(rho).txt")
    k = data[:, 1]
    Sk = data[:, 2]
    scatterlines!(ax2, k, rho*c_from_s(Sk, rho), label="ρ = $(rho)")
end

#zoom of peak

ax3 = Axis(f[2, 1], xlabel="k", ylabel="ρ c(k)", limits=(6,8,0,0.7), title="zoomed in on peak")

for rho in [0.8, 0.94]
    data = readdlm("Processed_Data/Sk/mean/Sk_rho_$(rho).txt")
    k = data[:, 1]
    Sk = data[:, 2]
    scatterlines!(ax3, k, rho*c_from_s(Sk, rho), label="ρ = $(rho)")
end

#zoom of second peak

ax4 = Axis(f[2, 2], xlabel="k", ylabel="ρ c(k)", limits=(11,15,0,0.3), title="Zoomed in on second peak")

for rho in [0.8, 0.94]
    data = readdlm("Processed_Data/Sk/mean/Sk_rho_$(rho).txt")
    k = data[:, 1]
    Sk = data[:, 2]
    scatterlines!(ax4, k, rho*c_from_s(Sk, rho), label="ρ = $(rho)")
end

display(f)