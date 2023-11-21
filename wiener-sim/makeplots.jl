using CairoMakie
using DelimitedFiles

set_my_theme!()

function CV_wiener(x0, L)
    return sqrt(2/3 * (L^2 + x0^2)/(L^2 - x0^2))
end

function MFPT_wiener(x0, L, s, r)
    k = sqrt(s^2 / r)
    d = sinh((L-x0) / k) + sinh((L+x0) / k)
    return 1/r * (sinh(2*L/k) / d - 1)
end
function MFPT_wiener_N(x0, L, s)
    k = L^2 - x0^2
    return k / s^2 / 2
end

function MFPT_min_norm(x0)
    all_p = [MFPT_wiener(x0, 1, 1, r)
                 for r in range(1e-7, 10, length=100)]
    push!(all_p, MFPT_wiener_N(x0, 1, 1))
    m = minimum(all_p)
    return m / MFPT_wiener_N(x0, 1, 1)
end

# manytr = [generate_tjectory(0.0, 1000, 1e-4) for i in 1:100]
manytr = readdlm("manytr.txt")


fig, ax = lines(manytr[:, 1], color=Cycled(1))
lines!(manytr[:, 3], color=Cycled(2))
lines!(manytr[:, 7], color=Cycled(3))
lines!(manytr[:, 9], color=Cycled(4))
ax.xlabel = L"$t$ - czas"
ax.ylabel = L"W_t"
save("only-tr.pdf", fig)


hlines!(ax, [0.1, -0.6], color=:black, linestyle=:dash, alpha=0.7, linewidth=3)
vlines!(ax, [findfirst(x -> x >= 0.1, manytr[:, 7])],
       color=Cycled(3), linestyle=:dash)
vlines!(ax, [findfirst(x -> x >= 0.1, manytr[:, 9])],
       color=Cycled(4), linestyle=:dash)

save("tr-with-bounds.pdf", fig)


fig, ax = lines(manytr[:, 3], color=Cycled(2))
lines!(1:250, manytr[1:250, 3], color=Cycled(1))
lines!(251:710, manytr[251:710, 3] .- manytr[251, 3], color=Cycled(1))
lines!(711:900, manytr[711:900, 3] .- manytr[711, 3], color=Cycled(1))
lines!(901:1000, manytr[901:1000, 3] .- manytr[901, 3], color=Cycled(1))

hlines!(ax, [0.1, -0.6], color=:black, linestyle=:dash, alpha=0.7, linewidth=3)

arrows!([250, 710, 900],
        [manytr[250, 3], manytr[710, 3] - manytr[251, 3], manytr[900, 3] .- manytr[711, 3]],
        [1, 1, 1],
        -[manytr[250, 3], manytr[710, 3] - manytr[251, 3], manytr[900, 3] .- manytr[711, 3]])

vlines!(ax, [710 + findfirst(x -> x >= 0.1,
                       manytr[711:900, 3] .- manytr[711, 3])],
       color=Cycled(1), linestyle=:dash)

save("tr-resetting.pdf", fig)

fig, ax = lines(0..0.7, MFPT_min_norm, label=L"\Lambda(x_0)",
                figure = (; resolution = (400, 300)))
lines!(ax, 0..0.7, x -> CV_wiener(x, 1), label=L"CV(x_0)")
ax.limits = (nothing, (0.7, 1.25))
ax.xlabel = L"x_0"
axislegend(ax, position = :lt)
save("CV-MFPT.pdf", fig)

