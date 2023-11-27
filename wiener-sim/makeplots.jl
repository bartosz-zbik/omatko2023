using CairoMakie
using DelimitedFiles

mytheme = merge(Theme(
        Lines=(cycle=Cycle([:color, :linestyle], covary=true),),
        Scatter=(cycle=Cycle([:color, :marker], covary=true),)
    ), theme_latexfonts())
set_theme!(mytheme)

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


fig, ax = lines(manytr[:, 1], color=Cycled(1), linewidth=2)
lines!(manytr[:, 3], color=Cycled(2), linewidth=2)
lines!(manytr[:, 7], color=Cycled(3), linewidth=2)
lines!(manytr[:, 9], color=Cycled(4), linewidth=2)
ax.xlabel = L"$t$ - czas"
ax.ylabel = L"W_t"
save("only-tr.pdf", fig)


hlines!(ax, [0.1, -0.6], color=:black, linestyle=:dash, alpha=0.7, linewidth=3)
vlines!(ax, [findfirst(x -> x >= 0.1, manytr[:, 7])],
       color=Cycled(3), linestyle=:dash, linewidth=2, alpha=0.9)
vlines!(ax, [findfirst(x -> x >= 0.1, manytr[:, 9])],
       color=Cycled(4), linestyle=:dash, linewidth=2, alpha=0.9)

save("tr-with-bounds.pdf", fig)


fig, ax = lines(manytr[:, 3], color=Cycled(2), linewidth=2)
lines!(1:250, manytr[1:250, 3], color=Cycled(1), linewidth=2)
lines!(251:710, manytr[251:710, 3] .- manytr[251, 3], color=Cycled(1), linewidth=2)
lines!(711:900, manytr[711:900, 3] .- manytr[711, 3], color=Cycled(1), linewidth=2)
lines!(901:1000, manytr[901:1000, 3] .- manytr[901, 3], color=Cycled(1), linewidth=2)

hlines!(ax, [0.1, -0.6], color=:black, linestyle=:dash, alpha=0.7, linewidth=4)

arrows!([250, 710, 900],
        [manytr[250, 3], manytr[710, 3] - manytr[251, 3], manytr[900, 3] .- manytr[711, 3]],
        [1, 1, 1],
        -[manytr[250, 3], manytr[710, 3] - manytr[251, 3], manytr[900, 3] .- manytr[711, 3]], linewidth=2, arrowsize=15)

vlines!(ax, [710 + findfirst(x -> x >= 0.1,
                       manytr[711:900, 3] .- manytr[711, 3])],
       color=Cycled(1), linestyle=:dash, linewidth=3, alpha=0.8)

save("tr-resetting.pdf", fig)

hspan!(ax, -0.5, 0.0, color = (:green, 0.2))
hspan!(ax, -0.6, -0.5, color = (:pink, 0.2))
hspan!(ax, 0.0, 0.1, color = (:pink, 0.2))
save("tr-resetting-reg.pdf", fig)


fig, ax = lines(0..0.7, MFPT_min_norm, label=L"\Lambda(x_0)",
                linewidth=3,
                figure = (; resolution = (400, 300)))
lines!(ax, 0..0.7, x -> CV_wiener(x, 1), label=L"CV(x_0)",
      linewidth=3)
ax.limits = (nothing, (0.7, 1.25))
ax.xlabel = L"x_0"
axislegend(ax, position = :lt, patchsize = (50, 25))
save("CV-MFPT.pdf", fig)


with_theme(mytheme, fontsize = 25) do
fig, ax = lines(6e-4..3, x -> MFPT_wiener(0.25, 0.35, 1, 1/x),
    linewidth=3, label="MFPT(r)")
hlines!(ax, [MFPT_wiener_N(0.25, 0.35, 1)], linewidth=3, color=:magenta, linestyle=:dash, label="MFPT, bez resetowania")
ax.xlabel = L"\mathbb{E} R"
ax.ylabel="MFPT"
ax.title = "Ucieczka z (-0.6, 0.1), proces Wienera W₀≡0"
axislegend(ax, position = :rb, patchsize = (50, 25))
ax.xscale = log10

save("MFPT-r.pdf", fig)
end

