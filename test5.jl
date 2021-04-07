using Arpack
using DelimitedFiles
using FastGaussQuadrature
using LinearAlgebra
using MatrixMarket
using SparseArrays
using SuiteSparse
using Printf


cd(@__DIR__)
include("algorithms/action_cauchy.jl")
include("algorithms/action_de.jl")
include("algorithms/action_gj1.jl")
include("algorithms/action_gj2.jl")
include("algorithms/action_gj2pre.jl")


function test5_spd()
    matname_list = ["poisson200", "s2rmt3m1", "fv3", "poisson200"]
    αx10_list = [2, 8]
    Data = Array{Any, 2}(undef, 9, 12)
    Data .= ""
    Data[1, 1:2] .= ["matrix", "α"]
    Data[1, 3:2:11] .= ["gj1", "gj2", "gj2pre", "de", "cauchy"]
    i = 1
    for matname in matname_list
        A = mmread("matrix/$(matname).mtx")
        n = size(A, 1)
        b = reshape(readdlm("matrix/test5_b_$(matname).txt"), n)
        ϵ = 1e-6
        for αx10 in αx10_list
            println("--- $(matname), α = 0.$(αx10) ---")
            i += 1
            α = αx10 / 10

            Data[i, 1] = matname
            Data[i, 2] = "0.$(αx10)"

            sttime = time_ns()
            f, m = powmvec_gj1_spd(A, b, α, ϵ=ϵ)
            entime = time_ns()
            cputime = (entime - sttime) * 1e-9
            Data[i, 3] = cputime
            Data[i, 4] = m
            @printf("   gj1: %6.2fsec (m = %4d)\n", cputime, m)

            sttime = time_ns()
            f, m = powmvec_gj2_spd(A, b, α, ϵ=ϵ)
            entime = time_ns()
            cputime = (entime - sttime) * 1e-9
            Data[i, 5] = cputime
            Data[i, 6] = m
            @printf("   gj2: %6.2fsec (m = %4d)\n", cputime, m)


            sttime = time_ns()
            f, m = powmvec_gj2pre_spd(A, b, α, ϵ=ϵ)
            entime = time_ns()
            cputime = (entime - sttime) * 1e-9
            Data[i, 7] = cputime
            Data[i, 8] = m
            @printf("gj2pre: %6.2fsec (m = %4d)\n", cputime, m)


            sttime = time_ns()
            f, m = powmvec_de_spd(A, b, α, ϵ=ϵ)
            entime = time_ns()
            cputime = (entime - sttime) * 1e-9
            Data[i, 9] = cputime
            Data[i, 10] = m
            @printf("    de: %6.2fsec (m = %4d)\n", cputime, m)


            sttime = time_ns()
            f, m = powmvec_cauchy_spd(A, b, α, ϵ=ϵ)
            entime = time_ns()
            cputime = (entime - sttime) * 1e-9
            Data[i, 11] = cputime
            Data[i, 12] = m
            @printf("cauchy: %6.2fsec (m = %4d)\n\n", cputime, m)
        end
    end

    writedlm("result/test5_spd.csv", Data, ',')
end


function test5_general()
    matname_list = ["cell1", "cell1", "TSOPF_RS_b9_c6", "circuit_3"]
    αx10_list = [2, 8]
    Data = Array{Any, 2}(undef, 9, 8)
    Data .= ""
    Data[1, 1:2] .= ["matrix", "α"]
    Data[1, 3:2:7] .= ["gj1", "gj2", "de"]
    i = 1
    for matname in matname_list
        A = mmread("matrix/$(matname).mtx")
        if matname == "cell1"
            A = A + 5e-8*I
        elseif matname == "TSOPF_RS_b9_c6"
            A = A + 40.35*I
        elseif matname == "circuit_3"
            A = A + 3*I
        end

        n = size(A, 1)
        b = reshape(readdlm("matrix/test5_b_$(matname).txt"), n)
        ϵ = 1e-6
        for αx10 in αx10_list
            println("--- $(matname), α = 0.$(αx10) ---")
            i += 1
            α = αx10 / 10

            Data[i, 1] = matname
            Data[i, 2] = "0.$(αx10)"

            sttime = time_ns()
            f, m = powmvec_gj1(A, b, α, ϵ=ϵ)
            entime = time_ns()
            cputime = (entime - sttime) * 1e-9
            Data[i, 3] = cputime
            Data[i, 4] = m
            @printf("   gj1: %6.2fsec (m = %4d)\n", cputime, m)

            sttime = time_ns()
            f, m = powmvec_gj2(A, b, α, ϵ=ϵ)
            entime = time_ns()
            cputime = (entime - sttime) * 1e-9
            Data[i, 5] = cputime
            Data[i, 6] = m
            @printf("   gj2: %6.2fsec (m = %4d)\n", cputime, m)

            sttime = time_ns()
            f, m = powmvec_de(A, b, α, ϵ=ϵ)
            entime = time_ns()
            cputime = (entime - sttime) * 1e-9
            Data[i, 7] = cputime
            Data[i, 8] = m
            @printf("    de: %6.2fsec (m = %4d)\n", cputime, m)

        end
    end
    writedlm("result/test5_general.csv", Data, ',')
end