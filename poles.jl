#using LinearAlgebra
using Roots
using DelimitedFiles
#using BenchmarkTools

# Constants
ħ = 1.0
m = 0.5

# Change directory if need be
listparams = readdlm("\\Users\\rigas\\Desktop\\Physics\\Διδακτορικό\\BCS Green functions\\Part II - Superconducting Fermi Liquid w SOC\\Code\\config.txt")

# Variables
k_F = listparams[1]
μ = (ħ*k_F)^2/(2m)
Δ = listparams[2]
α = listparams[3]
Β = listparams[4]

# Energy and η
E = listparams[5]
η = 0.001

ϵ(k) = (ħ*k)^2/(2m) - μ     # Single electron energy function

function energypoles(k,pm)
    return E - sqrt((ϵ(k))^2 + (α*k)^2 + Β^2 + Δ^2 + 2pm*sqrt((ϵ(k)*α*k)^2 + Β^2*(Δ^2 + (ϵ(k))^2)))
end

Eq_1(k) = energypoles(k,1)
Eq_2(k) = energypoles(k,-1)

# Finds the k-poles for the specific energy in order to be used for the integration
listzeros = zeros(0) # Initialize a 1-dimensional array
append!(listzeros, find_zeros(Eq_1,-10,10^10))
append!(listzeros, find_zeros(Eq_2,-10,10^10))

# Change directory if need be
file2 = open("\\Users\\rigas\\Desktop\\Physics\\Διδακτορικό\\BCS Green functions\\Part II - Superconducting Fermi Liquid w SOC\\Code\\kpoles.dat", "w")
for i in 1:length(listzeros)
    if listzeros[i] > 0.0
        println(file2, listzeros[i])
    end
end
close(file2);
