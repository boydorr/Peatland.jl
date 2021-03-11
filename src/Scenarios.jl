using Simulation.Units
using Distributions

function Invasive(eco::Ecosystem, timestep::Unitful.Time, rate::typeof(1.0/day))
    qual = eco.spplist.native .== false
    invasive = findall(qual)
    natives = findall(.!qual)
    eco.spplist.traits.var[invasive] .= maximum(eco.spplist.traits.var)
    invasive_abun = eco.spplist.abun[invasive]
    avgain = uconvert(NoUnits, rate * timestep)
    for i in eachindex(invasive)
        gains = rand(Poisson(avgain))
        pos = 1:size(eco.abenv.habitat.matrix, 1)
        smp = sample(pos, gains)
        eco.abundances.grid[invasive[i], :, :] .+=  map(x -> sum(smp .== x), pos)
    end
end
