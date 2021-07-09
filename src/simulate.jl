"""
simulate(::Simulation)

Simulation a TVVAR model.

The inputs m and Σ can be time-varying.

"""

function simulate(sim::Simulation)
    m, Σ = sim.m, sim.Σ
    N, d, p = size(m)
    println(N, " ", p, " ", d)
    d = convert(Integer, (d - 1)/p)
    y = zeros(N, p)
    for t in 1:N
        if t ≤ d
            y[1:d, :] = rand(MvNormal(m[t, 1, :], Σ[t, :, :]))
        else
            if d == 0
                F = getF()
            else
                F = getF(y[t - d:t - 1, :])
            end
            # Simulate
            y[t, :] = rand(MvNormal(m[t, :, :]' * F, Σ[t, :, :]))
        end
    end
    return y
end

sim = Simulation([5.0 0.0; 0.9 0.0; 0.0 0.0], [1.0 0.0; 0.0 2.0], 500);
sim = simulate(sim);

using Plots
plot(sim)

m = vcat(expandArray([1.0 1.0], 200), expandArray([0.0 10.0], 300));
Σ = vcat(expandArray([1.0 0.0; 0.0 2.0], 300), expandArray([5.0 0.0; 0.0 1.0], 200));
sim = Simulation(m, Σ);
sim = simulate(sim);

using Plots
plot(sim)