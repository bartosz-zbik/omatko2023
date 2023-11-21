
function generate_tjectory(x0::Real, nsteps::Integer, dt::Real)
        trajectory = Vector{Float64}(undef, nsteps)
        trajectory[1] = x0
        increment = sqrt(dt)
        for i in 2:length(trajectory)
                trajectory[i] = trajectory[i-1] + increment * randn()
        end
        return trajectory
end

function generate_tjectory!(trajectory::AbstractArray, x0::Real, dt::Real)
        trajectory[1] = x0
        increment = sqrt(dt)
        for i in 2:length(trajectory)
                trajectory[i] = trajectory[i-1] + increment * randn()
        end
        return trajectory
end

gt = generate_tjectory
gt! = generate_tjectory!

