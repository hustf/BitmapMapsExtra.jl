# Callees for 'spray_streamlines`
# concerning extraction of values from diffeq solutions.
#

################################################
# Sample streamlines from continuous to discrete  
################################################

function extract_discrete_points_on_streamline(sol, negy::NegateY, sol_density)
    if ! (sol.retcode == Success ||
        sol.retcode == Terminated) # || sol.retcode == DtLessThanMin)
        @warn sol.retcode
    end
    # Pre-allocate
    pts = CartesianIndex{2}[]
    # Circular buffer of recent points (as (i,j) tuples)
    buf = @MVector [(0, 0) for _ in 1:100]
    buf_len = 0          # how many valid entries we currently track (≤ 100)
    buf_idx = 1          # next position to overwrite (1-based, circular)
    # sol_density determines how far, in solution time, between examined solution points.
    trng = range(first(sol.t), last(sol.t), step = sol_density * sign( last(sol.t) - first(sol.t)   ))
    for t in trng
        x, y = sol(t) # Fast, interpolated lookup 
        ny = negy(y)
        i = Int(round(ny))
        j = Int(round(x))
        p = (i, j)
        # Skip if in recent buffer
        if !any(==(p), @view buf[1:buf_len])
            # Add to pts and to recent buffer
            push!(pts, CartesianIndex(p...))
            @inbounds buf[buf_idx] = p
            buf_idx = buf_idx == 100 ? 1 : buf_idx + 1
            buf_len = min(buf_len + 1, 100)
         end
    end
    pts
end