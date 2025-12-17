# Graphical output 
# This might belong in DrawAndSpray.jl,
# but also might not be worth the bloat.
# It's similar to mark_at!, but single-pixel only.

"""
    set_pixel_from_palette!(buf, pt, i::Int64,  palette)

If i > length(palette) will cycle back.
"""
function set_pixel_from_palette!(buf, pt, i::Int64,  palette)
    # Cycle if palette contains "too few colors"
    colono = mod(i - 1, length(palette)) + 1
    buf[pt] = palette[colono]
end
