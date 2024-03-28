using Luxor

Drawing(280, 280, "logo.png")
origin()

colors = (Luxor.julia_green, Luxor.julia_red, Luxor.julia_purple)
θ = [π / 6, 5π / 6, 3π / 2]
r = 100
centers = [Point(x, y) for (x, y) in zip(r * cos.(θ), r * sin.(θ))]

for (color, center) in zip(colors, centers)
    setcolor(color)
    circle(center, 40, :fill)
end

pts = ngon(Point(0, 0), 80, 21, vertices = true)
sethue(Luxor.julia_blue)
arrow(
    pts[4:7]...,
    :stroke,
    arrowheadlength = 35,
    linewidth = 10,
    startarrow = false,
    finisharrow = true,
)
arrow(
    pts[11:14]...,
    :stroke,
    arrowheadlength = 35,
    linewidth = 10,
    startarrow = false,
    finisharrow = true,
)
arrow(
    pts[18:21]...,
    :stroke,
    arrowheadlength = 35,
    linewidth = 10,
    startarrow = false,
    finisharrow = true,
)
finish()
preview()
