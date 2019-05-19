import sys
import turtle

name = sys.stdin.readline().split()[2]
print(name)

coordinates = []
with open("ExampleProblems/" + name + ".tsp") as tsplib_file:
    for line in tsplib_file:
        elements = line.strip().split(" ")
        try:
            index, x, y = [float(e) for e in elements]
        except ValueError:
            continue
        coordinates.append([x, y])

minX = min([x for (x, y) in coordinates])
maxX = max([x for (x, y) in coordinates])
minY = min([y for (x, y) in coordinates])
maxY = max([y for (x, y) in coordinates])

# print(f"{minX} - {maxX}, {minY} - {maxY}")

turtle.setup(800, 800, 0, 0)
turtle.setworldcoordinates(min(minX, minY) - 10, min(minX, minY) - 10, max(maxX, maxY) + 10, max(maxX, maxY) + 10)
turtle.shape('circle')
turtle.shapesize(0.2, 0.2, 0)
turtle.speed(0)
turtle.tracer(100000)
turtle.penup()
turtle.goto(coordinates[0])
turtle.pendown()

while True:
    line = sys.stdin.readline()
    if not line.startswith(" new tour: "):
        print(line, end="")
    if line.startswith(" new tour: "):
        turtle.clear()
        line = line[len(" new tour: "):]
        tour = [int(v) for v in line.split(", ")]
        for vertex in tour:
            turtle.goto(coordinates[vertex - 1])
            turtle.stamp()
        turtle.goto(coordinates[0])
        turtle.update()
    elif line.startswith("The best tour found by the heuristic"):
        break
