#! julia

basepath = realpath(joinpath((@__DIR__, "..")))

using Pkg
Pkg.activate(basepath)
Pkg.instantiate()

##

using LiveServer

run = true
while run
    include("make.jl")

    serve(dir=joinpath(@__DIR__, "build"))

    println("Run again? Enter! Exit witn 'q'.")
    if readline() == "q"
        global run = false
    end
end