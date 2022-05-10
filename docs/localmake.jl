#! julia

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
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