
cd(@__DIR__)

println("Test 1")
include("test1.jl")
generate_exact_solutions()
main()


println("\n\nTest 2")
include("test2.jl")
generate_exact_solutions()
main()


println("\n\nTest 3")
include("test3.jl")
generate_exact_solutions()
main()


println("\n\nTest 4")
include("test4.jl")
generate_exact_solutions()
main()


println("\n\nTest 5")
include("test5.jl")
test5_spd()
test5_general()
