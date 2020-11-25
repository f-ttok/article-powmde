
cd(@__DIR__)

include("test1.jl")
generate_exact_solutions()
main()


include("test2.jl")
generate_exact_solutions()
main()


include("test3.jl")
generate_exact_solutions()
main()


include("test4.jl")
generate_exact_solutions()
main()


include("test5.jl")
main()
