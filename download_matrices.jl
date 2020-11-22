using Pkg
Pkg.add("BinDeps")

using BinDeps


url_list = [
    "https://suitesparse-collection-website.herokuapp.com/MM/FIDAP/ex5.tar.gz",
    "https://suitesparse-collection-website.herokuapp.com/MM/HB/pores_1.tar.gz",
    "https://suitesparse-collection-website.herokuapp.com/MM/HB/nos4.tar.gz",
    "https://suitesparse-collection-website.herokuapp.com/MM/HB/bcsstk04.tar.gz",
    "https://suitesparse-collection-website.herokuapp.com/MM/HB/lund_b.tar.gz",
    "https://suitesparse-collection-website.herokuapp.com/MM/Cylshell/s2rmt3m1.tar.gz",
    "https://suitesparse-collection-website.herokuapp.com/MM/Norris/fv3.tar.gz",
    "https://suitesparse-collection-website.herokuapp.com/MM/Lucifora/cell1.tar.gz",
    "https://suitesparse-collection-website.herokuapp.com/MM/TSOPF/TSOPF_RS_b9_c6.tar.gz",
    "https://suitesparse-collection-website.herokuapp.com/MM/Bomhof/circuit_3.tar.gz"
]

matname_list = [split(split(url, "/")[end], ".")[1] for url in url_list]

for (url, matname) in zip(url_list, matname_list)
    println("$(matname)")
    println("\tdownloading...")
    Base.download(url, "matrix/$(matname).tar.gz")
    println("\textracting...")
    run(unpack_cmd("matrix/$(matname).tar.gz", "matrix", ".tar", ".gz"))
    mv("matrix/$(matname)/$(matname).mtx", "matrix/$(matname).mtx", force=true)
    rm("matrix/$(matname)", force=true, recursive=true)
    rm("matrix/$(matname).tar.gz", force=true, recursive=true)
end
