using Test
using ClimaSeaIce
using Aqua: Aqua
using ExplicitImports: ExplicitImports

function walk_submodules!(result, visited, mod::Module)
    for name in sort(names(mod; all=true, imported=false))
        isdefined(mod, name) || continue
        value = getproperty(mod, name)
        if value isa Module &&
            parentmodule(value) === mod &&
            !(value in visited) &&
            value !== mod

            push!(visited, value)
            push!(result, value)
            walk_submodules!(result, visited, value)
        end
    end
end
function get_submodules(mod::Module; self=true)
    result = self ? Module[mod] : Module[]
    visited = Set{Module}()

    walk_submodules!(result, visited, mod)
    return result
end

@testset "Aqua" begin
    Aqua.test_all(ClimaSeaIce; piracies=false)

    # `test_piracies` doesn't recurse in inner modules, so we have to test that manually.
    @testset "No type piracy in $(mod)" for mod in get_submodules(ClimaSeaIce)
        Aqua.test_piracies(mod)
    end
end

@testset "Explicit Imports" begin
    @test ExplicitImports.check_no_implicit_imports(ClimaSeaIce) === nothing
end

@testset "Import via Owner" begin
    @test ExplicitImports.check_all_explicit_imports_via_owners(ClimaSeaIce) === nothing
end

@testset "Stale Explicit Imports" begin
    @test ExplicitImports.check_no_stale_explicit_imports(ClimaSeaIce) === nothing
end

@testset "Qualified Accesses" begin
    @test ExplicitImports.check_all_qualified_accesses_via_owners(ClimaSeaIce) === nothing
end

@testset "Self Qualified Accesses" begin
    @test ExplicitImports.check_no_self_qualified_accesses(ClimaSeaIce) === nothing
end
