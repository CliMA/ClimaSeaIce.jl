using Test
using ClimaSeaIce
using ExplicitImports: ExplicitImports

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
