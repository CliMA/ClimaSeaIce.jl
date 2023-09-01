pushfirst!(LOAD_PATH, joinpath(@__DIR__, "..")) # add ClimaSeaIce to environment stack

using
  Documenter,
  Literate,
  # CairoMakie,  # so that Literate.jl does not capture precompilation output or warnings
  Glob,
  ClimaSeaIce

ENV["DATADEPS_ALWAYS_ACCEPT"] = "true"

#####
##### Generate examples
#####

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/generated")

example_scripts = [
    "freezing_bucket.jl",
]

for filename in example_scripts
    filepath = joinpath(EXAMPLES_DIR, filename)
    Literate.markdown(filepath, OUTPUT_DIR; flavor = Literate.DocumenterFlavor())
end

example_pages = [
    "Freezing bucket" => "generated/freezing_bucket.md",
]

#####
##### Build and deploy docs
#####

format = Documenter.HTML(
  collapselevel = 2,
     prettyurls = get(ENV, "CI", nothing) == "true",
      canonical = "https://clima.github.io/ClimaSeaIceDocumentation/dev/",
)

pages = [
    "Home" => "index.md",
    "Examples" => example_pages,

    "Library" => [ 
        "Contents"       => "library/outline.md",
        "Public"         => "library/public.md",
        "Private"        => "library/internals.md",
        "Function index" => "library/function_index.md",
        ],
]

makedocs(
   sitename = "ClimaSeaIce.jl",
    modules = [ClimaSeaIce],
     format = format,
      pages = pages,
    doctest = true,
     strict = true,
      clean = true,
  checkdocs = :exports
)

@info "Cleaning up temporary .jld2 and .nc files created by doctests..."

for file in vcat(glob("docs/*.jld2"), glob("docs/*.nc"))
    rm(file)
end

withenv("GITHUB_REPOSITORY" => "CliMA/ClimaSeaIceDocumentation") do
    deploydocs(        repo = "github.com/CliMA/ClimaSeaIceDocumentation.git",
                   versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"],
                  forcepush = true,
                  devbranch = "main",
               push_preview = true)
end

