using
  Documenter,
  Literate,
  ClimaSeaIce

ENV["DATADEPS_ALWAYS_ACCEPT"] = "true"

#####
##### Generate examples
#####

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/literated")

example_scripts = [
    "freezing_bucket.jl",
    "melting_in_spring.jl",
    "freezing_of_a_lake.jl",
    "ice_advected_by_anticyclone.jl",
    # "ice_advected_on_coastline.jl",
    "arctic_basin_seasonal_cycle.jl"
]

for filename in example_scripts
    filepath = joinpath(EXAMPLES_DIR, filename)
    Literate.markdown(filepath, OUTPUT_DIR; flavor = Literate.DocumenterFlavor())
end

example_pages = [
    "Freezing bucket" => "literated/freezing_bucket.md",
    "Melting in Spring" => "literated/melting_in_spring.md",
    "Freezing of a Lake" => "literated/freezing_of_a_lake.md",
    "Ice advected by anticyclone" => "literated/ice_advected_by_anticyclone.md",
    # "Ice advected on coastline" => "literated/ice_advected_on_coastline.md",
    "Arctic basin seasonal cycle" => "literated/arctic_basin_seasonal_cycle.md"
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
    warnonly = [:cross_references],
       clean = true,
   checkdocs = :exports
)

@info "Clean up temporary .jld2/.nc files created by doctests..."

"""
    recursive_find(directory, pattern)

Return list of filepaths within `directory` that contains the `pattern::Regex`.
"""
recursive_find(directory, pattern) =
    mapreduce(vcat, walkdir(directory)) do (root, dirs, files)
        joinpath.(root, filter(contains(pattern), files))
    end

files = []
for pattern in [r"\.jld2", r"\.nc"]
    global files = vcat(files, recursive_find(@__DIR__, pattern))
end

for file in files
    rm(file)
end

# Replace with below once https://github.com/JuliaDocs/Documenter.jl/pull/2692 is merged and available.
#  deploydocs(repo = "github.com/CliMA/ClimaSeaIce.jl",
#    deploy_repo = "github.com/CliMA/ClimaSeaIceDocumentation",
#    devbranch = "main",
#    push_preview = true)
if get(ENV, "GITHUB_EVENT_NAME", "") == "pull_request"
    deploydocs(repo = "github.com/CliMA/ClimaSeaIce.jl",
               repo_previews = "github.com/CliMA/ClimaSeaIceDocumentation",
               devbranch = "main",
               forcepush = true,
               push_preview = true,
               versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"])
else
    repo = "github.com/CliMA/ClimaSeaIceDocumentation"
    withenv("GITHUB_REPOSITORY" => repo) do
        deploydocs(repo = repo,
                   devbranch = "main",
                   forcepush = true,
                   versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"])
    end
end
