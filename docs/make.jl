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
    "ice_advected_by_anticyclone.jl",
    "ice_advected_on_coastline.jl",
]

for filename in example_scripts
    filepath = joinpath(EXAMPLES_DIR, filename)
    Literate.markdown(filepath, OUTPUT_DIR; flavor = Literate.DocumenterFlavor())
end

example_pages = [
    "Freezing bucket" => "literated/freezing_bucket.md",
    "Ice advected by anticyclone" => "literated/ice_advected_by_anticyclone.md",
    "Ice advected on coastline" => "literated/ice_advected_on_coastline.md",
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

deploydocs(repo = "github.com/CliMA/ClimaSeaIceDocumentation.git",
           versions = ["stable" => "v^", "dev" => "dev", "v#.#.#"],
           forcepush = true,
           push_preview = true,
           devbranch = "main")