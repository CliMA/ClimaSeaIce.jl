```@eval
using Markdown
using ClimaSeaIce

readme_path = joinpath(dirname(dirname(pathof(ClimaSeaIce))), "README.md")
Markdown.parse(read(readme_path, String))
```
