---
title: FLPKG NOTES
author: Iago Mosqueira, EC JRC G03
license: Creative Commons Attribution-ShareAlike 4.0 International
---

# Package structure

- README.md
- NEWS.md
- .travis.yml
- .Rbuildignore
- Makefile
- CITATION

## Vignette

Create at least one vignette, presenting some example uses of the package capabilities and, if neccessary, briefly documenting the new classes and methods in the package.

- [rmarkdown package_vignette format](http://rmarkdown.rstudio.com/package_vignette_format.html)
- [rmarkdown html_document format](http://rmarkdown.rstudio.com/html_document_format.html)
- [Bibliographies and Citations](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)

## /data

## tests/testthat

# WORKFLOWS

## Create new package (or make existing one follow this style)

- EDIT DESCRIPTION
- ADD NAMESPACE
- ADD R files following templates for roxygen documentation
- EDIT README.md
- EDIT inst/CITATION

## Add new code / fix a bug

- CHANGE code in R files
- ADD entry in NEWS.md
- RUN 'make test' and check test passes
- git add -A .
- git commit -m 'Message used in NEWS.md'
- EDIT Version and Date in DESCRIPTION
- RUN 'make README.md'
- git add .
- git commit -m 'Version `sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION`'
- git push

## Add or alter roxygen documentation

- CHANGE entry in R file(s)
- Run 'make roxygen'
- Run 'make gh-pages'
- git add -A .
- git commit -m 'Updated documentation for ...'
- git push


# HOW TO convert pkgs to use gh-pages using staticdocs

1. ADD README.md following template
2. COPY Makefile
3. CREATE an empty gh-pages branch

	git checkout --orphan gh-pages
	git rm -rf .

4. ADD a README.md to the branch

	touch README.md
	git add README.md

5. COMMIT and PUSH branch
	
	git commit -a -m "First gh-pages commit"
	git push origin gh-pages

# PKG CHECKLIST

- COPY folder and rename as per package name
- EDIT DESCRIPTION
    - Version: You can follow [Hadley Wickam's format](http://r-pkgs.had.co.nz/description.html#version)
    - 

# MORE INFO

- [R packages](http://r-pkgs.had.co.nz/)
