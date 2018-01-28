# Tutorials

This repository contains basic tutorials and walkthroughs on various
bioinformatics subjects:

## Dev

First install mkdocs in a virtual environment

```bash
virtualenv mkdocs_env
source mkdocs_env/bin/activate
pip install mkdocs
```

then the theme

```bash
pip install mkdocs-material
```

Then clone the directory

```bash
git clone https://github.com/HadrienG/tutorials.git
cd tutorials
```

For a live preview in your browser do

```bash
mkdocs serve &
```

## Deploy

The following command will build and push your website to a `gh-pages` branch.
Only do this if you want your own version of the website! If you are modifying
the original, please open a pull request.

```bash
mkdocs gh-deploy
```
