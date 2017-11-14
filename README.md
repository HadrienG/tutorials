# Tutorials

This repository contains basic tutorials and walkthroughs on various
bioinformatics subjects:

## Dev

First install mkdocs in a virtual environment

```
virtualenv mkdocs_env
source mkdocs_env/bin/activate
pip install mkdocs
```

then the theme

```
pip install mkdocs-material
```

Then clone the directory

```
git clone https://github.com/SGBC/cluster_doc.git
cd cluster_doc
```

For a live preview in your browser do

```
mkdocs serve &
```

## Deploy

The following command will build and push your website to a `gh-pages` branch.
Only do this if you want your own version of the website! If you are modifying
the original, please open a pull request.

```
mkdocs gh-deploy
```
