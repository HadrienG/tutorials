# Tutorials

This repository contains basic tutorials and walkthroughs on various
bioinformatics subjects:

## Dev

First clone the repository and install mkdocs and the theme using pipenv

```bash
git clone https://github.com/HadrienG/tutorials.git
cd tutorials
pipenv install
```

For a live preview in your browser do

```bash
pipenv run dev
```

## Deploy

The following command will build and push your website to a `gh-pages` branch.
Only do this if you want your own version of the website! If you are modifying
the original, please open a pull request.

```bash
pipenv run deploy
```
