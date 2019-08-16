# Lecture
This directory is for my lecture (+ What I learned outside of school) notes.

Please Look at [here!!](https://iwasakishuto.github.io/study/notes/)

## How to write
### markdown
1. Write down note as markdown format.
Note: You have to write this at the header.
ex.)
```
Title: <Lecture name> 第k回
Date: YYYY-MM-DD
Modified: YYYY-MM-DD
Slug: <Lecture name>-k
Tags: <Term>, <Lecture name>
Category: <Lecture name>
Authors: Shuto
Summary: <Summary of the notes>
FacebookImage: https://iwasakishuto.github.io/study/notes/theme/img/default_image.png
Password: <Password>
```

2. put it at `content/markdown/`

### notebook
1. put `sample.ipynb` at `content/notebook/`
2. make `content/notebook/sample.ipynb-meta` and write information above.

### both
Run `make html`
<b>WARNING: YOU SHOULD NOT RUN `make github`</b>

## Pelican theme
I use [elegant](https://github.com/Pelican-Elegant/elegant) as a pelican-theme. Special thanks to them.

<img src="https://raw.github.com/Pelican-Elegant/elegant/master/home-page-screenshot.png">
