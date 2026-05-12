# RD Packages Website

This repository hosts the public website for RD Packages:

https://rdpackages.github.io

The site is a lean entry point for regression discontinuity software, references, replication links, contributors, and funding acknowledgments. Package documentation lives in the individual package repositories.

## Packages

- [rdrobust](https://github.com/rdpackages/rdrobust): estimation, inference, and graphical procedures using local polynomial and partitioning regression methods.
- [rdhte](https://github.com/rdpackages/rdhte): estimation and inference for heterogeneous treatment effects.
- [rdlocrand](https://github.com/rdpackages/rdlocrand): finite-sample and large-sample estimation and inference using local randomization and related methods.
- [rddensity](https://github.com/rdpackages/rddensity): manipulation testing using local polynomial density methods.
- [rdpower](https://github.com/rdpackages/rdpower): power, sample size, and minimum detectable effects calculations using robust bias-corrected local polynomial inference.
- [rd2d](https://github.com/rdpackages/rd2d): estimation and inference for boundary discontinuity designs.
- [rdmulti](https://github.com/rdpackages/rdmulti): estimation, inference, RD plots, and extrapolation with multiple cutoffs and multiple scores.

The legacy package URLs under this site, such as `/rdrobust/`, are kept as redirect stubs to the corresponding GitHub repositories. Mirror-style aliases under `/rdpackages/<package>/` are also generated as redirects.

## Repository Structure

- `_config.yml`: Jekyll and GitHub Pages settings.
- `.editorconfig` and `.gitattributes`: editor and line-ending conventions.
- `Gemfile`: local GitHub Pages/Jekyll dependencies.
- `_includes/`: shared HTML fragments for metadata.
- `_layouts/`: page templates.
- `index.md`: homepage content.
- `<package>/index.md`: package redirect stubs.
- `rdpackages/<package>/index.md`: mirror-style package redirect aliases.
- `replication/index.md`: replication resources.
- `references/`: PDF references linked from the website.
- `public/`: CSS and static assets.

## Editing

Most content changes happen in Markdown files. Each rendered page needs YAML front matter at the top, for example:

```yaml
---
layout: page
title: Page title
permalink: /page-slug/
---
```

Use relative internal links when possible, such as `/rdrobust/` or `/references/file.pdf`, so the site remains portable across local previews and GitHub Pages.

## Local Preview

This is a Jekyll site. On a machine with Ruby and Bundler installed:

```bash
bundle install
bundle exec jekyll serve
```

Then open:

```text
http://127.0.0.1:4000
```

This checkout includes a `Gemfile` and `Gemfile.lock` for GitHub Pages/Jekyll dependencies.

Note: local Windows builds were tested with Ruby 3.3.10. Ruby 3.4 exposed compatibility issues in the pinned `github-pages`/Jekyll 3.10 stack, so use Ruby 3.3.x for local GitHub Pages previews until the dependency stack is upgraded.

## Deployment

The repository is intended for GitHub Pages at `rdpackages.github.io`. The current source branch is `main`.

## Queries and Requests

Please email [rdpackages@googlegroups.com](mailto:rdpackages@googlegroups.com).
