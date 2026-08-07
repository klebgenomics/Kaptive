---
title: Release Notes
author: Tom Stanton
comments: true
tags: [markdown, documentation, web]
icon: lucide/rocket
categories:
  - Development
---

# 🐋 Kaptive is containerized!
*Published on 2026-08-03*

## What's Changed
* Merge pull request #2 from klebgenomics/master by @tomdstanton in https://github.com/klebgenomics/Kaptive/pull/69


**Full Changelog**: https://github.com/klebgenomics/Kaptive/compare/v3.3.0...v3.3.1

---

# Database Decentralisation 📦
*Published on 2026-07-31*

## What's new?

### Decentralised databases 📦

- Kaptive v3.3.0 fully decouples databases from the tool.
- Fully versioned databases.
- Only install the databases you need.
- Easier for curators to maintain their databases.

### A new core 🍎
 
- Kaptive v3.3.0 is built around a new ultra fast API-first core.
- Allows for easier integration of new modules like plotting and annotation.
- Objects are structured with a data-oriented-design approach leveraging numpy for SoA (structure of arrays) batching.
- Designed for vectorised computation.
- Easier integration of machine learning in the future.
- Results are stored in frozen, slotted dataclasses - designed for efficiency, thread-safety and integration with web APIs. 
- Extremely fast with numba!

### Introduction of [rammappy](https://github.com/tomdstanton/rammappy) 🚀

- Kaptive is now 100% pip installable.
- No more background sub-processes, all computation happens directly in Kaptive.
- Safe for Web APIs.

**Full Changelog**: https://github.com/klebgenomics/Kaptive/compare/v3.2.2...v3.3.0

---

# v3.2.2
*Published on 2026-07-16*

## What's Changed
* Migrate to mkdocs. by @tomdstanton in https://github.com/klebgenomics/Kaptive/pull/54
* Bump urllib3 from 2.6.3 to 2.7.0 by @dependabot[bot] in https://github.com/klebgenomics/Kaptive/pull/57
* Bump tj-actions/changed-files from 42 to 46 in /.github/workflows by @dependabot[bot] in https://github.com/klebgenomics/Kaptive/pull/60
* Bump idna from 3.11 to 3.15 by @dependabot[bot] in https://github.com/klebgenomics/Kaptive/pull/59

## New Contributors
* @dependabot[bot] made their first contribution in https://github.com/klebgenomics/Kaptive/pull/57

**Full Changelog**: https://github.com/klebgenomics/Kaptive/compare/v3.2.1...v3.2.2

---

# v3.2.1
*Published on 2026-03-04*

## What's Changed
* Remove dna features viewer by @tomdstanton in https://github.com/klebgenomics/Kaptive/pull/52


**Full Changelog**: https://github.com/klebgenomics/Kaptive/compare/v3.2.0...v3.2.1

---

# v3.2.0
*Published on 2026-02-23*

## What's Changed
* Add missing links to Outputs.rst by @pvanheus in https://github.com/klebgenomics/Kaptive/pull/43
* Update annotations for KpSC K-locus genes by @tomdstanton in https://github.com/klebgenomics/Kaptive/pull/49
* Update KpSC K logic allele names by @tomdstanton in https://github.com/klebgenomics/Kaptive/pull/51

## New Contributors
* @pvanheus made their first contribution in https://github.com/klebgenomics/Kaptive/pull/43

**Full Changelog**: https://github.com/klebgenomics/Kaptive/compare/v3.1.0...v3.2.0

---

