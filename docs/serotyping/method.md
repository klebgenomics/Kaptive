---
title: Method
author: Tom Stanton
comments: true
tags: [serotyping]
icon: lucide/syringe
categories:
  - Serotyping
---

## :lucide-microscope: Serotyping Overview

The [`Serotyper`](../reference/serotyping/core.html#kaptive.serotyping.core.Serotyper) executes a highly-optimized, multi-phase pipeline to analyze a
genome assembly and produce a [SerotypingResult](../reference/serotyping/models.html#kaptive.serotyping.models.SerotypingResult). 


## :lucide-calculator: Locus Scoring Algorithm

This provides a detailed overview of the mathematical scoring algorithm that Kaptive uses to determine the best-matching locus for a given genome assembly. The algorithm is designed to be highly accurate and computationally efficient, utilizing dynamic programming alignment scores and strict completeness penalties to resolve locus identities.

### 1. Initial Gene Mapping

When Kaptive analyses a genome, it maps all genes from the reference database against the assembly using a highly optimized, minimap2-based banded dynamic programming aligner. For each alignment, Kaptive calculates two key metrics:

* **Alignment Score ($S_{DP}$):** The raw Smith-Waterman/Gotoh dynamic programming score with affine gap penalties, which accurately reflects the biological likelihood of the alignment.
* **Gene Coverage ($Cov_{gene}$):** The proportion of the reference gene length that was successfully aligned.

$$Cov_{gene} = \frac{\text{Alignment Length}}{\text{Expected Gene Length}}$$

### 2. Gene Selection

Since a single gene might map to multiple locations (e.g., duplicated elements or highly conserved core genes), Kaptive must determine the single "best" hit for each gene.

1. **Coverage Filter:** Any alignment covering less than 20% of the expected gene length is immediately discarded to prune noisy, spurious hits.
2. **Primary Hit Selection:** The remaining valid alignments for a gene are ranked by their **Alignment Score ($S_{DP}$)**, with **Gene Coverage ($Cov_{gene}$)** used as a tie-breaker. The highest-ranking alignment is selected as the primary hit.

!!! tip
    Relying on robust DP scoring for primary hit selection inherently prioritizes identical, full-length gene variants over distant homologs or fragments, ensuring the highest quality matches are carried forward into the locus scoring phase.

### 3. Accumulating Locus Scores

Once the primary hit for each gene is identified, Kaptive evaluates every possible locus in the database independently. 

It calculates a **Core Score** ($S_{core}$) for each locus by simply summing the Gene Coverage of the primary hits for all *expected* genes belonging to that locus. (Note: genes designated as optional or "extra" in the database are excluded from this core summation).

$$S_{core} = \sum_{gene \in \text{expected}} Cov_{gene}$$

!!! note
    In earlier versions of Kaptive, an *Inverse Document Frequency (IDF)* weighting system was used to penalize universally conserved "core" genes. The current algorithm replaces IDF weighting entirely; the combination of strict DP-based hit selection and the exponential completeness penalty (described below) naturally and robustly separates true locus matches from false positives driven by core genes.

### 4. Locus Completeness Penalty

A high Core Score alone is not enough to declare a match. A locus might accumulate a high score if a genome happens to contain highly-conserved core genes elsewhere in the assembly, but is missing the rare, specific genes that uniquely define that locus.

To strongly penalize loci with missing genes, Kaptive calculates a **Completeness** factor ($C$). Crucially, this completeness is evaluated strictly within the **reconstructed locus boundaries**—an expected gene is only counted as "found" if it physically overlaps the assembled locus region.

$$C = \frac{\text{Count}_{found\_inside}}{\text{Count}_{expected}}$$

### 5. Final Total Score and Selection

The final score used to rank and select the best locus is the product of the Core Score and the Completeness factor **cubed**.

$$\text{Score}_{total} = S_{core} \times C^3$$

!!! tip
    By multiplying the Core Score by $C^3$, Kaptive applies an **exponential penalty** to loci missing expected genes. For example, a locus that is only 50% complete will have its score scaled down by a factor of 0.125 ($0.5^3$). This ensures that complete or near-complete loci will consistently out-score fragmented false positives.

#### Summary of Selection

After calculating $\text{Score}_{total}$ for every possible locus in the database, Kaptive ranks them in descending order. The locus with the highest $\text{Score}_{total}$ is chosen as the **Best match locus**. Following this algorithmic selection, a secondary "Phenotype Evaluation" phase may be triggered (e.g., checking specific alleles or gene states) to refine the final serotype prediction.

## :lucide-flask-conical: Phenotype Evaluation Logic

While predicting the "best match locus" is heavily reliant on quantitative metrics (Core Score and Completeness), the final biological prediction of what the genome actually *does* or *looks like* requires qualitative logic. We refer to this as the **Phenotype Evaluation Phase**.

In bacterial serotyping, two genomes may share the exact same underlying locus backbone (e.g., KL1), but specific mutations, gene disruptions, or mobile genetic element insertions can alter the final expressed serotype or phenotype (e.g., adding an acetyl group, or turning off capsule production entirely).

Kaptive uses a highly-optimized Boolean rule engine to apply these nuanced biological rules.

### 1. Determining Gene States
Before any phenotype rules can be evaluated, Kaptive must assess the "health" of every gene found in the locus. 

Each gene is translated into its corresponding amino acid sequence (truncating at the first stop codon to prevent downstream frameshifts from skewing the identity score) and aligned to the database reference. Based on this alignment and coverage, it is assigned a **Gene State**:

* **NORMAL**: The gene is fully intact and its sequence identity exceeds the database's minimum identity threshold.
* **PARTIAL**: The gene was cut off by the edge of a fragmented assembly contig. Because we cannot know if the missing portion is intact, we assume it is functional for the sake of caution.
* **TRUNCATED**: The gene's translated protein sequence covers less than 90% of the expected reference protein length, indicating a premature stop codon, frameshift, or massive internal deletion.
* **NOVEL**: The gene is fully present, but its sequence identity falls below the database threshold, suggesting it may be a novel variant with potentially different function. Spurious homologous genes found *outside* the locus boundaries that fall below this identity threshold are completely ignored.

### 2. Defining "Active" Gene Clusters
In Kaptive's logic, a gene cluster (a group of homologous genes) is considered **Active** in the genome only if a corresponding gene is found in a `NORMAL` or `PARTIAL` state.

$$\text{Active}(g) = \begin{cases} \text{True} & \text{if state}(g) \in \{\text{NORMAL}, \text{PARTIAL}\} \\ \text{False} & \text{if state}(g) \in \{\text{TRUNCATED}, \text{NOVEL}, \text{MISSING}\} \end{cases}$$

### 3. The Boolean Rule Engine
Phenotype rules are defined in the database. Each rule contains three critical conditions that must be met simultaneously for the rule to trigger.

1. **Locus Match**: The rule must be explicitly valid for the currently predicted Best Match Locus.
2. **Presence Conditions**: A specific set of gene clusters must be *Active*.
3. **Absence Conditions**: A specific set of gene clusters must be *Inactive* (missing, broken, or novel).

Mathematically, a rule is valid if and only if:

$$\text{Valid}(R) = (\text{Locus} \in R_{loci}) \wedge \left( \forall g \in R_{required}, \text{Active}(g) \right) \wedge \left( \forall g \in R_{absent}, \neg\text{Active}(g) \right)$$

*By using highly-optimized bitwise operations on boolean masks, Kaptive evaluates these rules across thousands of genomes almost instantaneously.*

### 4. Applying the Rules
Every locus has a default "Base Phenotype" (for example, the locus `KL104` might have a default phenotype of `K104`). 

When one or more phenotype rules evaluate to True, they modify the Base Phenotype in two distinct ways depending on the type of the rule:

#### A. Replacement Rules
Replacement rules indicate a fundamental shift in the serotype. If a replacement rule is triggered, the Base Phenotype is completely overwritten by the new phenotype string.
* *Example:* If the genome is `KL30`, but a critical *wzy* gene is mutated/inactive, a rule might completely replace the phenotype with `Untypeable` or `K30-mutant`.
* *Conflict Resolution:* If multiple replacement rules trigger simultaneously, the rule with the highest predefined **Priority** wins.

#### B. Suffix Rules
Suffix rules indicate an additive change or a minor modification. They are appended to the end of the Base Phenotype.
* *Example:* If a mobile element inserts an extra acetyltransferase gene, a suffix rule might append `+Ac` to the phenotype, resulting in `K104+Ac`.
* *Conflict Resolution:* If multiple suffix rules trigger, they are all appended to the phenotype string in descending order of their **Priority**.

#### Final Output
Once all valid replacements and suffixes are processed, the final string is returned as the predicted **Phenotype** of the genome assembly.
