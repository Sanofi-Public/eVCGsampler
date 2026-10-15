# eVCGsampler - Help

The app draws a **Virtual Control Group (VCG)** from a historical control pool so
that it is balanced to the treated animals of your study, using energy-distance
based covariate balancing (R package **eVCGsampler**).

Whole, unweighted animals are selected. Nothing is imputed and nothing is
re-weighted, so the VCG can be analysed like any other control group.

> Balancing is a **final refinement step**. It reduces biological heterogeneity,
> but it does not replace careful pre-filtering of the historical control data,
> and it does not address technical heterogeneity caused by different collection
> methods, processing pipelines or experimental conditions.

Your data is processed temporarily inside the session. Nothing is uploaded,
stored or shared.

---

## 1. Data

Upload the historical control **POOL**. The treated animals (**TG**) can come
either way:

* **One file** - it already contains a variable with exactly two values, one for
  the pool animals and one for the treated animals.
* **Two files** - POOL first, TG second. The two are stacked and the group
  variable `POOL_vs_TG` is created for you.

**Supported formats:** `.csv`, `.txt`, `.xlsx`, `.xls`. For delimited files the
separator (`,` `;` tab) and the decimal mark are guessed from the first lines.

### Trying the app out

**Load the example data set** fills both slots with simulated data - 12
historical control studies (240 animals) plus 40 treated animals, with two
studies deliberately shifted on the outcome variables so every step has
something to show. An *example data* badge stays visible for the whole session
so a demo run cannot be mistaken for a real one. Uploading a POOL file switches
back to your own data.

### Combining check

When you upload two files, only variables that appear in **both** of them survive
the merge - that is how `combine_data()` works. The *Combining check* tab lists
every column and what happened to it:

| Status | Meaning |
|---|---|
| combined | present in both files, available for balancing |
| dropped | present in only one file, not available |
| type mismatch | present in both, but one side is text where the other is numeric |

A type mismatch is almost always a data problem worth fixing - typically text in
a numeric column (`"<0.1"`, `"n.d."`, a stray unit) or a decimal comma that was
not recognised.

### Variable overview

*Usable as* tells you which role each column can play:

* **group** - exactly two distinct non-missing values.
* **stratum** - at most 8 distinct values (eVCGsampler allows no more than 8
  strata in total).
* **covariate** - numeric, logical, or a category with at most 8 levels.

Free-text columns and per-animal identifiers cannot be *balanced*, so they are
not offered as covariates or strata. An identifier is still useful: pick it as
the **Animal ID** on the Setup step and it is carried through to the results.

---

## 2. Setup

Pick the animal ID, the group variable, the covariates to balance, optional
strata, and the VCG size. The selection is turned into the eVCGsampler formula
shown on the right:

```
POOL_vs_TG ~ weight + age + sex_score | study_phase
```

Covariates are joined with `+`; anything after `|` is a stratum. Balancing then
runs **inside each stratum** separately.

### Animal ID

The column that identifies an individual animal - without it you get a balanced
VCG but no way to tell *which* animals it is made of.

The first suitable column is picked automatically. To qualify, a column must
have no missing values and must not repeat a value **within** a group: POOL and
TG usually come from separate files and each may number its animals from 1, so
an ID that appears once in the POOL and once in the TG is still fine. Only POOL
animals can end up in the VCG, so that is enough to identify them.

The ID takes no part in the balancing - it is carried along untouched and put in
the **first column** of the Sampled data table, both exports and the report's
*Selected animals* section. Choose *(none)* if your data has no identifier; the
sampling works exactly the same, you just cannot trace the animals afterwards.

### VCG size

One number balances the whole data set, or is used in every stratum. Give one
number per stratum (comma-separated, e.g. `7, 5`) to use different sizes - the
order follows the strata shown in the sample size table.

### Covariate importance weights

Optional, one positive weight per covariate (e.g. `2, 1, 1`). A higher weight
pulls that covariate closer to the TG at the expense of the others. Fewer than
two weights are ignored, because scaling every covariate by the same factor
changes nothing.

### Sample sizes

The counts are taken **after** dropping rows with a missing value in any selected
variable, which is exactly the data the sampler and the tests will use. Watch the
status column:

| Status | What to do |
|---|---|
| ok | continue |
| size >= 75% of POOL | the sampler runs, but has little room to balance - widen the pool if you can |
| very small POOL | fewer than 10 pool animals in a stratum; the result will be unstable |
| size > POOL | not enough animals to draw from - reduce the size |
| no TG animals | a stratum has nothing to balance towards |
| size missing | give one size, or exactly one per stratum |

---

## 3. Pre-tests

### Energy test - POOL vs TG

Answers *"is there an imbalance worth correcting?"*.

The **energy distance** measures how far two multivariate distributions are
apart; 0 means identical. The permutation test reshuffles the group labels many
times to see how large that distance gets by chance alone, and compares your
observed distance against it. The test is one-sided.

* **small p-value** - the POOL and the TG really do differ. This is the situation
  VCG sampling is for.
* **large p-value** - no detectable imbalance to begin with; a plain random draw
  from the pool would do.

More permutations give a more stable p-value and take longer. 1000 is a good
default, 2000+ for a number you want to report.

With strata the test is run inside each stratum and reported per stratum; the
permutation histogram is only available without stratification.

### Outlier studies

Answers *"does any source study stand out from the rest?"*.

Historical control data pooled from several studies sometimes contains a study
whose **outcome** distribution sits far from all others. This test compares every
study against every other one by energy distance and flags those beyond the
chosen percentile of the permuted distances.

Use it on **outcome or technical variables**, not on the covariates you are
balancing - those are handled by the sampling itself. By default only the POOL
animals are compared.

A flagged study is not automatically wrong, but it is worth asking whether it
belongs in the pool at all. Excluding it is a decision for the Data step, not
something the sampler can fix.

---

## 4. Sampling

**Draw the VCG** runs the sampler.

* **Best subset** (default) - the selection that balances best.
* **Random subset sampling** - animals are drawn with a probability derived from
  their energy weights. Use this only when you genuinely need several different
  VCGs, e.g. for a proof-of-concept study. It adds randomness at the cost of
  balance quality, and it is *not* meant for generating several VCGs and then
  picking the nicest one.
* **Number of VCGs** (random subset sampling only) - more than one draws that
  many independent random VCGs in one go (`multiSampler()`). The Balance tab then
  shows how often each POOL animal was drawn (the overlap plot), the median (MAD)
  per VCG and a post-hoc test for every VCG; the Covariates tab has a selector
  to look at each VCG. The downloads hold every animal drawn at least once, with
  one `VCG_1` ... `VCG_N` column per VCG and a `Times_selected` count.

The **random seed** makes a run reproducible: both the random-subset option and
the permutation ellipse in the balance plot draw random numbers.

### Balance target plot

Each arrow starts at the POOL-vs-TG difference and points at the VCG-vs-TG
difference, in median (x axis) and MAD (y axis). The purple ellipse is the 95%
range of differences seen under random permutation, so an arrow tip **inside**
the ellipse means the remaining difference is ordinary noise.

The ellipse is interpolated between its two axes rather than computed in full -
it is a visual approximation, not a formal test.

### Post-hoc energy test - TG vs VCG

The same permutation test as before, now on the variable the sampler created
(`<group>_balanced`, TG vs the selected VCG). Here a **large** p-value is the
goal: it means the VCG is no longer distinguishable from the TG.

If it is still small, the pool simply does not contain animals close enough to
your TG. Options: widen the pool, balance fewer covariates, drop a covariate that
cannot be matched, or use weights to prioritise the covariates that matter most.

### Median (MAD) per group

`POOL` is the whole pool, `VCG` the selected subset, `TG` the treated animals.
The `VCG` row should sit closer to `TG` than the `POOL` row does.

### Sampled data

Only the **selected** animals are shown - the rows where `VCG` is `1`. The two
download buttons export exactly these rows. Columns, with the animal ID first:

* your ID column - which animals to pull from the pool.
* the strata and the balanced covariates.
* `VCG` - `1` throughout, since the table is filtered to the selection.
* `e_weights` - the energy weight this animal received.
* `<group>_balanced` - the group variable used for the post-hoc test.

Need the animals that were *not* picked as well - to check the pool, or to see
which animals came close? The Report step exports the complete data set.

---

## 5. Report

Produces one self-contained HTML file documenting the whole run: the data checks,
the settings, the pre-tests, the balance result and the selected animals, plus the
seed and the package versions needed to reproduce it.

Steps you skipped are simply left out.

### The two exports

| Where | What you get |
|---|---|
| Sampling step | the **selected** animals only (`VCG` = `1`) - the VCG itself |
| Report step | the **complete** data set - selected (`VCG` = `1`), unselected pool (`VCG` = `0`) and TG (`VCG` empty), plus the median (MAD) table and the sample sizes |

Both put the animal ID in the first column.

---

## Troubleshooting

**"This data set has no variable with exactly two distinct values."**
The sampler needs a POOL / TG group variable. Upload the treated animals as a
second file, or add a two-level group column.

**"No overlapping covariates were found."**
The two uploaded files share no column names. The covariates must be spelled
identically in both.

**A covariate is missing from the selection list.**
It is either free text, has only one distinct value, or is a category with more
than 8 levels.

**Too many strata.**
eVCGsampler allows at most 8 strata in total - that is the number of *combinations*
of the stratum variables, not the number of variables.

**Sampling is slow.**
Permutation tests dominate the runtime. Reduce the number of permutations while
you explore, and raise it for the final run.

---

## References

Schulz, A. *eVCGsampler: VCG Sampling using Energy-Based Covariate Balancing*.
R package, CRAN.
<https://cran.r-project.org/web/packages/eVCGsampler/index.html>

Székely, G. J., & Rizzo, M. L. (2016). Energy distance. *Wiley Interdisciplinary
Reviews: Computational Statistics*, 8(1), 27-38.
