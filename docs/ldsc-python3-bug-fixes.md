# LDSC Python 3 Compatibility Bug Fixes

**Repository:** [twillis209/ldsc](https://github.com/twillis209/ldsc) — branch `fix/annot-python3-scoping`

---

## Bug 1: List comprehension variable scoping (`--overlap-annot` crash)

**Symptom**

Running `ldsc.py` with `--overlap-annot` raised an `UnboundLocalError`:

```
UnboundLocalError: cannot access local variable 'fh' where it is not associated with a value
```

**Root cause**

In Python 2, loop variables in list comprehensions leaked into the enclosing scope. The `annot()` function in `ldscore/parse.py` relied on this:

```python
annot_suffix = ['.annot' for fh in fh_list]  # fh leaked into scope in Python 2
...
chrs = get_present_chrs(fh, num+1)            # fh undefined in Python 3
```

In Python 3, list comprehensions have their own scope, so `fh` was never defined at the point of use.

**Fix**

Reference the list directly instead of relying on the leaked variable:

```python
chrs = get_present_chrs(fh_list[0], num+1)
```

---

## Bug 2: Frq/annot SNP set mismatch (`--overlap-annot` ValueError)

**Symptom**

On some traits, `ldsc.py --overlap-annot` produced a confusing `TypeError`:

```
TypeError: cannot unpack non-iterable ValueError object
```

with the log showing `Error parsing annot matrix` / `Error parsing .annot file`.

**Root cause**

`annot_parser()` filters the annotation matrix to common SNPs (MAF 0.05–0.95) using a boolean mask derived from a separate frequency file. The two files (baselineLD v2.2 annotations and 1000G Phase 3 frq files) have different SNP sets per chromosome. The filtering code used positional boolean indexing:

```python
df_annot = df_annot[((.95 > df_frq.FRQ) & (df_frq.FRQ > 0.05)).values]
```

When `df_frq` and `df_annot` have different lengths, numpy raises `ValueError: Boolean array expected to be of length X, received Y`. LDSC's `_read_chr_split_files` catches `ValueError` and *returns* the exception object rather than raising it, causing the caller's tuple-unpack assignment to fail with `TypeError`.

An earlier attempt (`.values` on the mask) fixed the pandas index-alignment warning but introduced this length-mismatch crash for traits where the SNP counts differ between the two files.

**Fix**

Merge the annotation and frequency DataFrames on `SNP` before filtering, so differing SNP sets are handled correctly via an inner join:

```python
def annot_parser(fh, compression, frqfile_full=None, compression_frq=None):
    df_annot = read_csv(fh, header=0, compression=compression)
    if frqfile_full is not None:
        df_frq = frq_parser(frqfile_full, compression_frq)
        if 'SNP' in df_annot.columns:
            df_annot = df_annot.merge(df_frq[['SNP', 'FRQ']], on='SNP', how='inner')
            df_annot = df_annot[(.95 > df_annot['FRQ']) & (df_annot['FRQ'] > 0.05)]
            df_annot = df_annot.drop(['FRQ'], axis=1)
    df_annot = df_annot.drop(['SNP', 'CHR', 'BP', 'CM'], axis=1, errors='ignore').astype(float)
    return df_annot
```
