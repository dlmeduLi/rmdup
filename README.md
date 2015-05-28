# rmdup

`rmdup` is a python script which can be used to remove dumplicated alignments in the mapping file. `rmdup` is meant to be the upgraded version of [unique-sam](https://github.com/dlmeduLi/unique-sam). For the underlying library, it use [pysam](http://pysam.readthedocs.org/en/latest/) to sort and parse the sam/bam file.

## Install

`rmdup` is denpended on [pysam](http://pysam.readthedocs.org/en/latest/) library, you should install `pysam` library first:
```shell
pip install pysam
```
`rmdup.py` is a python script, simply run the script in the command line:

```python
rmdup.py [options] in.bam out.sam 
```

## Usage

The input file of `rmdup.py` must be in `.bam` format. 

## Unique Strategy

### Normal Mode (normal mode)

### Strict Mode (strict mode `-t`)

## Log File Specification



## Suggestions

If you have any problems or suggestions, please contact [dlmeduli@163.com](mailto:dlmeduli@163.com).

