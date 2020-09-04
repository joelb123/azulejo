# -*- coding: utf-8 -*-
"""Operations with pyarrow."""
# standard library imports
import functools
import io
import sys
import timeit
from pathlib import Path

# third-party imports
import pandas as pd

# first-party imports
from loguru import logger

# global constants
COMPRESSION_TYPES = [
    "SNAPPY",
    "GZIP",
    "BROTLI",
    "NONE",
    "LZ4",
    "LZO",
    "ZSTD",
    None,
]
N_TIMINGS = 10


def check_compression(compression):
    """Check to see if the specified compression is recognized."""
    if compression not in COMPRESSION_TYPES:
        logger.error(f"Unrecognized compression type {compression}")
        logger.error(f"Valid compressions are: {COMPRESSION_TYPES}")
        sys.exit(1)


def parquet_to_tsv(
    parquetfile,
    tsvfile,
    columns,
    index_name,
    col,
    no_index,
    head,
    tail,
    index_val,
    no_header,
    pretty,
    max_rows,
    max_cols,
    writefile,
):
    """Reads parquet file, writes tsv."""
    parquetpath = Path(parquetfile)
    write_index = not no_index
    write_header = not no_header
    if max_rows is not None:
        max_rows = int(max_rows)
    if max_cols is not None:
        max_cols = int(max_cols)
    if tsvfile == ():
        tsvfile = sys.stdout
    if writefile:
        tsvfile = parquetpath.stem + ".tsv"
    if not parquetpath.exists():
        logger.error(f'parquet file "{parquetpath} does not exist.')
        sys.exit(1)
    df = pd.read_parquet(parquetpath)
    if index_val != ():
        index_values = []
        for val in index_val:
            if val not in df.index:
                logger.error(f'name "{val}" not found in index')
                sys.exit(1)
            index_values.append(val)
        df = df[df.index.isin(index_values)]
    elif head > 0:
        df = df[:head]
    elif tail > 0:
        df = df[tail:]
    if col != ():
        col_list = []
        for c in col:
            if c not in df.columns:
                logger.error(f'name "{c}" not found in list of columns ')
                sys.exit(1)
            col_list.append(c)
        print(col_list)
        df = df[col_list]
    if pretty:
        with pd.option_context(
            "display.max_rows",
            max_rows,
            "display.max_columns",
            max_cols,
            "display.float_format",
            "{:,.2f}%".format,
        ):
            print(df)
    elif columns:
        print(df.dtypes)
    elif index_name:
        print(df.index.name)
    else:
        df.to_csv(tsvfile, sep="\t", index=write_index, header=write_header)


def change_compression(infile, compression):
    """Change the compression type of a Parquet file."""
    pd.read_parquet(infile).to_parquet(
        infile, compression=check_compression(compression)
    )


def tsv_to_parquet(tsvfile, parquetfile, compression):
    """Write a TSV file as a parquet file using specified compression."""
    pd.read_csv(tsvfile, sep="\t", index_col=0).to_parquet(
        parquetfile, compression=check_compression(compression)
    )


def parse_parquet(filebuf):
    """Use pandas to parse a Parquet file."""
    pd.read_parquet(filebuf)


def time_parquet_parsing(infile):
    """Calculate the time it takes to parse a buffered file"""
    filebuf = io.BytesIO(Path(infile).open("rb").read())
    t = timeit.Timer(functools.partial(parse_parquet, filebuf))
    time_s = t.timeit(N_TIMINGS)
    logger.info(f"File parsed in memory in {time_s:.2f} s")
