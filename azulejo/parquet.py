# -*- coding: utf-8 -*-
"""Operations with pyarrow."""
# standard library imports
import functools
import io
import sys
import timeit
from pathlib import Path

# third-party imports
import click
import pandas as pd

# first-party imports
from loguru import logger

# module imports
from . import cli
from . import click_loguru

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


@cli.command()
@click_loguru.init_logger(logfile=False)
@click.option(
    "--column_names",
    default=False,
    is_flag=True,
    show_default=False,
    help="Print the names of columns.",
)
@click.option(
    "--index_name",
    default=False,
    is_flag=True,
    show_default=False,
    help="Print the name of the index.",
)
@click.option(
    "--no_index",
    default=False,
    is_flag=True,
    show_default=False,
    help="Do not write index column.",
)
@click.option(
    "--no_header",
    "-h",
    default=False,
    is_flag=True,
    show_default=False,
    help="Do not write the header.",
)
@click.option(
    "--column",
    "-c",
    default=None,
    multiple=True,
    show_default=True,
    help="Write only the named column.",
)
@click.option(
    "--index_val",
    "-i",
    default=None,
    multiple=True,
    show_default=True,
    help="Write only the row with this index value.",
)
@click.option(
    "--first_n",
    default=0,
    show_default=False,
    help="Write only the first N rows.",
)
@click.option(
    "--last_n",
    default=0,
    show_default=False,
    help="Write only the last N rows.",
)
@click.argument("parquetfile", type=click.Path(exists=True))
@click.argument("tsvfile", nargs=-1)
def parquet_to_tsv(
    parquetfile,
    tsvfile,
    column_names,
    index_name,
    column,
    no_index,
    first_n,
    last_n,
    index_val,
    no_header,
):
    """Reads parquet file, writes tsv."""
    parquetpath = Path(parquetfile)
    write_index = not no_index
    write_header = not no_header
    if tsvfile == ():
        tsvfile = sys.stdout
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
    elif first_n > 0:
        df = df[:first_n]
    elif last_n > 0:
        df = df[last_n:]
    if column_names:
        for column in df.columns:
            print(column)
    elif index_name:
        print(df.index.name)
    else:
        if column != ():
            columns = []
            for col in column:
                if col not in df.columns:
                    logger.error(
                        f'name "{column}" not found in list of columns '
                    )
                    sys.exit(1)
                columns.append(col)
            df.to_csv(
                tsvfile,
                sep="\t",
                columns=columns,
                index=write_index,
                header=write_header,
            )
        else:
            df.to_csv(
                tsvfile, sep="\t", index=write_index, header=write_header
            )


@cli.command()
@click_loguru.init_logger(logfile=False)
@click.argument("compression")
@click.argument("infile", type=click.Path(exists=True))
def change_compression(infile, compression):
    """Change the compression type of a Parquet file."""
    pd.read_parquet(infile).to_parquet(
        infile, compression=check_compression(compression)
    )


@cli.command()
@click_loguru.init_logger(logfile=False)
@click.option(
    "--compression",
    "-c",
    default=None,
    show_default=True,
    help="Specify the compression to be used.",
)
@click.argument("tsvfile", type=click.Path(exists=True))
@click.argument("parquetfile", type=click.Path())
def tsv_to_parquet(tsvfile, parquetfile, compression):
    """Write a TSV file as a parquet file using specified compression."""
    pd.read_csv(tsvfile, sep="\t", index_col=0).to_parquet(
        parquetfile, compression=check_compression(compression)
    )


def parse_parquet(filebuf):
    """Use pandas to parse a Parquet file."""
    pd.read_parquet(filebuf)


@cli.command()
@click_loguru.init_logger(logfile=False)
@click.argument("infile", type=click.Path(exists=True))
def time_parquet_parsing(infile):
    """Calculate the time it takes to parse a buffered file"""
    filebuf = io.BytesIO(Path(infile).open("rb").read())
    t = timeit.Timer(functools.partial(parse_parquet, filebuf))
    time_s = t.timeit(N_TIMINGS)
    logger.info(f"File parsed in memory in {time_s:.2f} s")
