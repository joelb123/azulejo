# -*- coding: utf-8 -*-
"""King Phylip came over for great spaghetti."""

PREFIXES = {
    "giga": -5,
    "mega": -4,
    "magn": -4,
    "hyper": -3,
    "mir": -3,
    "super": -2,
    "epi": -1,
    "": 0,
    "sub": 1,
    "infra": 2,
    "parv": 3,
    "micro": 3,
    "nan": 4,
    "nano": 4,
    "min": 4,
}
RANKS = {
    "domain": 10,
    "realm": 10,
    "kingdom": 20,
    "phylum": 30,
    "class": 40,
    "legion": 50,
    "cohort": 60,
    "order": 70,
    "family": 80,
    "tribe": 90,
    "genus": 100,
    "section": 110,
    "series": 120,
    "species": 130,
    "group": 140,
    "varietas": 150,
    "variety": 150,
    "morph": 150,
    "forma": 160,
    "serovar": 160,
    "cultivar": 160,
    "strain": 170,
    "landrace": 170,
    "grex": 170,
}


def rankname_to_number(rankname):
    """Convert a rank name to a number."""
    for rank in RANKS:
        if rankname.endswith(rank):
            rankval = RANKS[rank]
            rankval += prefix_to_number(rankname[: -len(rank)])
            return rankval
    raise ValueError(f'rank "{rankname}" not found in list of ranks')


def prefix_to_number(prefix):
    """Return the number of the prefix."""
    if prefix in PREFIXES:
        return PREFIXES[prefix]
    raise ValueError(f'prefix "{prefix}" not found in list of prefixes')


def print_taxonomic_ranks():
    """Prints the allowed values of taxonomic ranks."""
    print("Allowed taxonomic rank names are, in order of depth:")
    for rank in RANKS:
        print(f"\t{rank}\t{RANKS[rank]}")
    print(
        "Each taxonomic rank may be modified by one of the following prefixes:"
    )
    for prefix in PREFIXES:
        if prefix == "":
            prefixname = "[None]"
        else:
            prefixname = prefix
        print(f"\t{prefixname}\t{PREFIXES[prefix]}")
