# -*- coding: utf-8 -*-
"""Protein sequence checking and sanitization."""

# standard library imports
import zlib

# third-party imports
from Bio.Data import IUPACData
from loguru import logger

# global constants
AMBIGUOUS_CHARACTER = "X"
ALPHABET = IUPACData.protein_letters + AMBIGUOUS_CHARACTER + "-"


class Sanitizer(object):

    """Count and clean up problems with protein sequence.

    Problems recognized are:
          alphabet:  if not in IUPAC set, changed to 'X'
            dashes:    optional, removed if remove_dashes=True
         ambiguous:

    """

    def __init__(self, remove_dashes=False):
        """Initialize counters."""
        self.remove_dashes = remove_dashes
        self.seqs_sanitized = 0
        self.seqs_out = 0
        self.chars_in = 0
        self.chars_removed = 0
        self.chars_fixed = 0
        self.endchars_removed = 0
        self.chars_out = 0
        self.ambiguous = 0

    def char_remover(self, s, character):
        """Remove positions with a given character.

        :param s: mutable sequence
        :return: sequence with characters removed
        """
        removals = [i for i, j in enumerate(s) if j == character]
        self.chars_removed += len(removals)
        [s.pop(pos - k) for k, pos in enumerate(removals)]  # pylint: disable=expression-not-assigned
        return s

    def fix_alphabet(self, s):
        """Replace everything out of alphabet with AMBIGUOUS_CHARACTER.

        :param s: mutable sequence, upper-cased
        :return: fixed sequence
        """
        fix_positions = [pos for pos, char in enumerate(s) if char not in ALPHABET]
        self.chars_fixed += len(fix_positions)
        [s.__setitem__(pos, AMBIGUOUS_CHARACTER) for pos in fix_positions]  # pylint: disable=expression-not-assigned
        return s

    def remove_char_on_ends(self, s, character):
        """Remove leading/trailing characters..

        :param s: mutable sequence
        :return: sequence with characters removed from ends
        """
        in_len = len(s)
        while s[-1] == character:
            s.pop()
        while s[0] == character:
            s.pop(0)
        self.endchars_removed += in_len - len(s)
        return s

    def sanitize(self, s):
        """Sanitize potential problems with sequence.

        Remove dashes, change non-IUPAC characters to
        ambiguous, and remove ambiguous characters on ends.
        :param s: mutable sequence
        :return: sanitized sequence
        """
        self.seqs_sanitized += 1
        self.chars_in += len(s)
        if not len(s):
            raise ValueError("zero-length sequence")
        if self.remove_dashes:
            s = self.char_remover(s, "-")
        if not len(s):
            raise ValueError("zero-length sequence after dashes removed")
        s = self.fix_alphabet(s)
        s = self.remove_char_on_ends(s, AMBIGUOUS_CHARACTER)
        if not len(s):
            raise ValueError("zero-length sequence after ends trimmed")
        self.chars_out += len(s)
        self.seqs_out += 1
        return s

    def count_ambiguous(self, s):
        """Count ambiguous residues.

        :param s: sequence
        :return: Number of ambiguous residues
        """
        ambig = sum([i == AMBIGUOUS_CHARACTER for i in s])
        self.ambiguous += ambig
        return ambig

    def file_stats(self):
        """Return a dictionary of file stats."""
        return {
            "seqs_in": self.seqs_sanitized,
            "seqs_out": self.seqs_out,
            "residues_in": self.chars_in,
            "residues_out": self.chars_out,
            "dashes_removed": self.chars_removed,
            "chars_fixed": self.chars_fixed,
            "ends_trimmed": self.endchars_removed,
            "ambiguous_residues": self.ambiguous,
        }


class DuplicateSequenceIndex(object):

    """Count duplicated sequences."""

    def __init__(self, concat_names=False):
        self.match_index = 0
        self.hash_set = set()
        self.duplicates = {}
        self.match_count = {}

    def exact(self, s):
        "Test and count if exact duplicate."
        seq_hash = zlib.adler32(bytearray(str(s), "utf-8"))
        if seq_hash not in self.hash_set:
            self.hash_set.add(seq_hash)
            return ""
        if seq_hash not in self.duplicates:
            self.duplicates[seq_hash] = self.match_index
            self.match_count[self.match_index] = 1
            self.match_index += 1
        else:
            self.match_count[self.duplicates[seq_hash]] += 1
        return str(self.duplicates[seq_hash])

    def counts(self, index):
        """Return the number of counts for a match index."""
        self.match_count[int(index)]
