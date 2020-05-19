# -*- coding: utf-8 -*-
"""Protein sequence checking and sanitization."""

# standard library imports
import zlib

# first-party imports
from Bio.Data import IUPACData

# global constants
AMBIGUOUS_CHARACTER = "X"
ALPHABET = IUPACData.protein_letters + AMBIGUOUS_CHARACTER + "-"


class Sanitizer:

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

    def char_remover(self, seq, character):
        """Remove positions with a given character.

        :param seq: mutable sequence
        :return: sequence with characters removed
        """
        removals = [i for i, j in enumerate(seq) if j == character]
        self.chars_removed += len(removals)
        for k, pos in enumerate(removals):
            seq.pop(pos - k)
        return seq

    def fix_alphabet(self, seq):
        """Replace everything out of alphabet with AMBIGUOUS_CHARACTER.

        :param seq: mutable sequence, upper-cased
        :return: fixed sequence
        """
        fix_positions = [
            pos for pos, char in enumerate(seq) if char not in ALPHABET
        ]
        self.chars_fixed += len(fix_positions)
        for pos in fix_positions:
            seq.__setitem__(pos, AMBIGUOUS_CHARACTER)
        return seq

    def remove_char_on_ends(self, seq, character):
        """Remove leading/trailing characters..

        :param seq: mutable sequence
        :return: sequence with characters removed from ends
        """
        in_len = len(seq)
        while seq[-1] == character:
            seq.pop()
        while seq[0] == character:
            seq.pop(0)
        self.endchars_removed += in_len - len(seq)
        return seq

    def sanitize(self, seq):
        """Sanitize potential problems with sequence.

        Remove dashes, change non-IUPAC characters to
        ambiguous, and remove ambiguous characters on ends.
        :param seq: mutable sequence
        :return: sanitized sequence
        """
        self.seqs_sanitized += 1
        self.chars_in += len(seq)
        if len(seq) == 0:
            raise ValueError("zero-length sequence")
        if self.remove_dashes:
            seq = self.char_remover(seq, "-")
        if len(seq) == 0:
            raise ValueError("zero-length sequence after dashes removed")
        seq = self.fix_alphabet(seq)
        seq = self.remove_char_on_ends(seq, AMBIGUOUS_CHARACTER)
        if len(seq) == 0:
            raise ValueError("zero-length sequence after ends trimmed")
        self.chars_out += len(seq)
        self.seqs_out += 1
        return seq

    def count_ambiguous(self, seq):
        """Count ambiguous residues.

        :param seq: sequence
        :return: Number of ambiguous residues
        """
        ambig = sum([i == AMBIGUOUS_CHARACTER for i in seq])
        self.ambiguous += ambig
        return ambig

    def file_stats(self):
        """Return a dictionary of file stats."""
        return {
            "seqs": self.seqs_out,
            "resids": self.chars_out,
            "seqs_in": self.seqs_sanitized,
            "resids_in": self.chars_in,
            "dashes": self.chars_removed,
            "fixed": self.chars_fixed,
            "trimmed": self.endchars_removed,
            "ambig": self.ambiguous,
        }


class DuplicateSequenceIndex:

    """Count duplicated sequences."""

    def __init__(self):
        """Save stats."""
        self.match_index = 0
        self.hash_set = set()
        self.duplicates = {}
        self.match_count = {}

    def exact(self, seq):
        """Test and count if exact duplicate."""
        seq_hash = zlib.adler32(bytearray(str(seq), "utf-8"))
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
        return self.match_count[int(index)]
