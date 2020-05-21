# -*- coding: utf-8 -*-
"""Protein sequence checking and sanitization."""

# standard library imports
import zlib

# third-party imports
from Bio.Data import IUPACData

# global constants
AMBIGUOUS = "X"
STOP = "*"
DASH = "-"
START_CHARS = ("M",)
ALPHABET = IUPACData.protein_letters + AMBIGUOUS


class Sanitizer:

    """Count and clean up problems with protein sequence.

    Problems recognized are:
          alphabet:  if not in IUPAC set, changed to 'X'
            dashes:    optional, removed if remove_dashes=True
         ambiguous:

    """

    def __init__(self, remove_stops=True, dashes_ok=False, start_chars=None):
        """Initialize counters."""
        self.remove_stops = remove_stops
        if start_chars is None:
            self.starts = START_CHARS
        else:
            self.starts = tuple(start_chars)
        if dashes_ok:
            self.alphabet = ALPHABET + DASH
            self.remove_dashes = False
        else:
            self.alphabet = ALPHABET
            self.remove_dashes = True
        self.improper_starts = 0
        self.seqs_sanitized = 0
        self.seqs_out = 0
        self.resid_removed = 0
        self.resid_fixed = 0
        self.resid_out = 0
        self.ambiguous = 0
        self.stops = 0

    def char_remover(self, seq, character):
        """Remove positions with a given character.

        :param seq: mutable sequence
        :return: sequence with characters removed
        """
        removals = [i for i, j in enumerate(seq) if j == character]
        self.resid_removed += len(removals)
        for k, pos in enumerate(removals):
            seq.pop(pos - k)
        return seq

    def fix_alphabet(self, seq):
        """Replace everything out of alphabet with AMBIGUOUS.

        :param seq: mutable sequence, upper-cased
        :return: fixed sequence
        """
        fix_positions = [
            pos for pos, char in enumerate(seq) if char not in self.alphabet
        ]
        self.resid_fixed += len(fix_positions)
        for pos in fix_positions:
            seq.__setitem__(pos, AMBIGUOUS)
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
        self.resid_removed += in_len - len(seq)
        return seq

    def stops_properly(self, seq):
        """Remove a single stop on the end, if it exists."""
        stops = seq[-1] == STOP
        if stops and self.remove_stops:
            seq.pop()
        return seq, stops

    def sanitize(self, seq):
        """Sanitize potential problems with sequence.

        Remove dashes, change non-IUPAC characters to
        ambiguous, and remove ambiguous characters on ends.
        :param seq: mutable sequence
        :return: sanitized sequence
        """
        self.seqs_sanitized += 1
        if len(seq) == 0:
            raise ValueError("zero-length sequence")
        seq, has_stop = self.stops_properly(seq)
        if self.remove_dashes:
            seq = self.char_remover(seq, DASH)
        if len(seq) == 0:
            raise ValueError("zero-length sequence after stop removed")
        seq = self.fix_alphabet(seq)
        seq = self.remove_char_on_ends(seq, AMBIGUOUS)
        if len(seq) == 0:
            raise ValueError("zero-length sequence after ends trimmed")
        self.resid_out += len(seq)
        self.stops += int(has_stop)
        n_ambig = self.count_ambiguous(seq)
        self.ambiguous += n_ambig
        bad_start = self.starts_improperly(seq)
        self.improper_starts += int(bad_start)
        self.seqs_out += 1
        return seq, has_stop, bad_start, n_ambig

    def count_ambiguous(self, seq):
        """Count ambiguous residues.

        :param seq: sequence
        :return: Number of ambiguous residues
        """
        return sum([i == AMBIGUOUS for i in seq])

    def starts_improperly(self, seq):
        """Return True if start char not in list."""
        return seq[0] not in self.starts

    def file_stats(self):
        """Return a dictionary of file stats."""
        return {
            "idx": 0,
            "seqs.n": self.seqs_out,
            "seqs.rmv": self.seqs_sanitized - self.seqs_out,
            "seqs.stp": self.stops,
            "seqs.nostrt": self.improper_starts,
            "res.n": self.resid_out,
            "res.fix": self.resid_fixed,
            "res.rmv": self.resid_removed,
            "res.ambg": self.ambiguous,
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
