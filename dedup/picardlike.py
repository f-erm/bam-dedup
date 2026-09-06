#!/usr/bin/env python3
"""
pymarkdup.py  --  a faithful, dependency-light re-implementation of the core of
Picard MarkDuplicates (coordinate-sorted, Illumina paired-end, default options).

It reproduces Picard's decisions about WHICH records get the 0x400 (PCR/optical
duplicate) flag, using the exact same algorithm:

    1. Reduce every primary mapped read to a "ReadEnds": the 5' *unclipped*
       coordinate, orientation, library, and (optionally) the cluster location
       parsed from the read name for optical-duplicate detection.
    2. Pairs (both mates mapped) are collapsed into one ReadEnds keyed by the
       two 5' ends; single/half-mapped reads stay as fragment ReadEnds.
    3. Sort each list and walk it, grouping records that share
       (library, ref, coord, orientation [, mate ref, mate coord]).
    4. In each group keep the highest-scoring representative (default score =
       sum of base qualities >= 15) and flag the rest as duplicates.
    5. Rewrite the BAM setting/clearing the duplicate flag by index-in-file.

The only non-stdlib dependency is **pysam** (a C/htslib wrapper -- no JVM). This
removes the Java dependency, which is what you asked for.

SCOPE / faithful to Picard for:
    * coordinate-sorted input
    * DUPLICATE_SCORING_STRATEGY = SUM_OF_BASE_QUALITIES (Picard's default)
    * standard Illumina read names for optical detection
    * single or multiple read groups / libraries
    * REMOVE_DUPLICATES / REMOVE_SEQUENCING_DUPLICATES

barcode_tag (mark_duplicates/deduplicate keyword, --barcode-tag on the CLI) is
an intentional *departure* from Picard, not a port of it: it groups duplicate
candidates by (library, barcode_tag value, position, orientation [, mate]), so
records are only flagged as duplicates of one another if they also share the
same tag value -- e.g. a cell-barcode tag, so reads from different cells at
the same genomic position are never collapsed. Picard's own BARCODE_TAG does
UMI-aware marking via UmiUtil.getTopStrandNormalizedUmi, which in testing
silently finds zero duplicates whenever the tag value is constant across a
large group of reads (e.g. a cell barcode with no true per-molecule UMI
attached) -- not what most callers actually want, so it isn't reproduced here.

NOT (yet) handled -- these fall back to sane behavior but are not bit-identical:
    * queryname-sorted / query-grouped input (Picard marks unmapped mates &
      secondary/supplementary differently in that mode)
    * Picard's own BARCODE_TAG / UMI-aware marking semantics, DUPLEX_UMI
    * flow-based (FLOW_MODE) scoring
    * TOTAL_MAPPED_REFERENCE_LENGTH / RANDOM scoring strategies
    * DT / DI / DS tagging policies (easy to add; see notes)

Reference: picard/sam/markduplicates/MarkDuplicates.java and
htsjdk DuplicateScoringStrategy / OpticalDuplicateFinder.
"""

import argparse
import math
import sys

import pysam

# Optional Cython acceleration for the profiled hot spots. Falls back to pure
# Python (bit-identical) when the compiled module is not present.
try:
    from dedup import _fast
    _HAVE_FAST = True
except ImportError:  # pragma: no cover - pure-Python fallback
    _fast = None
    _HAVE_FAST = False

# ----------------------------------------------------------------------------
# Orientation constants -- mirror picard ReadEnds.java exactly.
# ----------------------------------------------------------------------------
F, R, FF, FR, RR, RF = 0, 1, 2, 3, 4, 5

# CIGAR ops that are clips (consume query but not reference for S; neither for H)
_SOFT_CLIP = 4
_HARD_CLIP = 5

SHORT_MAX = 32767       # Short.MAX_VALUE
SHORT_MIN = -32768      # Short.MIN_VALUE

DEFAULT_OPTICAL_PIXEL_DISTANCE = 100


def orientation_byte(read1_neg, read2_neg):
    """ReadEnds.getOrientationByte -- encode a pair's strand combination."""
    if read1_neg:
        return RR if read2_neg else RF
    else:
        return FR if read2_neg else FF


# ----------------------------------------------------------------------------
# ReadEnds: the compact per-read(-pair) record MarkDuplicates sorts and groups.
# ----------------------------------------------------------------------------
class ReadEnds:
    __slots__ = (
        "library_id", "orientation",
        "read1_ref", "read1_coord", "read2_ref", "read2_coord",
        "read1_index", "read2_index",
        "score",
        "read_group", "orientation_optical",
        "tile", "x", "y",
        "is_optical_duplicate",
        "barcode_value",
    )

    def __init__(self):
        self.library_id = -1
        self.orientation = -1
        self.read1_ref = -1
        self.read1_coord = -1
        self.read2_ref = -1
        self.read2_coord = -1
        self.read1_index = -1
        self.read2_index = -1
        self.score = 0
        self.read_group = -1
        self.orientation_optical = -1
        self.tile = -1
        self.x = -1
        self.y = -1
        self.is_optical_duplicate = False
        # "" when barcode_tag is not in use, so grouping is unaffected by default.
        self.barcode_value = ""

    @property
    def is_paired(self):
        # Matches ReadEnds.isPaired(): read2ReferenceIndex != -1
        return self.read2_ref != -1

    def has_location(self):
        return self.tile != -1

    def clone(self):
        c = ReadEnds()
        c.library_id = self.library_id
        c.orientation = self.orientation
        c.read1_ref = self.read1_ref
        c.read1_coord = self.read1_coord
        c.read2_ref = self.read2_ref
        c.read2_coord = self.read2_coord
        c.read1_index = self.read1_index
        c.read2_index = self.read2_index
        c.score = self.score
        c.read_group = self.read_group
        c.orientation_optical = self.orientation_optical
        c.tile = self.tile
        c.x = self.x
        c.y = self.y
        c.is_optical_duplicate = self.is_optical_duplicate
        c.barcode_value = self.barcode_value
        return c

    # Sort key == picard ReadEndsMDComparator.compare, extended with an
    # optional barcode_value (constant "" -- and thus a no-op -- unless
    # barcode_tag is passed to mark_duplicates/deduplicate).
    def sort_key(self):
        return (
            self.library_id,
            self.barcode_value,
            self.read1_ref,
            self.read1_coord,
            self.orientation,
            self.read2_ref,
            self.read2_coord,
            self.tile,
            self.x,
            self.y,
            self.read1_index,
            self.read2_index,
        )


# ----------------------------------------------------------------------------
# Scoring -- htsjdk DuplicateScoringStrategy, SUM_OF_BASE_QUALITIES.
# ----------------------------------------------------------------------------
def sum_of_base_qualities(read):
    """htsjdk getSumOfBaseQualities: sum of base qualities that are >= 15."""
    quals = read.query_qualities
    if quals is None:
        return 0
    if _HAVE_FAST:
        return _fast.sum_of_base_qualities(quals)
    s = 0
    for q in quals:
        if q >= 15:
            s += q
    return s


def compute_duplicate_score(read):
    """
    htsjdk computeDuplicateScore for SUM_OF_BASE_QUALITIES.
    Capped to fit a signed short; vendor-quality-failing reads are pushed to the
    lowest possible score so they are never chosen as the representative.
    """
    score = min(sum_of_base_qualities(read), SHORT_MAX - SHORT_MIN)  # cap
    # keep it in short range (the raw sum is already small in practice)
    if score > SHORT_MAX:
        score = SHORT_MAX
    if read.is_qcfail:
        # Java: (short) Math.max(score + Short.MIN_VALUE/2, Short.MIN_VALUE + 1)
        score = max(score + (SHORT_MIN // 2), SHORT_MIN + 1)
    return score


# ----------------------------------------------------------------------------
# 5' unclipped coordinate -- what Picard groups on.
# ----------------------------------------------------------------------------
def unclipped_5prime_coord(read):
    """
    Picard buildReadEnds:
        read1Coordinate = negStrand ? getUnclippedEnd() : getUnclippedStart()
    getUnclippedStart = alignmentStart  - leading  soft+hard clips
    getUnclippedEnd   = alignmentEnd    + trailing soft+hard clips
    Absolute value need not match Picard's 1-based value; only the *relative*
    equality between reads matters, and that is preserved.
    """
    cig = read.cigartuples
    if cig is None:
        return read.reference_start
    if not read.is_reverse:
        clip = 0
        for op, length in cig:
            if op == _SOFT_CLIP or op == _HARD_CLIP:
                clip += length
            else:
                break
        return read.reference_start - clip
    else:
        clip = 0
        for op, length in reversed(cig):
            if op == _SOFT_CLIP or op == _HARD_CLIP:
                clip += length
            else:
                break
        return read.reference_end + clip


# ----------------------------------------------------------------------------
# Optical-duplicate location parsing -- ReadNameParser default regex.
# ----------------------------------------------------------------------------
def parse_location(read_name):
    if _HAVE_FAST:
        return _fast.parse_location(read_name)
    return _parse_location_py(read_name)


def _parse_location_py(read_name):
    """
    Emulate ReadNameParser default parsing: split on ':', require the read name
    to have exactly 5 or 7 colon-separated fields, and take the LAST THREE as
    (tile, x, y). Returns (tile, x, y) or None if it doesn't match.
    Trailing non-digits in a field are ignored (rapidParseInt behavior).
    """
    fields = read_name.split(":")
    n = len(fields)
    if n != 5 and n != 7:
        return None
    try:
        tile = _rapid_parse_int(fields[-3])
        x = _rapid_parse_int(fields[-2])
        y = _rapid_parse_int(fields[-1])
    except ValueError:
        return None
    return tile, x, y


def _rapid_parse_int(s):
    """Parse leading (optionally negative) digits, stop at first non-digit."""
    i = 0
    n = len(s)
    neg = False
    if n > 0 and s[0] == "-":
        i = 1
        neg = True
    val = 0
    has = False
    while i < n and s[i].isdigit():
        val = val * 10 + (ord(s[i]) - 48)
        has = True
        i += 1
    if not has:
        raise ValueError(s)
    return -val if neg else val


# ----------------------------------------------------------------------------
# Optical duplicate finding -- OpticalDuplicateFinder (fast + union-find graph).
# Only affects optical *counts* / metrics / REMOVE_SEQUENCING_DUPLICATES, never
# the library duplicate flag itself.
# ----------------------------------------------------------------------------
def _close_enough(a, b, dist):
    return (a is not b and a.has_location() and b.has_location()
            and a.read_group == b.read_group
            and a.tile == b.tile
            and abs(a.x - b.x) <= dist
            and abs(a.y - b.y) <= dist)


def find_optical_duplicates(ends, keeper, dist):
    """Return a boolean list flagging which entries are optical duplicates."""
    length = len(ends)
    flags = [False] * length
    if length < 2:
        return flags

    actual_keeper = keeper if (keeper is not None and keeper.has_location()
                               and keeper in ends) else None

    if length >= (3 if actual_keeper is None else 4):
        return _optical_with_graph(ends, actual_keeper, flags, dist)
    return _optical_fast(ends, actual_keeper, flags, dist)


def _optical_fast(ends, keeper, flags, dist):
    length = len(ends)
    if keeper is not None:
        for i in range(length):
            flags[i] = _close_enough(keeper, ends[i], dist)
    for i in range(length):
        lhs = ends[i]
        if lhs is keeper:
            continue
        for j in range(i + 1, length):
            rhs = ends[j]
            if rhs is keeper:
                continue
            if flags[i] and flags[j]:
                continue
            if _close_enough(lhs, rhs, dist):
                idx = i if flags[j] else j
                flags[idx] = True
    return flags


def _optical_with_graph(ends, keeper, flags, dist):
    if _HAVE_FAST:
        keeper_index = -1
        if keeper is not None:
            for i, e in enumerate(ends):
                if e is keeper:
                    keeper_index = i
                    break
        xs = [e.x for e in ends]
        ys = [e.y for e in ends]
        tiles = [e.tile for e in ends]
        rgs = [e.read_group for e in ends]
        return _fast.optical_flags_graph(xs, ys, tiles, rgs, keeper_index, dist)

    # Union-find over reads within pixel distance (same read group + tile).
    parent = list(range(len(ends)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    tile_rg = {}
    keeper_index = -1
    for i, e in enumerate(ends):
        if e is keeper:
            keeper_index = i
        if e.has_location():
            key = ((e.read_group & 0xFFFF) << 16) + (e.tile & 0xFFFF)
            tile_rg.setdefault(key, []).append(i)

    for group in tile_rg.values():
        for a in range(len(group)):
            ia = group[a]
            for b in range(a + 1, len(group)):
                ib = group[b]
                if (abs(ends[ia].x - ends[ib].x) <= dist
                        and abs(ends[ia].y - ends[ib].y) <= dist):
                    union(ia, ib)

    cluster_rep = {}
    keeper_cluster = None
    if keeper_index >= 0:
        keeper_cluster = find(keeper_index)
        cluster_rep[keeper_cluster] = keeper_index

    for i in range(len(ends)):
        c = find(i)
        if c in cluster_rep and i != keeper_index:
            rep = ends[cluster_rep[c]]
            cur = ends[i]
            in_keeper_cluster = (keeper_index >= 0 and c == keeper_cluster)
            if (not in_keeper_cluster and
                    (cur.x < rep.x or (cur.x == rep.x and cur.y < rep.y))):
                flags[cluster_rep[c]] = True
                cluster_rep[c] = i
            else:
                flags[i] = True
        else:
            cluster_rep[c] = i
    return flags


def track_optical_duplicates(ends, keeper, dist):
    """
    AbstractMarkDuplicatesCommandLineProgram.trackOpticalDuplicates: partition by
    orientationForOpticalDuplicates when both FR and RF are present, then flag.
    Sets is_optical_duplicate on entries and returns the optical-duplicate count.
    """
    has_fr = any(e.orientation_optical == FR for e in ends)
    has_rf = any(e.orientation_optical == RF for e in ends)

    def _flag(sublist):
        f = find_optical_duplicates(sublist, keeper, dist)
        c = 0
        for i, is_opt in enumerate(f):
            if is_opt:
                sublist[i].is_optical_duplicate = True
                c += 1
        return c

    if has_fr and has_rf:
        fr = [e for e in ends if e.orientation_optical == FR]
        rf = [e for e in ends if e.orientation_optical == RF]
        return _flag(fr) + _flag(rf)
    return _flag(ends)


# ----------------------------------------------------------------------------
# Phase 1: build sorted read-end lists (single pass over coordinate-sorted BAM).
# ----------------------------------------------------------------------------
def build_read_ends(bam, read_name_regex_enabled, optical_dist, barcode_tag=None):
    header = bam.header
    # library-id: map library (LB) name -> small int; RG id -> library, RG ordinal
    rg_records = header.get("RG", [])
    rg_to_library = {}
    rg_to_ordinal = {}
    for ordinal, rg in enumerate(rg_records):
        rg_to_library[rg["ID"]] = rg.get("LB", "Unknown Library")
        rg_to_ordinal[rg["ID"]] = ordinal

    library_ids = {}
    next_library_id = [1]

    def library_id_for(read):
        try:
            rg = read.get_tag("RG")
        except KeyError:
            rg = None
        lib = rg_to_library.get(rg, "Unknown Library") if rg is not None else "Unknown Library"
        lid = library_ids.get(lib)
        if lid is None:
            lid = next_library_id[0]
            next_library_id[0] += 1
            library_ids[lib] = lid
        return lid, lib

    def make_end(read, index):
        e = ReadEnds()
        e.read1_ref = read.reference_id
        e.read1_coord = unclipped_5prime_coord(read)
        e.orientation = R if read.is_reverse else F
        e.read1_index = index
        e.score = compute_duplicate_score(read)
        if read.is_paired and not read.mate_is_unmapped:
            e.read2_ref = read.next_reference_id
        e.library_id, _ = library_id_for(read)
        if barcode_tag is not None:
            try:
                e.barcode_value = read.get_tag(barcode_tag)
            except KeyError:
                e.barcode_value = ""
        # optical location
        if read_name_regex_enabled:
            loc = parse_location(read.query_name)
            if loc is not None:
                e.tile, e.x, e.y = loc
                try:
                    rg = read.get_tag("RG")
                except KeyError:
                    rg = None
                e.read_group = rg_to_ordinal.get(rg, 0) if rg is not None else 0
        return e

    frag_list = []
    pair_list = []
    pending = {}   # (mate_ref, rg+readname) -> ReadEnds waiting for its mate

    index = 0
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped:
            # coordinate-sorted: trailing unmapped (ref==-1) reads carry no info
            if read.reference_id == -1:
                # they still need to be copied out later; just stop collecting.
                index += 1
                continue
            index += 1
            continue
        if read.is_secondary or read.is_supplementary:
            index += 1
            continue

        frag = make_end(read, index)
        frag_list.append(frag)

        if read.is_paired and not read.mate_is_unmapped:
            try:
                rg = read.get_tag("RG")
            except KeyError:
                rg = ""
            key = "{}{}".format(rg, read.query_name)
            paired = pending.pop((read.reference_id, key), None)
            if paired is None:
                # First mate seen: stash a clone keyed by THIS read's mate ref.
                paired = frag.clone()
                pending[(read.next_reference_id, key)] = paired
            else:
                _combine_mate(paired, frag, read, optical=True)
                paired.score += compute_duplicate_score(read)
                pair_list.append(paired)
        index += 1

    return frag_list, pair_list, library_ids


def _combine_mate(paired, frag, read, optical):
    """
    Merge the second mate (frag/read) into the stashed first mate (paired),
    reproducing MarkDuplicates.buildSortedReadEndLists lines 575-620:
    order the two ends so read1 <= read2, and set the pair orientation.
    `paired.orientation` at entry holds the FIRST mate's single orientation (F/R).
    """
    mates_ref = frag.read1_ref
    mates_coord = frag.read1_coord

    # orientationForOpticalDuplicates: always first-end then second-end strands.
    if read.is_read1:
        paired.orientation_optical = orientation_byte(read.is_reverse,
                                                      paired.orientation == R)
    else:
        paired.orientation_optical = orientation_byte(paired.orientation == R,
                                                      read.is_reverse)

    if (mates_ref > paired.read1_ref or
            (mates_ref == paired.read1_ref and mates_coord >= paired.read1_coord)):
        paired.read2_ref = mates_ref
        paired.read2_coord = mates_coord
        paired.read2_index = frag.read1_index
        paired.orientation = orientation_byte(paired.orientation == R,
                                              read.is_reverse)
        # Undefined RF at identical position -> force FR (see Picard comment).
        if (paired.read2_ref == paired.read1_ref and
                paired.read2_coord == paired.read1_coord and
                paired.orientation == RF):
            paired.orientation = FR
    else:
        paired.read2_ref = paired.read1_ref
        paired.read2_coord = paired.read1_coord
        paired.read2_index = paired.read1_index
        paired.read1_ref = mates_ref
        paired.read1_coord = mates_coord
        paired.read1_index = frag.read1_index
        paired.orientation = orientation_byte(read.is_reverse,
                                              paired.orientation == R)


# ----------------------------------------------------------------------------
# Phase 2: generate duplicate indexes from the sorted lists.
# ----------------------------------------------------------------------------
def _comparable(lhs, rhs, compare_read2):
    if lhs.library_id != rhs.library_id:
        return False
    if lhs.barcode_value != rhs.barcode_value:
        return False
    if not (lhs.read1_ref == rhs.read1_ref and
            lhs.read1_coord == rhs.read1_coord and
            lhs.orientation == rhs.orientation):
        return False
    if compare_read2:
        return (lhs.read2_ref == rhs.read2_ref and
                lhs.read2_coord == rhs.read2_coord)
    return True


def generate_duplicate_indexes(frag_list, pair_list, index_optical, optical_dist, rmlog):
    duplicate_indexes = set()
    optical_indexes = set()
    optical_cluster_count = 0

    logfile = open("w",rmlog)

    # ---- pairs ----
    pair_list.sort(key=ReadEnds.sort_key)
    for chunk in _chunks(pair_list, compare_read2=True):
        if len(chunk) > 1:
            logfile.write(str(len(chunk))+"\n")
            optical_cluster_count += _mark_pairs(
                chunk, duplicate_indexes, optical_indexes,
                index_optical, optical_dist)

    # ---- fragments ----
    frag_list.sort(key=ReadEnds.sort_key)
    first = None
    contains_pairs = False
    contains_frags = False
    current = []
    for nxt in frag_list:
        if first is not None and _comparable(first, nxt, compare_read2=False):
            current.append(nxt)
            contains_pairs = contains_pairs or nxt.is_paired
            contains_frags = contains_frags or (not nxt.is_paired)
        else:
            if len(current) > 1 and contains_frags:
                logfile.write(str(len(current))+"\n")
                _mark_fragments(current, contains_pairs, duplicate_indexes)
            current = [nxt]
            first = nxt
            contains_pairs = nxt.is_paired
            contains_frags = not nxt.is_paired
    if len(current) > 1 and contains_frags:
        logfile.write(str(len(current))+"\n")
        _mark_fragments(current, contains_pairs, duplicate_indexes)

    logfile.close()
    return duplicate_indexes, optical_indexes, optical_cluster_count


def _chunks(sorted_list, compare_read2):
    """Yield maximal runs of consecutive comparable ReadEnds."""
    first = None
    current = []
    for nxt in sorted_list:
        if first is not None and _comparable(first, nxt, compare_read2):
            current.append(nxt)
        else:
            if current:
                yield current
            current = [nxt]
            first = nxt
    if current:
        yield current


def _best(chunk):
    """Highest score wins; ties broken by sort order (strict > like Picard)."""
    best = None
    max_score = 0
    for e in chunk:
        if best is None or e.score > max_score:
            max_score = e.score
            best = e
    return best


def _mark_pairs(chunk, duplicate_indexes, optical_indexes,
                index_optical, optical_dist):
    best = _best(chunk)
    n_optical = 0
    # optical detection over the whole chunk (keeper = best)
    n_optical = track_optical_duplicates(chunk, best, optical_dist)

    for e in chunk:
        if e is best:
            continue
        duplicate_indexes.add(e.read1_index)
        if e.read2_index != e.read1_index:
            duplicate_indexes.add(e.read2_index)
        if e.is_optical_duplicate and index_optical:
            optical_indexes.add(e.read1_index)
            if e.read2_index != e.read1_index:
                optical_indexes.add(e.read2_index)
    return 1 if n_optical > 0 else 0


def _mark_fragments(chunk, contains_pairs, duplicate_indexes):
    if contains_pairs:
        # Any unpaired fragment sharing a start with a pair is a duplicate.
        for e in chunk:
            if not e.is_paired:
                duplicate_indexes.add(e.read1_index)
    else:
        best = _best(chunk)
        for e in chunk:
            if e is not best:
                duplicate_indexes.add(e.read1_index)


# ----------------------------------------------------------------------------
# Phase 3: rewrite BAM setting/clearing the duplicate flag by index-in-file.
# ----------------------------------------------------------------------------
def write_output(in_path, out_path, duplicate_indexes, optical_indexes,
                 remove_duplicates, remove_sequencing_duplicates):
    with pysam.AlignmentFile(in_path, "rb") as bam_in:
        with pysam.AlignmentFile(out_path, "wb", template=bam_in) as bam_out:
            index = 0
            n_dup = 0
            for read in bam_in.fetch(until_eof=True):
                is_dup = index in duplicate_indexes
                is_optical = index in optical_indexes
                read.is_duplicate = is_dup
                if is_dup:
                    n_dup += 1
                index += 1
                if remove_duplicates and read.is_duplicate:
                    continue
                if remove_sequencing_duplicates and is_optical:
                    continue
                bam_out.write(read)
    return n_dup


# ----------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------
def mark_duplicates(input_bam, output_bam,
                    remove_duplicates=False,
                    remove_sequencing_duplicates=False,
                    read_name_regex_enabled=True,
                    optical_pixel_distance=DEFAULT_OPTICAL_PIXEL_DISTANCE,
                    metrics_file=None,
                    barcode_tag=None,
                    rmlog=None):
    """
    barcode_tag: when set, reads/pairs are additionally grouped by the value of
    this tag before duplicate detection, so records that would otherwise look
    like duplicates (same library, position, orientation) are only flagged as
    duplicates of each other if they also share the same barcode_tag value.
    Reads missing the tag are grouped together under a single "no barcode"
    bucket. This is a deliberate, well-defined design -- not a port of Picard's
    BARCODE_TAG (which does UMI-aware duplicate marking via
    UmiUtil.getTopStrandNormalizedUmi and, in testing, silently finds zero
    duplicates when the tag value is constant across a large group of reads,
    e.g. a single-cell barcode with no true per-molecule UMI).
    """
    with pysam.AlignmentFile(input_bam, "rb") as bam:
        so = bam.header.get("HD", {}).get("SO", "unknown")
        if so != "coordinate":
            sys.stderr.write(
                "WARNING: input SO='{}'. This implementation is faithful only "
                "for coordinate-sorted input.\n".format(so))
        frag_list, pair_list, _ = build_read_ends(
            bam, read_name_regex_enabled, optical_pixel_distance, barcode_tag)

    index_optical = remove_sequencing_duplicates or (metrics_file is not None)
    dup_idx, opt_idx, opt_clusters = generate_duplicate_indexes(
        frag_list, pair_list, index_optical, optical_pixel_distance, rmlog)

    n_dup = write_output(input_bam, output_bam, dup_idx, opt_idx,
                         remove_duplicates, remove_sequencing_duplicates)

    sys.stderr.write(
        "Marked {} records as duplicates ({} optical-duplicate clusters).\n"
        .format(n_dup, opt_clusters))
    return dup_idx, opt_idx


def deduplicate(input_bam, output_bam, **kwargs):
    """
    Convenience wrapper: produce a *deduplicated* BAM (duplicate records removed,
    not just flagged). Equivalent to mark_duplicates(..., remove_duplicates=True).
    All other keyword arguments are forwarded to :func:`mark_duplicates`.
    """
    kwargs["remove_duplicates"] = True
    return mark_duplicates(input_bam, output_bam, **kwargs)


def main(argv=None):
    p = argparse.ArgumentParser(description="Minimal Picard MarkDuplicates clone.")
    p.add_argument("-i", "--input", required=True, help="coordinate-sorted BAM")
    p.add_argument("-o", "--output", required=True, help="output BAM")
    p.add_argument("--remove-duplicates", action="store_true")
    p.add_argument("--remove-sequencing-duplicates", action="store_true")
    p.add_argument("--no-optical", action="store_true",
                   help="disable optical-duplicate detection (READ_NAME_REGEX=null)")
    p.add_argument("--optical-pixel-distance", type=int,
                   default=DEFAULT_OPTICAL_PIXEL_DISTANCE)
    p.add_argument("--metrics-file", default=None)
    p.add_argument("--barcode-tag", default=None,
                   help="group duplicate candidates by this tag's value in addition to "
                        "position/orientation (e.g. a cell barcode tag); reads sharing a "
                        "position but carrying different tag values are never collapsed")
    args = p.parse_args(argv)

    mark_duplicates(
        args.input, args.output,
        remove_duplicates=args.remove_duplicates,
        remove_sequencing_duplicates=args.remove_sequencing_duplicates,
        read_name_regex_enabled=not args.no_optical,
        optical_pixel_distance=args.optical_pixel_distance,
        metrics_file=args.metrics_file,
        barcode_tag=args.barcode_tag,
    )


if __name__ == "__main__":
    main()
