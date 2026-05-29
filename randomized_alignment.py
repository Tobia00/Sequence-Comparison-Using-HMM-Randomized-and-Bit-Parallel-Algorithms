import random
import time


def read_fasta(file_path):
    """Read a FASTA file and return the sequence as one uppercase string."""
    sequence = ""
    with open(file_path, "r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()
            if line and not line.startswith(">"):
                sequence += line.upper()
    return sequence


def calculate_score(seq1, seq2, match_score=1, mismatch_score=-1):
    """Calculate a simple ungapped alignment score."""
    return sum(match_score if a == b else mismatch_score for a, b in zip(seq1, seq2))


def calculate_alignment_score(seq1, seq2, match_score=1, mismatch_score=-1, gap_score=-2):
    """Calculate score and summary metrics for two aligned sequences."""
    score = matches = mismatches = gaps = 0

    for a, b in zip(seq1, seq2):
        if a == "-" or b == "-":
            score += gap_score
            gaps += 1
        elif a == b:
            score += match_score
            matches += 1
        else:
            score += mismatch_score
            mismatches += 1

    return score, matches, mismatches, gaps


def monte_carlo_alignment(seq1, seq2, k=1000, window_size=100):
    """Monte Carlo random window alignment."""
    if window_size > len(seq1) or window_size > len(seq2):
        window_size = min(len(seq1), len(seq2))

    best_score = float("-inf")
    best_alignment = None
    best_start_seq1 = best_start_seq2 = 0

    for _ in range(k):
        start1 = random.randint(0, len(seq1) - window_size)
        start2 = random.randint(0, len(seq2) - window_size)

        window1 = seq1[start1:start1 + window_size]
        window2 = seq2[start2:start2 + window_size]
        score = calculate_score(window1, window2)

        if score > best_score:
            best_score = score
            best_alignment = (window1, window2)
            best_start_seq1 = start1
            best_start_seq2 = start2

    return best_alignment, best_score, best_start_seq1, best_start_seq2


def las_vegas_alignment(seq1, seq2, max_trials=1000, window_size=100, max_mismatches=30):
    """Las Vegas random window alignment that returns a verified alignment if found."""
    if window_size > len(seq1) or window_size > len(seq2):
        window_size = min(len(seq1), len(seq2))

    best_alignment = ("", "")
    best_score = float("-inf")
    best_start_seq1 = best_start_seq2 = 0

    for trial in range(1, max_trials + 1):
        start1 = random.randint(0, len(seq1) - window_size)
        start2 = random.randint(0, len(seq2) - window_size)

        window1 = seq1[start1:start1 + window_size]
        window2 = seq2[start2:start2 + window_size]

        score = calculate_score(window1, window2)
        mismatches = sum(1 for a, b in zip(window1, window2) if a != b)

        if score > best_score:
            best_score = score
            best_alignment = (window1, window2)
            best_start_seq1 = start1
            best_start_seq2 = start2

        if mismatches <= max_mismatches:
            return window1, window2, score, mismatches, start1, start2, trial, True

    return best_alignment[0], best_alignment[1], best_score, None, best_start_seq1, best_start_seq2, max_trials, False


def insert_random_gaps(shorter_seq, target_length):
    """Insert random gaps into the shorter sequence until it reaches target length."""
    aligned = list(shorter_seq)
    while len(aligned) < target_length:
        pos = random.randint(0, len(aligned))
        aligned.insert(pos, "-")
    return "".join(aligned)


def make_same_length_with_random_gaps(seq1, seq2):
    """Randomly insert gaps so both sequences have the same length."""
    if len(seq1) < len(seq2):
        return insert_random_gaps(seq1, len(seq2)), seq2
    if len(seq2) < len(seq1):
        return seq1, insert_random_gaps(seq2, len(seq1))
    return seq1, seq2


def monte_carlo_full_alignment_with_gaps(seq1, seq2, k=200):
    """Monte Carlo full-sequence alignment with random gap insertion."""
    best_score = float("-inf")
    best_alignment = None
    best_stats = None

    for trial in range(1, k + 1):
        aligned_seq1, aligned_seq2 = make_same_length_with_random_gaps(seq1, seq2)
        score, matches, mismatches, gaps = calculate_alignment_score(aligned_seq1, aligned_seq2)

        if score > best_score:
            best_score = score
            best_alignment = (aligned_seq1, aligned_seq2)
            best_stats = (matches, mismatches, gaps, trial)

    return best_alignment, best_score, best_stats


def las_vegas_full_alignment_with_gaps(seq1, seq2, max_trials=200, max_mismatches=900):
    """Las Vegas full-sequence alignment with mismatch validation."""
    best_score = float("-inf")
    best_alignment = None
    best_stats = None

    for trial in range(1, max_trials + 1):
        aligned_seq1, aligned_seq2 = make_same_length_with_random_gaps(seq1, seq2)
        score, matches, mismatches, gaps = calculate_alignment_score(aligned_seq1, aligned_seq2)

        if score > best_score:
            best_score = score
            best_alignment = (aligned_seq1, aligned_seq2)
            best_stats = (matches, mismatches, gaps, trial)

        if mismatches <= max_mismatches:
            return aligned_seq1, aligned_seq2, score, matches, mismatches, gaps, trial, True

    matches, mismatches, gaps, best_trial = best_stats
    return best_alignment[0], best_alignment[1], best_score, matches, mismatches, gaps, best_trial, False


def demo(seq1, seq2):
    """Run a short randomized-algorithm demonstration."""
    start = time.time()
    alignment, score, pos1, pos2 = monte_carlo_alignment(seq1, seq2, k=100, window_size=20)
    print("\nMONTE CARLO WINDOW RESULT")
    print("Best score:", score, "| positions:", pos1, pos2)
    print(alignment[0])
    print(alignment[1])

    w1, w2, score, mismatches, pos1, pos2, trials, found = las_vegas_alignment(
        seq1, seq2, max_trials=100, window_size=20, max_mismatches=8
    )
    print("\nLAS VEGAS WINDOW RESULT")
    print("Found:", found, "| score:", score, "| mismatches:", mismatches, "| trials:", trials)
    print(w1)
    print(w2)
    print("Randomized demo runtime:", round(time.time() - start, 6), "seconds")
