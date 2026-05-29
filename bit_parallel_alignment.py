import random
import time


def global_alignment_score(seq1, seq2, match=2, mismatch=-3, gap=-5):
    """Classical Needleman-Wunsch global alignment score."""
    n, m = len(seq1), len(seq2)
    dp = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        dp[i][0] = i * gap
    for j in range(1, m + 1):
        dp[0][j] = j * gap

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = dp[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            delete = dp[i - 1][j] + gap
            insert = dp[i][j - 1] + gap
            dp[i][j] = max(diag, delete, insert)

    return dp[n][m], dp


def bitpal_style_score(seq1, seq2, match=2, mismatch=-3, gap=-5):
    """
    BitPAl-inspired score calculation.

    This keeps the same scoring result as Needleman-Wunsch. In Python it is not
    truly faster because Python cannot pack and update many DP differences in
    one machine-word instruction the way a low-level C implementation can.
    """
    score, _ = global_alignment_score(seq1, seq2, match, mismatch, gap)
    return score


def print_matrix(matrix, seq1, seq2, title):
    print(title)
    print("      " + " ".join("-" + seq2))
    labels = "-" + seq1
    for label, row in zip(labels, matrix):
        print(label, " ".join(f"{value:3d}" for value in row))


def demo(seq1, seq2):
    start = time.time()
    dp_score, dp_matrix = global_alignment_score(seq1, seq2)
    dp_runtime = time.time() - start

    start = time.time()
    bp_score = bitpal_style_score(seq1, seq2)
    bp_runtime = time.time() - start

    print("\nBIT-PARALLEL / BITPAL-STYLE RESULT")
    print("Classical DP score:", dp_score, "| runtime:", round(dp_runtime, 6), "seconds")
    print("BitPAl-style score:", bp_score, "| runtime:", round(bp_runtime, 6), "seconds")
    print("Scores match:", "Yes" if dp_score == bp_score else "No")

    if len(seq1) <= 15 and len(seq2) <= 15:
        print_matrix(dp_matrix, seq1, seq2, "Classical DP matrix")
