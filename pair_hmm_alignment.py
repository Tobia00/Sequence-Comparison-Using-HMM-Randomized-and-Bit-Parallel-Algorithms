import math
import time

NEG_INF = -10**18


def read_fasta_sequence(file_path):
    """Read a FASTA file and return the sequence as one uppercase string."""
    sequence = ""
    with open(file_path, "r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()
            if line and not line.startswith(">"):
                sequence += line.upper()
    return sequence


def safe_log(value):
    if value <= 0:
        return NEG_INF
    return math.log(value)


def pair_hmm_viterbi(seq1, seq2):
    """Pair-HMM Viterbi alignment using Match, Insertion, and Deletion states."""
    n, m = len(seq1), len(seq2)

    delta = 0.05
    epsilon = 0.80

    transitions = {
        "M": {"M": 1 - 2 * delta, "I": delta, "D": delta},
        "I": {"M": 1 - epsilon, "I": epsilon, "D": 0},
        "D": {"M": 1 - epsilon, "I": 0, "D": epsilon},
    }

    match_emit = 0.90
    mismatch_emit = 0.10 / 3
    gap_emit = 0.20

    M = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    I = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    D = [[NEG_INF] * (m + 1) for _ in range(n + 1)]

    back_M = [[None] * (m + 1) for _ in range(n + 1)]
    back_I = [[None] * (m + 1) for _ in range(n + 1)]
    back_D = [[None] * (m + 1) for _ in range(n + 1)]

    M[0][0] = 0

    for i in range(1, n + 1):
        D[i][0] = (D[i - 1][0] if i > 1 else M[i - 1][0]) + safe_log(gap_emit)
        back_D[i][0] = "D" if i > 1 else "M"

    for j in range(1, m + 1):
        I[0][j] = (I[0][j - 1] if j > 1 else M[0][j - 1]) + safe_log(gap_emit)
        back_I[0][j] = "I" if j > 1 else "M"

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            emit = safe_log(match_emit if seq1[i - 1] == seq2[j - 1] else mismatch_emit)

            candidates_m = [
                (M[i - 1][j - 1] + safe_log(transitions["M"]["M"]), "M"),
                (I[i - 1][j - 1] + safe_log(transitions["I"]["M"]), "I"),
                (D[i - 1][j - 1] + safe_log(transitions["D"]["M"]), "D"),
            ]
            M[i][j], back_M[i][j] = max(candidates_m, key=lambda x: x[0])
            M[i][j] += emit

            candidates_i = [
                (M[i][j - 1] + safe_log(transitions["M"]["I"]), "M"),
                (I[i][j - 1] + safe_log(transitions["I"]["I"]), "I"),
                (D[i][j - 1] + safe_log(transitions["D"]["I"]), "D"),
            ]
            I[i][j], back_I[i][j] = max(candidates_i, key=lambda x: x[0])
            I[i][j] += safe_log(gap_emit)

            candidates_d = [
                (M[i - 1][j] + safe_log(transitions["M"]["D"]), "M"),
                (I[i - 1][j] + safe_log(transitions["I"]["D"]), "I"),
                (D[i - 1][j] + safe_log(transitions["D"]["D"]), "D"),
            ]
            D[i][j], back_D[i][j] = max(candidates_d, key=lambda x: x[0])
            D[i][j] += safe_log(gap_emit)

    final_score, state = max([(M[n][m], "M"), (I[n][m], "I"), (D[n][m], "D")], key=lambda x: x[0])

    aligned1, aligned2, hidden_path = [], [], []
    i, j = n, m

    while i > 0 or j > 0:
        hidden_path.append(state)

        if state == "M":
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            prev = back_M[i][j]
            i -= 1
            j -= 1
        elif state == "I":
            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            prev = back_I[i][j]
            j -= 1
        else:
            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            prev = back_D[i][j]
            i -= 1

        state = prev if prev is not None else "M"

    return "".join(reversed(aligned1)), "".join(reversed(aligned2)), list(reversed(hidden_path)), final_score


def calculate_metrics(aligned_seq1, aligned_seq2):
    matches = mismatches = insertions = deletions = 0

    for a, b in zip(aligned_seq1, aligned_seq2):
        if a == "-":
            insertions += 1
        elif b == "-":
            deletions += 1
        elif a == b:
            matches += 1
        else:
            mismatches += 1

    alignment_length = len(aligned_seq1)
    identity = (matches / alignment_length) * 100 if alignment_length else 0
    return matches, mismatches, insertions, deletions, identity, alignment_length


def print_alignment(aligned_seq1, aligned_seq2, width=80):
    for i in range(0, len(aligned_seq1), width):
        s1 = aligned_seq1[i:i + width]
        s2 = aligned_seq2[i:i + width]
        mid = "".join("|" if a == b else " " for a, b in zip(s1, s2))
        print(s1)
        print(mid)
        print(s2)
        print()


def demo(seq1, seq2):
    start = time.time()
    aligned1, aligned2, hidden_path, score = pair_hmm_viterbi(seq1, seq2)
    matches, mismatches, insertions, deletions, identity, length = calculate_metrics(aligned1, aligned2)

    print("\nPAIR-HMM VITERBI RESULT")
    print("Alignment length:", length)
    print("Matches:", matches, "| mismatches:", mismatches, "| insertions:", insertions, "| deletions:", deletions)
    print("Identity:", round(identity, 2), "%")
    print("Best log-score:", round(score, 4))
    print("Runtime:", round(time.time() - start, 6), "seconds")
    print_alignment(aligned1, aligned2, width=80)
    print("Hidden path preview:", " ".join(hidden_path[:50]))
